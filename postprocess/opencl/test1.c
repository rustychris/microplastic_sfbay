#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// RH
#define CL_TARGET_OPENCL_VERSION 220

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
 
#define MAX_SOURCE_SIZE (0x100000)

float time_diff(struct timespec start, struct timespec end);

 
int main(void) {
    // Create the two input vectors
    struct timespec time1, time2;
    int platform_index,repeat;
    int i;
    const int LIST_SIZE = 50*1024000;
    float *A = (float*)malloc(sizeof(float)*LIST_SIZE);
    float *B = (float*)malloc(sizeof(float)*LIST_SIZE);
    float *C = (float*)malloc(sizeof(float)*LIST_SIZE);

    for(i = 0; i < LIST_SIZE; i++) {
        A[i] = 1+i;
        B[i] = 1+LIST_SIZE - i;
    }
 
    // Load the kernel source code into the array source_str
    FILE *fp;
    char *source_str;
    size_t source_size;
 
    //fp = fopen("vector_add_kernel.cl", "r");
    fp = fopen("vector_log_kernel.cl", "r");
    if (!fp) {
        fprintf(stderr, "Failed to load kernel.\n");
        exit(1);
    }
    source_str = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp);
    fclose( fp );
 
    // Get platform and device information
    // cl_platform_id platform_id = NULL;
    cl_platform_id platform_ids[4];
    
    cl_device_id device_id = NULL;   
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret = clGetPlatformIDs(4, platform_ids, &ret_num_platforms);
    printf("%d platforms\n",ret_num_platforms);
    // CL_DEVICE_TYPE_DEFAULT.
    // once the debs were installed, then this works with platform_ids[0]
    // and GPU

    for( platform_index=0;platform_index<ret_num_platforms;platform_index++) {
      printf("Trying platform %d of %d\n",platform_index,ret_num_platforms);
      
      ret = clGetDeviceIDs( platform_ids[platform_index], CL_DEVICE_TYPE_DEFAULT, 1, 
                            &device_id, &ret_num_devices);
      
      if(ret!=CL_SUCCESS) { // failing for CL_DEVICE_TYPE_GPU
        printf("Error: Failed to query platforms! (%d)\n", ret);
        return EXIT_FAILURE;      
      }
      
      // Create an OpenCL context
      cl_context context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);
      
      // Create a command queue
      cl_command_queue command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
      
      // Create memory buffers on the device for each vector 
      cl_mem a_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY, 
                                        LIST_SIZE * sizeof(float), NULL, &ret);
      cl_mem b_mem_obj = clCreateBuffer(context, CL_MEM_READ_ONLY,
                                        LIST_SIZE * sizeof(float), NULL, &ret);
      cl_mem c_mem_obj = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 
                                        LIST_SIZE * sizeof(float), NULL, &ret);
      
      // Copy the lists A and B to their respective memory buffers
      ret = clEnqueueWriteBuffer(command_queue, a_mem_obj, CL_TRUE, 0,
                                 LIST_SIZE * sizeof(float), A, 0, NULL, NULL);
      ret = clEnqueueWriteBuffer(command_queue, b_mem_obj, CL_TRUE, 0, 
                                 LIST_SIZE * sizeof(float), B, 0, NULL, NULL);
      
      // Create a program from the kernel source
      cl_program program = clCreateProgramWithSource(context, 1, 
                                                     (const char **)&source_str, (const size_t *)&source_size, &ret);
      
      // Build the program
      ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
      
      // Create the OpenCL kernel
      cl_kernel kernel = clCreateKernel(program, "vector_add", &ret);
      
      // Set the arguments of the kernel
      ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&a_mem_obj);
      ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&b_mem_obj);
      ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&c_mem_obj);
      
      // Measure time just around the actual computation. But check back in case
      // it's nonblocking.
      // Careful, though. This counts time across cores!
      //clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
      clock_gettime(CLOCK_REALTIME, &time1);
      
      // Execute the OpenCL kernel on the list
      size_t global_item_size = LIST_SIZE; // Process the entire lists
      size_t local_item_size = 64; // Divide work items into groups of 64
      for(repeat=0;repeat<50;repeat++) {
        ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, 
                                     &global_item_size, &local_item_size, 0, NULL, NULL);
      }

      // This comes back *very* fast, before any work is done.
      // and no immediate option to make the kernel blocking.
      
      // Read the memory buffer C on the device to the local variable C
      ret = clEnqueueReadBuffer(command_queue, c_mem_obj, CL_TRUE, 0, 
                                LIST_SIZE * sizeof(float), // LIST_SIZE * sizeof(int),
                                C, 0, NULL, NULL);
      
      clock_gettime(CLOCK_REALTIME, &time2);
      printf("Computation and fetch: %.3fs\n",time_diff(time1,time2));
      
      // Display the result to the screen
      // LIST_SIZE
      for(i = 0; i < 10; i++)
        printf("%f + %f = %f\n", A[i], B[i], C[i]);
      
      // Clean up
      ret = clFlush(command_queue);
      ret = clFinish(command_queue);
      ret = clReleaseKernel(kernel);
      ret = clReleaseProgram(program);
      ret = clReleaseMemObject(a_mem_obj);
      ret = clReleaseMemObject(b_mem_obj);
      ret = clReleaseMemObject(c_mem_obj);
      ret = clReleaseCommandQueue(command_queue);
      ret = clReleaseContext(context);
    }
    free(A);
    free(B);
    free(C);
    return 0;
}

float time_diff(struct timespec start, struct timespec end)
{
    float seconds;
    struct timespec temp;
    
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    seconds=temp.tv_sec + temp.tv_nsec/1e9;
    return seconds;
}

// multiplying floats: 1.6x speedup.
// multiplying floats by exp( log+log): about the same
// This was still not a great comparison.
// Making the task a bit heavier for each data point, nearing reality,
// we get a 3x or 4x speedup.
