__kernel void vector_add(__global const float *A, __global const float *B, __global float *C) {
  int loops=10;
  
    // Get the index of the current element to be processed
    int i = get_global_id(0);

    for(;loops>0;loops--){
    // Do the operation
    //C[i] = exp( log(A[i]) + log(B[i]));
    //C[i] = native_exp( native_log(A[i]) + native_log(B[i]));
    //C[i] = half_exp( half_log(A[i]) + half_log(B[i]));
      C[i] = native_exp( 5*native_log(A[i]) + 0.1*native_log(B[i]) +  native_log(A[i]) + native_log(B[i]));
    }    
}
