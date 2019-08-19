import os
from stompy import utils
import pandas as pd
import logging as log
import calc_settling
import numpy as np

##
try:
    cwd=os.path.dirname(__file__)
except NameError:
    cwd="." # assume interactive run is from this directory.

if 0:
    # earlier dataset, manually extracted to csv
    effluent_df=pd.read_csv(os.path.join(cwd,"2019-04-08_Datafile_Effluent.csv"))
    effluent_df['pathway']='effluent'
    storm_df=pd.read_csv(os.path.join(cwd,"2019-04-08_Data_Stormwater.csv"))
    storm_df['pathway']='stormwater'
    dfs=[effluent_df,storm_df]
else:
    xlsx_fn=os.path.join(cwd,'2019-08-19_DataExport_ForRH.xlsx')
    # also Manta, Sediment, Fish

    # code needs Category_Final {Fiber, Film, Fiber Bundle, Fragment, Sphere, Foam
    #   that all seems fine.
    # PlasticType_Final -- rename.
    # Length.mm, Width.mm - these are fine.
    def load_sheet(sheet_name):
        csv_fn=os.path.join(cwd,sheet_name+".csv")
        if utils.is_stale(csv_fn,[xlsx_fn,__file__]):
            log.info("Reading from xlsx")
            df=pd.read_excel(xlsx_fn,sheet_name=sheet_name)
            df.to_csv(csv_fn)
        else:
            log.info("Reading from csv")
            df=pd.read_csv(csv_fn)
        return df
            
    log.info("Reading effluent sheet")
    effluent_df=load_sheet("Effluent")
    log.info("Reading stormwater sheet")
    storm_df=load_sheet("Stormwater")
    log.info("Reading sediment sheet")
    sediment_df=load_sheet("Sediment")
    log.info("Reading manta sheet")
    manta_df=load_sheet("Manta")
    log.info("Reading fish sheet")
    fish_df=load_sheet("Fish")

    dfs=[effluent_df,storm_df,sediment_df,manta_df,fish_df]

    matrix_to_pathway={'eff':'effluent',
                       'sw':'stormwater',
                       'sed':'sediment',
                       'manta':'manta',
                       'fish':'fish'}
    # copy these to simplify processing
    for df in dfs:
        df['PlasticType_Final']=df['PlasticType_TrueFinal']
        rubber_idxs=[]
        for idx,row in df.iterrows():
            if row['Rubbery']=='Yes':
                rubber_types= ['Not Identified',
                               'Unknown Potentially Rubber',
                               'Anthropogenic (synthetic)',
                               'Unknown','Anthropogenic (unknown base)']
                if row['PlasticType_Final'] in rubber_types:
                    rubber_idxs.append(idx)
        rubber_before=(df.PlasticType_Final=='Rubber').sum()
        df.loc[rubber_idxs,'PlasticType_Final']='Rubber'
        rubber_now=(df.PlasticType_Final=='Rubber').sum()
                    
        log.info(f"{len(rubber_idxs)} particles assumed to be Rubber ({rubber_before}->{rubber_now})")
                                                
        df['pathway']=[matrix_to_pathway[m] for m in df['MatrixID'].values]

        if ('SampleDate' in df) and np.issubdtype(df['SampleDate'].dtype,np.integer):
            # convert to proper dates:
            as_date=np.datetime64("2017-08-21") + (df['SampleDate']-42968)*np.timedelta64(1,'D')
            df['SampleDate']=as_date

log.info("calculating settling velocities")        
for df in dfs:
    w_s=[calc_settling.record_to_ws(row)
         for i,row in df.iterrows()]
    w_s=np.array(w_s)
    df['w_s']=w_s

combined=pd.concat(dfs,sort=False)

