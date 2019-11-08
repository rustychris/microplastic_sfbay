
def make_load_table_oldrun(con,load_name,load_fn,msg=log.debug):
    """
    bring concentrations into db as a temporary, in-memory
    table.
    This version assumes the old runs, where only a small
    number of stormwater sources are included.
    
    the new table is in the schema load, with table name = load_name
    """
    load_ds=xr.open_dataset(load_fn)
    
    curs=con.cursor()
    curs.execute("select name from ptm_group group by name")
    group_names=curs.fetchall()

    stormwater_scale=1/0.33
    wastewater_scale=1/0.70

    behavior_to_ws={'down50000':0.05,
                    'down5000':0.005,
                    'down500':0.0005,
                    'none':0.0,
                    'up500':-0.0005,
                    'up5000':-0.005,
                    'up50000':-0.05}

    # almost everything is stormwater, so just create a default
    # dicut and explicitly name the non-stormwater.
    source_map=defaultdict(lambda:'stormwater')
    source_map['cccsd']='CCCSD'
    source_map['sunnyvale']='SUNN'
    source_map['fs']='FSSD'
    source_map['palo_alto']='PA'
    source_map['san_jose']='SJ'
    source_map['src000']='EBDA'
    source_map['src001']='EBMUD'
    source_map['src002']='SFPUC'

    # these we don't actually use
    source_map['petaluma']='SKIP'

    # These shouldn't be used, but including just to be sure
    # that if they somehow show up, they won't contaminate
    # stormwater.
    source_map['SacRiver']='DELTA'
    source_map['SJRiver']='DELTA'

    # tuples of group_name, particles/m3.
    # can omit rows with 0 concentration.
    conc_rows=[] 

    for group_name, in group_names:
        m=re.match(r'(.*)_(down\d+|up\d+|none)(_rel.*)?',group_name)
        source=m.group(1)
        behavior=m.group(2)
        rel_time=m.group(3) # may be missing
        w_s=behavior_to_ws[behavior]

        source_name=source_map[source]
        if source_name in ['DELTA','SKIP']:
            conc=0.0
        else:
            conc=load_ds.conc.sel(source=source_name,w_s=w_s).item()

        if source_name=='stormwater':
            conc*=stormwater_scale
        else:
            conc*=wastewater_scale

        # load netcdf is in particles/l, but we want to set the calculation
        # up as particle/m3. Updated 2019-11-17
        conc*=1000
        if conc>0.0:
            conc_rows.append( (group_name,conc) )
        msg(f"{source:15s}  {behavior:9s} => {source_name:15s} {conc:.4f}")

    # --- load that into a table    
    try:
        curs.execute("DETACH load")
    except sql.OperationalError:
        print("memory table wasn't there.  no worries.")
    curs.execute("ATTACH ':memory:' as load")

    curs.execute(f"""CREATE TABLE load.{load_name} (
       group_name text,
       part_per_m3 double);""")
    curs.executemany(f"""INSERT into load.{load_name} (group_name,part_per_m3) VALUES (?,?)""",
                    conc_rows)
    con.commit()
    load_ds.close()
