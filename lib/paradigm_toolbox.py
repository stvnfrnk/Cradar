


#####################################
# Functions for Paradigm data import

def read_readme(file_readme):
    with open(file_readme) as f:
        text = f.read()
        start_time      = text.split('Start time: ')[1].split('T')[1].split('.00\n')[0].replace(" ", "")
        stop_time       = text.split('Stop time: ')[1].split('T')[1].split('.00\n')[0].replace(" ", "")
        num_traces      = text.split('Number of traces:')[1].split('\n')[0].replace(" ", "")
        sample_interval = text.split('Resampled data sample interval in ns:')[1].split('\n')[0].replace(" ", "")
        twt_trace       = text.split('TWT of resampled full trace in ms:')[1].split('\n')[0].replace(" ", "")
        try:
            num_samples     = text.split('Number of resampled SGY file samples:')[1].split('\n')[0].replace(" ", "")
        except:
            print("'Number of samples per trace:' string not found... calculating num_samples instead...")
            num_samples = int( (float(twt_trace) * 1000) / float(sample_interval) + 1 )


        return num_traces, sample_interval, twt_trace, num_samples



def write_gin(sample_interval, num_samples, line_label, label_suffix, line_label_coords, segment, seissrv_sgy_path, gin_filename, file_sgy, line_folder):

    t_length = float(sample_interval) * float(num_samples)

    if int(num_samples) <= 1024:
        new_num_samples = 1024
    else:
        i = (int(num_samples) // 1024) + 1
        new_num_samples = 1024 * i

    new_sample_interval = int(t_length // new_num_samples)
    new_t_length        = int(new_num_samples * new_sample_interval)

    ts_buffer = " " * int(9 - len(str(int(t_length))))

    with open(line_folder + '/' + gin_filename, 'w') as gin:
        gin.write('*JOB    s/p_rada{}\n'.format(line_label))
        gin.write('** FLIGHT No/Segment {}\n'.format(segment))
        gin.write('** {}\n'.format(file_sgy))
        gin.write('** dt = {}        samples = {}\n'.format(sample_interval, num_samples))
        gin.write('**\n')
        gin.write('*CALL   GIN     {}{}{:.04f}          SHOT\n'.format(int(t_length), ts_buffer, float(sample_interval)))
        gin.write('TAPEOPT -tapefile {}\n'.format(seissrv_sgy_path))
        gin.write('DEFINE  SHOT    JPHYSIN\n')
        gin.write('REEL    1                               \n')
        gin.write('**\n')
        # gin.write('*CALL   RESAMP  {}       {}\n'.format(new_sample_interval, new_t_length))
        # gin.write('**\n')
        # gin.write('*CALL   HDRMATH\n')
        # gin.write('HCADD   shot    0       chan\n')
        # gin.write('**\n')
        # gin.write('*CALL   HDRMATH\n')
        # gin.write('HCMUL   hour    3600    desc\n')
        # gin.write('HCMUL   minute  60      min60\n')
        # gin.write('HHADD   desc    min60   dsec\n')
        # gin.write('HHADD   dsec    second  dsec\n')
        # gin.write('**\n')
        # gin.write('*CALL   HDRMATH\n')
        # gin.write('HCMUL   hour    10000   TIME\n')
        # gin.write('HCMUL   minute  100     min100\n')
        # gin.write('HHADD   TIME    min100  TIME\n')
        # gin.write('HHADD   TIME    second  TIME\n')
        # gin.write('**\n')
        # gin.write('*CALL   HDRMATH\n')
        # gin.write('HCMUL   WDEPTHRC-0.3    statcor\n')
        # gin.write('HCSUB   statcor 1000    statcor\n')
        # gin.write('**\n')
        gin.write('*CALL   DSOUT   OVERWRT\n')
        gin.write('LABEL   {}{}\n'.format(line_label_coords, label_suffix))
        gin.write('**\n')
        gin.write('*END')




def write_scale(sample_interval, num_samples, line_label, label_suffix, line_label_coords, segment, scale_filename, file_sgy, line_folder, scale_factor):

    with open(line_folder + '/' + scale_filename, 'w') as gin:
        gin.write('*JOB    s/p_rada{}\n'.format(line_label))
        gin.write('** FLIGHT No/Segment {}\n'.format(segment))
        gin.write('** {}\n'.format(file_sgy))
        gin.write('** dt = {}        samples = {}\n'.format(sample_interval, num_samples))
        gin.write('*CALL   DSIN\n')
        gin.write('LABEL   {}{}\n'.format(line_label_coords, label_suffix))
        gin.write('**\n')
        gin.write('*CALL   SCALE                   {}\n'.format(scale_factor))
        gin.write('**\n')
        gin.write('*CALL   DSOUT   OVERWRT\n')
        gin.write('LABEL   {}{}_scaled\n'.format(line_label_coords, label_suffix))   #####
        gin.write('**\n')
        gin.write('*END\n')



def write_agc(line_label, label_suffix, line_label_coords, agc_filename, line_folder, agc_gate):

    with open(line_folder + '/' + agc_filename, 'w') as gin:
        gin.write('*JOB    s/p_rada{}\n'.format(line_label))
        gin.write('*CALL   DSIN\n')
        gin.write('LABEL   {}{}_scaled\n'.format(line_label_coords, label_suffix))
        gin.write('**\n')
        gin.write('*CALL   AGC     {}\n'.format(agc_gate))
        gin.write('**\n')
        gin.write('*CALL   DSOUT   OVERWRT\n')
        gin.write('LABEL   {}{}_agc\n'.format(line_label_coords, label_suffix))   #####
        gin.write('**\n')
        gin.write('*END\n')



def get_ice_surf(file_ll, line_label_coords, line_folder, campaign):
    import pandas as pd

    df = pd.read_csv(file_ll, sep='\s+')

    df["Surface Name"] = "Surface_autotracked"
    df["Line Name"]    = line_label_coords
    df["CMP Label"]    = df["NR"]
    df["Shot Label"]   = df["NR"]
    df["X Coordinate"] = df["LONGITUDE"]
    df["Y Coordinate"] = df["LATITUDE"]
    df["Value"]        = df["SURF_TWT"] * 1000 * 1000 * 1000
    df["Segment ID"]   = line_label_coords

    df = df[["Surface Name", "Line Name", "CMP Label", "Shot Label", "X Coordinate", "Y Coordinate", "Value", "Segment ID"]]

    df.to_csv(line_folder + '/' + line_label_coords + '_icesurf.csv', sep='\t', index=False)

    return str(line_folder + '/' + line_label_coords + '_icesurf.csv')



def create_coord_file(file_ll, line_label, line_label_coords, line_folder):
    import pandas as pd

    df = pd.read_csv(file_ll, sep='\s+')
    df.insert(0, 'LINE_LABEL', line_label_coords)
    df = df[['LINE_LABEL', 'LINE', 'NR', 'LONGITUDE', 'LATITUDE']]
    
    df['LONGITUDE'] = df['LONGITUDE'].map(lambda x: '{0: .9f}'.format(float(x)))
    df['LATITUDE']  = df['LATITUDE'].map(lambda x: '{0: .9f}'.format(float(x)))

    df['NR'] = df['NR'].astype(str)
    df['NR'] = df['NR'].str.rjust(5, " ")

    df.to_csv(line_folder + '/' + line_label_coords + '_coords.csv', sep='\t', index=False)

    return str(line_folder + '/' + line_label_coords + '_coords.csv')



#########################################
# Functions for Horizon conversion

def paradigm_picks2csv(dir_paradigm_picks, dir_csv_picks, picks_format, layer_group, layer, filter, uwbm_only):

    '''
    ..........
    
    '''

    import pandas as pd
    import os, glob
    import copy


    if uwbm_only == True:
        suffix = "_uwbm_only"
    else:
        suffix = ""


    if layer == "Base_master" or layer == "Surface_master" or layer == "Base_from_ice":
        out_dir = "{}\\{}".format(dir_csv_picks, layer)
    else:
        out_dir = "{}\\IRHs\\{}".format(dir_csv_picks, layer)

    pick_files = sorted(glob.glob("{}\\{}\\{}*{}*{}{}.txt".format(dir_paradigm_picks, layer_group, layer, filter, picks_format, suffix)))

    for file in pick_files:
        print(file)
        
        if picks_format == "UKOOA":
            df               = pd.read_csv(file, sep="\s+",  header=None)
            df.columns       = ["paradigm_id", "trace", "longitude", "latitude", "twt"]
            df["season"]     = "antr" + str(df["paradigm_id"].iloc[0])[0:4]
            df["profile_id"] = df["paradigm_id"]
            
            df["twt"]    = df["twt"] / 1000 / 1000 / 1000
            df           = df.sort_values(by=["season", "paradigm_id", "trace"])

            df['longitude'] = df['longitude'].astype(float).apply(lambda x: "{:.7f}".format(x))
            df['latitude']  = df['latitude'].astype(float).apply(lambda x: "{:.7f}".format(x))
            df['twt']       = df['twt'].astype(float).apply(lambda x: "{:.12f}".format(x))
            df['trace']       = df['trace'].astype(int).apply(lambda x: "{:5d}".format(x))
            
            df = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]
        
        if picks_format == "GQC7":
            df           = pd.read_csv(file, sep="\s+", skiprows=2, header=None)
            df           = df[:-1]
            df.columns   = ["longitude", "latitude", "dm1", "dm2", "twt", "dm3", "dm4", "trace", "dm5", "id"]
            df           = df[df["longitude"].astype(str).str.contains("EOD|PROFILE|SNAPPING") == False]
            
            
            
            # UWB antr2025
            if "antr2025_20257" in df["id"].iloc[0]:
                xs = df["id"].str.split("antr2025_20257_", expand = True)
                xs.columns   = ["season_nom", "profile_id"]

                df["season"]      = "antr2025"
                df["prefix"]      = "20257_"
                df["profile_id"] = xs["profile_id"]
                df["paradigm_id"]  = df["prefix"] + df["profile_id"]
            
            
            elif "arkr2024_20247" in df["id"].iloc[0]:
                xs = df["id"].str.split("arkr2024_20247_", expand = True)
                xs.columns   = ["season_nom", "profile_id"]

                df["season"]      = "arkr2024"
                df["prefix"]      = "20247_"
                df["profile_id"] = xs["profile_id"]
                df["paradigm_id"]  = df["prefix"] + df["profile_id"]
            
            # UWB antr2024
            elif "antr2024_20247" in df["id"].iloc[0]:
                xs = df["id"].str.split("antr2024_20247_", expand = True)
                xs.columns   = ["season_nom", "profile_id"]

                df["season"]      = "antr2024"
                df["prefix"]      = "20247_"
                df["profile_id"] = xs["profile_id"]
                df["paradigm_id"]  = df["prefix"] + df["profile_id"]


            # UWB arkr2022
            elif "arkr2022_20227" in df["id"].iloc[0]:
                xs = df["id"].str.split("arkr2022_20227_", expand = True)
                xs.columns   = ["season_nom", "profile_id"]

                df["season"]      = "arkr2022"
                df["prefix"]      = "20227_"
                df["profile_id"] = xs["profile_id"]
                df["paradigm_id"]  = df["prefix"] + df["profile_id"]
                
            # UWB antr2019
            elif "antr2019_20197" in df["id"].iloc[0]:
                xs = df["id"].str.split("antr2019_20197_", expand = True)
                xs.columns   = ["season_nom", "profile_id"]

                df["season"]      = "antr2019"
                df["prefix"]      = "20197_"
                df["profile_id"] = xs["profile_id"]
                df["paradigm_id"]  = df["prefix"] + df["profile_id"]

            # UWBM arkr2023
            elif "arkr2023_20237" in df["id"].iloc[0]:
                xs = df["id"].str.split("arkr2023_20237_", expand = True)
                xs.columns   = ["season_nom", "profile_id"]

                df["season"]      = "arkr2023"
                df["prefix"]      = "20237_"
                df["profile_id"] = xs["profile_id"]
                df["paradigm_id"]  = df["prefix"] + df["profile_id"]

            # UWBM antr2023
            elif "antr2023_20238" in df["id"].iloc[0] and uwbm_only == True:
                xs = df["id"].str.split("antr2023_20238_", expand = True)
                xs.columns   = ["season_nom", "profile_id"]

                df["season"]      = "antr2023"
                df["prefix"]      = "20238_"
                df["profile_id"] = xs["profile_id"]
                df["paradigm_id"]  = df["prefix"] + df["profile_id"]

            ### antr2017 contains emr and uwb
            elif "antr2017" in df["id"].iloc[0]:

                # split emr and uwb part
                mask   = df["id"].str.contains("20177")
                df_uwb = df[mask]
                df_emr = df[~mask]

                if len(df_uwb) > 0:
                    # handle uwb part
                    xs_uwb = df_uwb["id"].str.split("antr2017_20177_", expand = True)
                    xs_uwb.columns   = ["season_nom", "profile_id"]

                    df_uwb["season"]      = "antr2017"
                    df_uwb["prefix"]      = "20177_"
                    df_uwb["profile_id"] = xs_uwb["profile_id"]
                    df_uwb["paradigm_id"]  = df_uwb["prefix"] + df_uwb["profile_id"]

                    # handle emr part
                    xs_emr         = df_emr["id"].str.split("_", expand = True)
                    xs_emr.columns = ["season", "paradigm_id"]

                    df_emr["season"]      = "antr2017"
                    df_emr["profile_id"]  = xs_emr["paradigm_id"]
                    df_emr["paradigm_id"] = xs_emr["paradigm_id"]

                    df = pd.concat([df_emr, df_uwb]).reset_index(drop=True)
                
                else:
                    xs                = df["id"].str.split("_", expand = True)
                    xs.columns        = ["season", "paradigm_id"]
                    df["season"]      = "antr2017"
                    df["profile_id"]  = xs["paradigm_id"]
                    df["paradigm_id"] = xs["paradigm_id"]

            #### antr2023 contains emr and uwbm
            elif "antr2023" in df["id"].iloc[0] and uwbm_only == False:

                # split emr and uwb part
                mask   = df["id"].str.contains("20238")
                df_uwb = df[mask]
                df_emr = df[~mask]

                if len(df_uwb) > 0:
                    # handle uwb part
                    xs_uwb = df_uwb["id"].str.split("antr2023_20238_", expand = True)
                    xs_uwb.columns   = ["season_nom", "profile_id"]

                    df_uwb["season"]      = "antr2023"
                    df_uwb["prefix"]      = "20238_"
                    df_uwb["profile_id"] = xs_uwb["profile_id"]
                    df_uwb["paradigm_id"]  = df_uwb["prefix"] + df_uwb["profile_id"]

                    # handle emr part
                    xs_emr         = df_emr["id"].str.split("_", expand = True)
                    xs_emr.columns = ["season", "paradigm_id"]

                    df_emr["season"]      = "antr2023"
                    df_emr["profile_id"] = "None"
                    df_emr["paradigm_id"]  = xs_emr["paradigm_id"]

                    df = pd.concat([df_emr, df_uwb]).reset_index(drop=True)
                
                else:
                    xs                = df["id"].str.split("_", expand = True)
                    xs.columns        = ["season", "paradigm_id"]
                    df["season"]      = "antr2023"
                    df["profile_id"]  = xs["paradigm_id"]
                    df["paradigm_id"] = xs["paradigm_id"]

            else:
                xs = df["id"].str.split("_", expand = True)

                try:
                    xs.columns   = ["season", "paradigm_id", "subline"]
                except:
                    xs["subline"] = "None"
                    xs.columns    = ["season", "paradigm_id", "subline"]

                    df["season"] = xs["season"].iloc[0]
                    df["paradigm_id"]   = xs["paradigm_id"].map(str) + "_" + xs["subline"].map(str)
                    df["paradigm_id"]   = df["paradigm_id"].str.replace("_None", "")
                    df["profile_id"]   = copy.copy(df["paradigm_id"])

            if "antr2017" in df["season"].iloc[0]:
                df = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]

            elif "arkr2022" in df["season"].iloc[0]:
                df = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]

            elif "antr2024" in df["season"].iloc[0]:
                df = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]

            elif "antr2025" in df["season"].iloc[0]:
                df = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]
            
            elif "arkr2023" in df["season"].iloc[0]:
                df = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]
            
            elif "antr2023" in df["season"].iloc[0]:
                df = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]
                
            elif "antr2019" in df["season"].iloc[0]:
                df = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]

            else:
                df["profile_id"]   = copy.copy(df["paradigm_id"])
                df                 = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]
            
            df["twt"]    = df["twt"] / 1000 / 1000 / 1000
            df           = df.sort_values(by=["season", "paradigm_id", "trace"])

            df['longitude'] = df['longitude'].astype(float).apply(lambda x: "{:.7f}".format(x))
            df['latitude']  = df['latitude'].astype(float).apply(lambda x: "{:.7f}".format(x))
            df['twt']       = df['twt'].astype(float).apply(lambda x: "{:.12f}".format(x))
            df['trace']       = df['trace'].astype(int).apply(lambda x: "{:5d}".format(x))
            df              = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]


        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        group_list = []
        groups = df.groupby("profile_id")

        for name, group in groups:
            #try:
            season = group["season"].iloc[0]
            line   = name


            
            if "antr2024" in season:
                profile_id = group["profile_id"].iloc[0]
                line        = profile_id

            if "antr2025" in season:
                profile_id = group["profile_id"].iloc[0]
                line        = profile_id

            elif "arkr2022" in season:
                profile_id = group["profile_id"].iloc[0]
                line        = profile_id
                
            elif "antr2019" in season:
                profile_id = group["profile_id"].iloc[0]
                line        = profile_id

            elif "arkr2023" in season:
                profile_id = group["profile_id"].iloc[0]
                line        = profile_id


            elif "antr2017" in season:
                if group["profile_id"].iloc[0] == "None":
                    line = name
                    del group["profile_id"]
                else:
                    profile_id = group["profile_id"].iloc[0]
                    line        = profile_id

            elif "antr2023" in season:
                if group["profile_id"].iloc[0] == "None":
                    line = name
                    del group["profile_id"]
                else:
                    profile_id = group["profile_id"].iloc[0]
                    line        = profile_id
        
            group_list.append(group)
            print("Saving: {}\\{}_{}.csv".format(out_dir, layer, line))
            group.to_csv("{}\\{}_{}.csv".format(out_dir, layer, line), sep="\t", index=False)