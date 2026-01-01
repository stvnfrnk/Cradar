


#####################################
# Functions for Paradigm data import

def read_readme(file_readme):
    with open(file_readme) as f:
        text = f.read()
        start_time      = text.split('Start time: ')[1].split('T')[1].split('.00\n')[0].replace(" ", "")
        stop_time       = text.split('Stop time: ')[1].split('T')[1].split('.00\n')[0].replace(" ", "")

        try:
            num_traces      = text.split('Number of traces:')[1].split('\n')[0].replace(" ", "")
        except:
            text.split('Number of traces:')[1].split('\n')[0].replace("         ", "")

        try:
            sample_interval = text.split('Resampled data sample interval in ns:')[1].split('\n')[0].replace(" ", "")
        except:
            sample_interval = text.split('Raw data Sample interval in ns:')[1].split('\n')[0].replace(" ", "")

        try:
            twt_trace       = text.split('TWT of resampled full trace in ms:')[1].split('\n')[0].replace(" ", "")
        except:
            twt_trace       = text.split('TWT of full trace in ms:')[1].split('\n')[0].replace(" ", "")
        
        try:
            try:
                num_samples     = text.split('Number of resampled SGY file samples:')[1].split('\n')[0].replace(" ", "")
            except:
                num_samples     = text.split('Number of raw data samples:         ')[1].split('\n')[0]#.replace("         ", "")
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
        gin.write('*CALL   DSOUT   OVERWRT\n')
        gin.write('LABEL   {}{}\n'.format(line_label_coords, label_suffix))
        gin.write('**\n')
        gin.write('*END')


def write_gin_agc(sample_interval, num_samples, line_label, label_suffix, line_label_coords, segment, seissrv_segy_agc_path, gin_agc_filename, file_sgy_agc, line_folder):

    t_length = float(sample_interval) * float(num_samples)

    if int(num_samples) <= 1024:
        new_num_samples = 1024
    else:
        i = (int(num_samples) // 1024) + 1
        new_num_samples = 1024 * i

    new_sample_interval = int(t_length // new_num_samples)
    new_t_length        = int(new_num_samples * new_sample_interval)

    ts_buffer = " " * int(9 - len(str(int(t_length))))

    with open(line_folder + '/' + gin_agc_filename, 'w') as gin:
        gin.write('*JOB    s/p_rada{}\n'.format(line_label))
        gin.write('** FLIGHT No/Segment {}\n'.format(segment))
        gin.write('** {}\n'.format(file_sgy_agc))
        gin.write('** dt = {}        samples = {}\n'.format(sample_interval, num_samples))
        gin.write('**\n')
        gin.write('*CALL   GIN     {}{}{:.04f}          SHOT\n'.format(int(t_length), ts_buffer, float(sample_interval)))
        gin.write('TAPEOPT -tapefile {}\n'.format(seissrv_segy_agc_path))
        gin.write('DEFINE  SHOT    JPHYSIN\n')
        gin.write('REEL    1                               \n')
        gin.write('**\n')
        gin.write('*CALL   DSOUT   OVERWRT\n')
        gin.write('LABEL   {}{}_AGC_mapgis\n'.format(line_label_coords, label_suffix))
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


def get_ice_base(file_ll, line_label_coords, line_folder, campaign):
    import pandas as pd

    df = pd.read_csv(file_ll, sep='\s+')

    df["Surface Name"] = "Base_autotracked"
    df["Line Name"]    = line_label_coords
    df["CMP Label"]    = df["NR"]
    df["Shot Label"]   = df["NR"]
    df["X Coordinate"] = df["LONGITUDE"]
    df["Y Coordinate"] = df["LATITUDE"]
    df["Value"]        = df["BED_TWT"] * 1000 * 1000 * 1000
    df["Segment ID"]   = line_label_coords

    df = df[["Surface Name", "Line Name", "CMP Label", "Shot Label", "X Coordinate", "Y Coordinate", "Value", "Segment ID"]]

    df.to_csv(line_folder + '/' + line_label_coords + '_icebase.csv', sep='\t', index=False)

    return str(line_folder + '/' + line_label_coords + '_icebase.csv')



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

def paradigm_picks2csv(dir_paradigm_picks, dir_csv_picks, dir_ll_files, picks_format, layer_group, layer, filter, AWI_radar_metadata):

    '''
    
    
    '''

    import pandas as pd
    import numpy as np
    import os, glob
    import copy

    pd.options.mode.chained_assignment = None  # default='warn'

    if layer == "Base_master" or layer == "Surface_master" or layer == "Base_from_ice":
        out_dir = "{}\\{}".format(dir_csv_picks, layer)
    else:
        out_dir = "{}\\IRHs\\{}".format(dir_csv_picks, layer)


    pick_files   = sorted(glob.glob("{}\\{}\\{}*{}*{}.txt".format(dir_paradigm_picks, layer_group, layer, filter, picks_format)))
    df_metadata  = pd.read_csv(AWI_radar_metadata, sep=";")

    for file in pick_files:
        print("")
        print("Loading: {}".format(file))
        
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
            df           = pd.read_csv(file, sep="\s+", skiprows=2, header=None, low_memory=False)
            df           = df[:-1]
            df.columns   = ["longitude", "latitude", "dm1", "dm2", "twt", "dm3", "dm4", "trace", "dm5", "id"]
            df           = df[df["longitude"].astype(str).str.contains("EOD|PROFILE|SNAPPING") == False]
            

            #############################################################
            
            if "EMR" in df["id"].iloc[0]:
                season = str(df["id"].iloc[0][5:13])
                year   = str(df["id"].iloc[0][0:4])
            else:
                season = str(df["id"].iloc[0][0:8])
                year   = str(df["id"].iloc[0][9:13])
                        
            # get EMR
            mask      = df["id"].str.contains("{}_{}2".format(season, year))
            df_emr600 = df[mask]
            mask      = df["id"].str.contains("{}_{}3".format(season, year))
            df_emr60  = df[mask]
            df_emr    = pd.concat([df_emr600, df_emr60])
            
            # get EMR (MODAMS)
            mask      = df["id"].str.contains("EMR".format(season, year))
            df_modams = df[mask]
            
            # get ACCU
            mask      = df["id"].str.contains("{}_{}4_".format(season, year))
            df_accu   = df[mask]
            
            # get SNOW
            mask      = df["id"].str.contains("{}_{}5_".format(season, year))
            df_snow   = df[mask]
            
            # get UWB
            mask     = df["id"].str.contains("{}_{}7_".format(season, year))
            df_uwb   = df[mask]
            
            # get UWBM
            mask     = df["id"].str.contains("{}_{}8_".format(season, year))
            df_uwbm  = df[mask]
            
            list_df  = []

            # handle emr part
            if len(df_emr) > 0:
                xs_emr                = df_emr["id"].str.split("_", expand = True)
                # xs_emr.to_csv("C:\\Users\\sfranke\\Desktop\\xs_emr.csv", sep="\t", index=False)
                xs_emr.columns        = ["season_nom", "profile_id"]
                df_emr["season"]      = season
                df_emr["profile_id"]  = xs_emr["profile_id"]
                df_emr["paradigm_id"] = xs_emr["profile_id"]
                list_df.append(df_emr)
                
            # handle emr (modams) part
            if len(df_modams) > 0:
                xs_modams                = df_modams["id"].str.split("_", expand = True)
                xs_modams.columns        = ["year", "season_nom", "flight_id", "pulse", "date", "part"]
                df_modams["season"]      = season
                df_modams["profile_id"]  = xs_modams["season_nom"] + "_" + xs_modams["flight_id"] + "_" + xs_modams["pulse"] + "_" + xs_modams["date"] + "_" + xs_modams["part"]
                df_modams["paradigm_id"] = xs_modams["season_nom"] + "_" + xs_modams["flight_id"] + "_" + xs_modams["pulse"] + "_" + xs_modams["date"] + "_" + xs_modams["part"]
                list_df.append(df_modams)
                
            # handle accu part
            if len(df_accu) > 0:
                xs_accu                = df_accu["id"].str.split("{}_{}4_".format(season, year), expand = True)
                # xs_accu.to_csv("C:\\Users\\sfranke\\Desktop\\xs_accu.csv", sep="\t", index=False)
                xs_accu.columns        = ["season_nom", "profile_id"]
                df_accu["season"]      = season
                df_accu["prefix"]      = "{}4_".format(year)
                df_accu["profile_id"]  = xs_accu["profile_id"]
                df_accu["paradigm_id"] = df_accu["prefix"] + df_accu["profile_id"]
                df_accu.to_csv("C:\\Users\\sfranke\\Desktop\\df_accu.csv", sep="\t", index=False)
                list_df.append(df_accu)
                
            # handle snow part
            if len(df_snow) > 0:
                xs_snow                = df_snow["id"].str.split("{}_{}5_".format(season, year), expand = True)
                xs_snow.columns        = ["season_nom", "profile_id"]
                df_snow["season"]      = season
                df_snow["prefix"]      = "{}5_".format(year)
                df_snow["profile_id"]  = xs_snow["profile_id"]
                df_snow["paradigm_id"] = df_snow["prefix"] + df_snow["profile_id"]
                list_df.append(df_snow)
                
            # handle uwb part
            if len(df_uwb) > 0:
                xs_uwb                = df_uwb["id"].str.split("{}_{}7_".format(season, year), expand = True)
                xs_uwb.columns        = ["season_nom", "profile_id"]
                df_uwb["season"]      = season
                df_uwb["prefix"]      = "{}7_".format(year)
                df_uwb["profile_id"]  = xs_uwb["profile_id"]
                df_uwb["paradigm_id"] = df_uwb["prefix"] + df_uwb["profile_id"]
                list_df.append(df_uwb)
            
            # handle uwbm part
            if len(df_uwbm) > 0:
                xs_uwbm                = df_uwbm["id"].str.split("{}_{}8_".format(season, year), expand = True)
                xs_uwbm.columns        = ["season_nom", "profile_id"]
                df_uwbm["season"]      = season
                df_uwbm["prefix"]      = "{}8_".format(year)
                df_uwbm["profile_id"]  = xs_uwbm["profile_id"]
                df_uwbm["paradigm_id"] = df_uwbm["prefix"] + df_uwbm["profile_id"]
                list_df.append(df_uwbm)
                
                
            # combine data frames
            df = pd.concat(list_df).reset_index(drop=True)
            

            # # list of UWB or UWBM seasons
            # list_UWB_M_seasons = ["antr2017", "antr2019", "antr2023", "antr2024", "antr2025",
            #                       "arkr2012", "arkr2016", "arkr2018", "arkr2021", "arkr2022", "arkr2023", "arkr2024"]

            # # profile_id is different to paradigm_id (UWB, UWBM, ACCU, ASIRAS data)
            # if any(x in df["season"].iloc[0] for x in list_UWB_M_seasons):
            #     df = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]

            # # profile_id is the same as paradigm_id (EMR data)
            # else:
            #     df["profile_id"]   = copy.copy(df["paradigm_id"])
            #     df                 = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]
            
            df["twt"]       = df["twt"] / 1000 / 1000 / 1000
            df              = df.sort_values(by=["season", "paradigm_id", "trace"])
            df              = df[["season", "paradigm_id", "profile_id", "trace", "longitude", "latitude", "twt"]]

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        group_list = []
        groups     = df.groupby("paradigm_id")

        for name, group in groups:
            try:
                line = name
                # group_list.append(group)
                if len(df[df.duplicated(subset=["profile_id","trace"], keep=False)]) != 0:
                    group.drop_duplicates(subset=["profile_id", "trace"], keep="first", inplace=True)
                    
                # get lon lat from corresponding .ll file
                # res_id = int(group["paradigm_id"].iloc[0][4])
                pid    = str(group["profile_id"].iloc[0])
                

                ll_file = np.array(df_metadata["LL File Path Win"].loc[df_metadata.index[df_metadata["Profile ID"] == pid]])[0]                                
                
                # if "EMR" in pid:
                #     res_system = "MODAMS"
                #     ll_file    = glob.glob("{}\\{}\\*\\*{}*.ll".format(dir_ll_files, res_system, pid))[0]
                # if res_id == 2:
                #     res_system = "EMR"
                #     ll_file    = glob.glob("{}\\{}\\*\\*{}*.ll".format(dir_ll_files, res_system, pid))[0]
                # if res_id == 3:
                #     res_system = "EMR"
                #     ll_file    = glob.glob("{}\\{}\\*\\*{}*.ll".format(dir_ll_files, res_system, pid))[0]
                # if res_id == 4:
                #     res_system = "ACCU"
                #     # print("{}\\{}\\*\\ACCU_{}*.ll".format(dir_ll_files, res_system, pid))
                #     # group.to_csv("C:\\Users\\sfranke\\Desktop\\group.csv", sep="\t", index=False)
                #     ll_file    = glob.glob("{}\\{}\\*\\ACCU_{}*.ll".format(dir_ll_files, res_system, pid))[0]
                # if res_id == 5:
                #     res_system = "SNOW"
                #     ll_file    = glob.glob("{}\\{}\\*\\SNOW_{}*.ll".format(dir_ll_files, res_system, pid))[0]
                # if res_id == 7:
                #     res_system = "UWB"
                #     ll_file    = glob.glob("{}\\{}\\*_standard\\*\\*{}*.ll".format(dir_ll_files, res_system, pid))[0]
                # if res_id == 8:
                #     res_system = "UWBM"
                #     ll_file    = glob.glob("{}\\{}\\*\\*\\*{}*.ll".format(dir_ll_files, res_system, pid))[0]
                
                df_ll               = pd.read_csv(ll_file, sep="\s+")
                df_ll["profile_id"] = pid
                df_ll["NR"]         = df_ll["NR"].astype(int)
                
                group["profile_id"] = group["profile_id"].astype(str)
                group["trace"]      = group["trace"].astype(int)
                
                # df_ll.to_csv("C:\\Users\\sfranke\\Desktop\\tmp\\{}_ll_file.csv".format(line), sep="\t", index=False)
                # group.to_csv("C:\\Users\\sfranke\\Desktop\\tmp\\{}_group.csv".format(line), sep="\t", index=False)
                
                #print(df_ll.dtypes)
                #print(group.dtypes)
                
                # Merge df1 with df2 on 'profile_id' and 'trace'/'NR'
                df_out              = group.merge(df_ll[['profile_id', 'NR', 'LONGITUDE', 'LATITUDE']], left_on=['profile_id', 'trace'], right_on=['profile_id', 'NR'], how='left')
                # df_out.to_csv("C:\\Users\\sfranke\\Desktop\\tmp\\{}_merged.csv".format(line), sep="\t", index=False)
                df_out['longitude'] = df_out['LONGITUDE']
                df_out['latitude']  = df_out['LATITUDE']
                df_out              = df_out.drop(columns=['LONGITUDE', 'LATITUDE', 'NR'])
                
                df_out['longitude'] = df_out['longitude'].astype(float).apply(lambda x: "{:.9f}".format(x))
                df_out['latitude']  = df_out['latitude'].astype(float).apply(lambda x: "{:.9f}".format(x))
                df_out['twt']       = df_out['twt'].astype(float).apply(lambda x: "{:.12f}".format(x))
                df_out['trace']     = df_out['trace'].astype(int).apply(lambda x: "{:5d}".format(x))

                print("Saving: {}\\{}_{}.csv".format(out_dir, layer, line))
                df_out.to_csv("{}\\{}_{}.csv".format(out_dir, layer, line), sep="\t", index=False)
                
                
                # print("Bla Bla Bla: {}\\{}_{}.csv".format(out_dir, layer, line))
                
            except:
                print("Problem with: {}".format(pid))