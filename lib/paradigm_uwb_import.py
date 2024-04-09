

def read_readme(sgy_path, file_readme):
    with open(file_readme) as f:
        text = f.read()
        start_time      = text.split('Start time: ')[1].split('T')[1].split('.00\n')[0].replace(" ", "")
        stop_time       = text.split('Stop time: ')[1].split('T')[1].split('.00\n')[0].replace(" ", "")
        num_traces      = text.split('Number of traces:')[1].split('\nSample')[0].replace(" ", "")
        sample_interval = text.split('Sample interval in ns:')[1].split('\n')[0].replace(" ", "")
        twt_trace       = text.split('TWT of full trace in ms:')[1].split('\n')[0].replace(" ", "")
        num_samples     = text.split('Number of samples per trace:')[1].split('\n')[0].replace(" ", "")


        return start_time, stop_time, num_traces, sample_interval, twt_trace, num_samples



def write_gin(sample_interval, num_samples, line_label, label_suffix, line_label_coords, segment, seissrv_sgy_path, gin_filename, file_sgy, line_folder):

    t_length = float(sample_interval) * float(num_samples)

    if int(num_samples) <= 1024:
        new_num_samples = 1024
    else:
        i = (int(num_samples) // 1024) + 1
        new_num_samples = 1024 * i

    new_sample_interval = int(t_length // new_num_samples)
    new_t_length        = int(new_num_samples * new_sample_interval)

    with open(line_folder + '/' + gin_filename, 'w') as gin:
        gin.write('*JOB    s/p_rada{}\n'.format(line_label))
        gin.write('** FLIGHT No/Segment {}\n'.format(segment))
        gin.write('** {}\n'.format(file_sgy))
        gin.write('** dt = {}        samples = {}\n'.format(sample_interval, num_samples))
        gin.write('**\n')
        gin.write('*CALL   GIN     {}     {:.04f}          SHOT\n'.format(int(t_length), float(sample_interval)))
        gin.write('TAPEOPT -tapefile {}\n'.format(seissrv_sgy_path))
        gin.write('DEFINE  SHOT    JPHYSIN\n')
        gin.write('REEL    1                               \n')
        gin.write('**\n')
        # gin.write('*CALL   RESAMP  {}       {}\n'.format(new_sample_interval, new_t_length))
        # gin.write('**\n')
        gin.write('*CALL   HDRMATH\n')
        gin.write('HCADD   shot    0       chan\n')
        gin.write('**\n')
        gin.write('*CALL   HDRMATH\n')
        gin.write('HCMUL   hour    3600    desc\n')
        gin.write('HCMUL   minute  60      min60\n')
        gin.write('HHADD   desc    min60   dsec\n')
        gin.write('HHADD   dsec    second  dsec\n')
        gin.write('**\n')
        gin.write('*CALL   HDRMATH\n')
        gin.write('HCMUL   hour    10000   TIME\n')
        gin.write('HCMUL   minute  100     min100\n')
        gin.write('HHADD   TIME    min100  TIME\n')
        gin.write('HHADD   TIME    second  TIME\n')
        gin.write('**\n')
        gin.write('*CALL   HDRMATH\n')
        gin.write('HCMUL   WDEPTHRC-0.3    statcor\n')
        gin.write('HCSUB   statcor 1000    statcor\n')
        gin.write('**\n')
        gin.write('*CALL   DSOUT   OVERWRT\n')
        gin.write('LABEL   {}\n'.format(line_label_coords))
        gin.write('**\n')
        gin.write('*END')




# def write_gin(sample_interval, num_samples, line_label, label_suffix, line_label_coords, segment, seissrv_sgy_path, gin_filename, file_sgy, line_folder):

#     import numpy as np

#     t_length = float(sample_interval) * float(num_samples)

#     if int(num_samples) <= 1024:
#         new_num_samples = 1024
#     else:
#         i = (int(num_samples) // 1024) + 1
#         new_num_samples = 1024 * i

#     new_sample_interval = int(t_length // new_num_samples)
#     new_t_length        = int(new_num_samples * new_sample_interval)

#     with open(line_folder + '/' + gin_filename, 'w') as gin:
#         gin.write('*JOB    s/p_rada{}\n'.format(line_label))
#         gin.write('** FLIGHT No/Segment {}\n'.format(segment))
#         gin.write('** {}\n'.format(file_sgy))
#         gin.write('** dt = {}\t\tsamples = {}\n'.format(sample_interval, num_samples))
#         gin.write('**\n')
#         gin.write('*CALL\tGIN\t\t{}\t{:.04f}\t\t\tSHOT\n'.format(int(t_length), float(sample_interval)))
#         gin.write('TAPEOPT -tapefile {}\n'.format(seissrv_sgy_path))
#         gin.write('DEFINE\tSHOT\tJPHYSIN\n')
#         gin.write('REEL\t1\t\t\t\t\t\t\t\t{}\n'.format(line_label_coords))
#         gin.write('**\n')
#         gin.write('*CALL\tRESAMP\t{}\t\t{}\n'.format(new_sample_interval, new_t_length))
#         gin.write('**\n')
#         gin.write('*CALL\tHDRMATH\n')
#         gin.write('HCADD\tshot\t0\t\tchan\n')
#         gin.write('**\n')
#         gin.write('*CALL\tHDRMATH\n')
#         gin.write('HCMUL\thour\t3600\tdesc\n')
#         gin.write('HCMUL\tminute\t60\t\tmin60\n')
#         gin.write('HHADD\tdesc\tmin60\tdsec\n')
#         gin.write('HHADD\tdsec\tsecond\tdsec\n')
#         gin.write('**\n')
#         gin.write('*CALL\tHDRMATH\n')
#         gin.write('HCMUL\thour\t10000\tTIME\n')
#         gin.write('HCMUL\tminute\t100\t\tmin100\n')
#         gin.write('HHADD\tTIME\tmin100\tTIME\n')
#         gin.write('HHADD\tTIME\tsecond\tTIME\n')
#         gin.write('**\n')
#         gin.write('*CALL\tHDRMATH\n')
#         gin.write('HCMUL\tWDEPTHRC-0.3\tstatcor\n')
#         gin.write('HCSUB\tstatcor\t1000\tstatcor\n')
#         gin.write('**\n')
#         gin.write('*CALL\tDSOUT\tOVERWRT\n')
#         gin.write('LABEL\t{}_{}\n'.format(line_label_coords, label_suffix))
#         gin.write('**\n')
#         gin.write('*END')



def write_scale(sample_interval, num_samples, line_label, label_suffix, line_label_coords, segment, scale_filename, file_sgy, line_folder, scale_factor):

    with open(line_folder + '/' + scale_filename, 'w') as gin:
        gin.write('*JOB    s/p_rada{}\n'.format(line_label))
        gin.write('** FLIGHT No/Segment {}\n'.format(segment))
        gin.write('** {}\n'.format(file_sgy))
        gin.write('** dt = {}        samples = {}\n'.format(sample_interval, num_samples))
        gin.write('*CALL   DSIN\n')
        gin.write('LABEL   {}\n'.format(line_label_coords))
        gin.write('**\n')
        gin.write('*CALL   SCALE                   {}\n'.format(scale_factor))
        gin.write('**\n')
        gin.write('*CALL   DSOUT   OVERWRT\n')
        gin.write('LABEL   {}_scaled\n'.format(line_label_coords))   #####
        gin.write('**\n')
        gin.write('*END\n')



def write_agc(line_label, label_suffix, line_label_coords, agc_filename, line_folder):

    with open(line_folder + '/' + agc_filename, 'w') as gin:
        gin.write('*JOB    s/p_rada{}\n'.format(line_label))
        gin.write('*CALL   DSIN\n')
        gin.write('LABEL   {}_scaled\n'.format(line_label_coords))
        gin.write('**\n')
        gin.write('*CALL   AGC     100\n')
        gin.write('**\n')
        gin.write('*CALL   DSOUT   OVERWRT\n')
        gin.write('LABEL   {}_agc\n'.format(line_label_coords))   #####
        gin.write('**\n')
        gin.write('*END\n')



def get_ice_surf(file_ll, line_label, sample_interval, line_folder):
    import pandas as pd

    df = pd.read_csv(file_ll, delim_whitespace=True)
    df['LINE_LABEL'] = line_label
    df['RT_SURF']    = df['RT_SURF'].astype(float) * float(sample_interval) * -1
    df['RT_BED']    = df['RT_BED'].astype(float) * float(sample_interval) * -1
    df = df[['LINE_LABEL', 'NR', 'LONGITUDE', 'LATITUDE', 'RT_SURF', 'RT_BED']]

    df['NR'] = df['NR'].astype(str)
    df['NR'] = df['NR'].str.rjust(5, " ")

    df['LONGITUDE'] = df['LONGITUDE'].map(lambda x: '{0:.9f}'.format(float(x)))
    df['LATITUDE']  = df['LATITUDE'].map(lambda x: '{0:.9f}'.format(float(x)))
    df['RT_SURF']   = df['RT_SURF'].map(lambda x: '{0:.2f}'.format(float(x)))
    df['RT_BED']    = df['RT_BED'].map(lambda x: '{0:.2f}'.format(float(x)))

    df.to_csv(line_folder + '/' + line_label + '_icesurf.csv', sep='\t', index=False)



def create_coord_file(file_ll, line_label, line_label_coords, line_folder):
    import pandas as pd

    df = pd.read_csv(file_ll, delim_whitespace=True)
    df.insert(0, 'LINE_LABEL', line_label_coords)
    df = df[['LINE_LABEL', 'LINE', 'NR', 'LONGITUDE', 'LATITUDE']]
    
    df['LONGITUDE'] = df['LONGITUDE'].map(lambda x: '{0: .9f}'.format(float(x)))
    df['LATITUDE']  = df['LATITUDE'].map(lambda x: '{0: .9f}'.format(float(x)))

    df['NR'] = df['NR'].astype(str)
    df['NR'] = df['NR'].str.rjust(5, " ")

    df.to_csv(line_folder + '/' + line_label_coords + '_coords.csv', sep='\t', index=False)

    return str(line_folder + '/' + line_label_coords + '_coords.csv')