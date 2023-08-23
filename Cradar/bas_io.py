def correct_chirp_data(nc_file):

    import xarray as xr
    import numpy as np
    import pandas as pd


    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]


    ds = xr.open_dataset(nc_file, decode_times=False)

    Time         = ds['fast_time'].values.astype(float) * 10e-7

    # Data_Pulse   = ds.variables['pulse_data'].values
    Data_Pulse   = ds.variables['polarised_pulse_SPHV_data'].values
    # Data_Chirp   = ds.variables['chirp_data'].values
    Data_Chirp   = ds.variables['polarised_chirp_SSHH_data'].values

    Traces_Chirp = ds.traces_chirp.values
    Traces_Pulse = ds.traces_pulse.values

    Longitude    = ds.variables['longitude_layerData'].values
    Latitude     = ds.variables['latitude_layerData'].values
    Elevation    = ds.variables['aircraft_altitude_layerData'].values
    GPS_time     = ds.variables['UTC_time_layerData'].values

    Surface_idx  = ds.variables['surface_pick_layerData']
    Surface_m    = ds.variables['surface_altitude_layerData'].values
    Bed_idx      = ds.variables['bed_pick_layerData']
    Bed_m        = ds.variables['bed_altitude_layerData'].values

    pri_pulse    = ds.variables['PriNumber_pulse'].values
    pri_chirp    = ds.variables['PriNumber_chirp'].values

    ### data frame with relevant variables ###

    df                 = pd.DataFrame(Traces_Chirp)
    df['Traces_Pulse'] = pd.DataFrame(Traces_Pulse)
    df['Longitude']    = pd.DataFrame(Longitude)
    df['Latitude']     = pd.DataFrame(Latitude)
    df['Elevation']    = pd.DataFrame(Elevation)
    df['GPS_time']     = pd.DataFrame(GPS_time)
    df['Surface_idx']  = pd.DataFrame(Surface_idx)
    df['Surface_m']    = pd.DataFrame(Surface_m)
    df['Bed_idx']      = pd.DataFrame(Bed_idx)
    df['Bed_m']        = pd.DataFrame(Bed_m)
    df['pri_pulse']    = pd.DataFrame(pri_pulse)

    df.columns         = ['Traces_Chirp',
                         'Traces_Pulse',
                         'Longitude', 
                         'Latitude',
                         'Elevation',
                         'GPS_time',
                         'Surface_idx',
                         'Surface_m',
                         'Bed_idx',
                         'Bed_m',
                         'pri_pulse']

    pri_chirp_list    = []

    for i in range(len(pri_chirp)):
        closest_pulse_value = find_nearest(pri_pulse, pri_chirp[i])
        pri_chirp_list.append(closest_pulse_value)

    ### dataframe whith pri chirp and closest pri pulse ###

    df_pri                      = pd.DataFrame(pri_chirp)
    df_pri['closest_pri_pulse'] = pd.DataFrame(np.array(pri_chirp_list))
    df_pri.columns              = ['pri_chirp', 'closest_pri_pulse']

    df_comb = df.merge(df_pri, how='inner', left_on='pri_pulse', right_on='closest_pri_pulse')

    traces_chirp_index = np.array(df_comb['Traces_Chirp'] - 1)
    traces_pulse_index = np.array(df_comb['Traces_Pulse'] - 1)

    Data_Chirp = Data_Chirp[:, traces_chirp_index]
    Data_Pulse = Data_Pulse[:, traces_pulse_index]

    data_chirp  = pd.DataFrame(Data_Chirp)
    time        = Time
    longitude   = df_comb['Longitude'].values
    latitude    = df_comb['Latitude'].values
    elevation   = df_comb['Elevation'].values
    gps_time    = df_comb['GPS_time'].values
    surface_idx = df_comb['Surface_idx'].values
    surface_m   = df_comb['Surface_m'].values
    bed_idx     = df_comb['Bed_idx'].values
    bed_m       = df_comb['Bed_m'].values

    return data_chirp, time, longitude, latitude, elevation, gps_time, surface_idx, surface_m, bed_idx, bed_m