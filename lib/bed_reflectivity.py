
def get_bed_reflectivity(crd_object, envelope=200, fixed_envelope_bins=50, n_noise=200, percent=10000, min_width=10, max_width=75, mode='flexible_envelope'):

    import numpy as np
    import pandas as pd
    import copy


    # needs better way to determine noise average
    # only works for continuous bits
    # ==> needs to be divided into continuous sets

    #######################
    # get Bedrock 
    # index, trace, window 
    # and lat lon
    #######################
    data = np.transpose(crd_object.Data)
    lon  = crd_object.Longitude
    lat  = crd_object.Latitude

    bed_trace  = crd_object.Layer['Base']['trace']
    bed_twt    = crd_object.Layer['Base']['value']
    bed_index  = crd_object.Layer['Base']['value_idx']

    surf_twt   = crd_object.Layer['Surface']['value']

    bed_index_list  = []
    bed_trace_list  = []
    bed_win_list    = []
    dB_max_before   = []
    mean_noise_list = []
    bed_twt_list    = []
    surf_twt_list   = []
    lon_list        = []
    lat_list        = []

    # iterate over bedrock picks
    for i in range(len(bed_trace)):
        trace       = bed_trace[i]           # trace number in bedpick
        bed_idx     = bed_index[i]           # time in 

        bed_twt_list.append(bed_twt[i])
        surf_twt_list.append(surf_twt[i])
        lon_list.append(lon[i])
        lat_list.append(lat[i])

        # get index where surface reflection is located
        # idx = (np.abs(twt - bed_twt)).argmin()

        # get single trace of radargram
        trace = data[trace - 1]

        # get average value for noise based on lowest n values
        mean_noise = np.sort(trace)[0:n_noise].mean()

        # delete values around bedrock reflection
        lower_lim = int(bed_idx - envelope)
        upper_lim = int(bed_idx + envelope)


        # reduce data to bed section
        try:
            data_e1   = trace[lower_lim:upper_lim]
            dB_max_bf = data_e1.max()
        except:
            print("")
            print("Problem in line 64: data_e1   = trace[lower_lim:upper_lim]")
            print("Problem in line 65: dB_max_bf = data_e1.max()")
            data_e1   = np.repeat(np.nan, upper_lim-lower_lim)
            dB_max_bf = np.nan

        if 0:
            # find maximum pixel in reduced section
            # and shift around maximum
            max_idx   = data_e1.argmax()
            lower_lim = int(max_idx - int(envelope/2))
            upper_lim = int(max_idx + int(envelope/2))
            data_e2      = data_e1[lower_lim:upper_lim]

        # append data
        bed_index_list.append(bed_idx)
        bed_trace_list.append(bed_trace[i])
        bed_win_list.append(data_e1)
        mean_noise_list.append(mean_noise)
        # Bed_TWT.append(bed_twt)
        # Surf_TWT.append(surf_twt)
        # Longitude.append(lon)
        # Latitude.append(lat)
        dB_max_before.append(dB_max_bf)

    # convert to array or dataframes
    bedrock_index = np.array(bed_index_list)
    bed_traces    = np.array(bed_trace_list)
    bedrock_win   = bed_win_list
    dB_max_before = np.array(dB_max_before)
    mean_noise_floor = np.array(mean_noise_list)

    # create dataframe 
    df_near_bed   = pd.DataFrame(bedrock_win)
    df_no_average = copy.copy(df_near_bed)
    # df            = df_near_bed.rolling(50, center=True, win_type='hamming').mean().T
    df        = df_near_bed

    # print(df_near_bed)
    # print(df)

    # empty lists
    dB          = []
    dB_max_     = []
    x_min_      = []
    x_max_      = []
    width_      = []
    noise_floor = []


    ##################################################
    #
    ##################################################

    if mode == 'fixed_envelope':

        # print(df)
        # print(df.shape)

        df = df.T

        # iterate over columns in dataframe (for real now)
        for col in df:
            # print(df[col])

            # find maximum value in column
            dB_max = df[col].max()
            dB_max_.append(dB_max)
            db_max_idx = df[col].idxmax()

            llim = int( ( len(df[col])/2 ) - fixed_envelope_bins  / 4 )
            ulim = int( ( len(df[col])/2 ) + fixed_envelope_bins )

            try:
                dB_integral   = 20 * np.log10(np.trapz(np.array(df[col][llim:ulim])))
            except:
                dB_integral   = np.nan

            

            dB.append(dB_integral)
            x_min_.append( int(bed_index_list[col] - (fixed_envelope_bins / 4 ) ) )
            x_max_.append( int(bed_index_list[col] + fixed_envelope_bins) )

        # print(x_min_)

    ##################################################
    #
    ##################################################

    if mode == 'flexible_envelope':

        factor  = (percent/100)
        
        df = df.T

        # iterate over columns in dataframe (for real now)
        for i in range(len(df)):

            # find maximum value in column
            dB_max = 20 * np.log10( df.iloc[i].max() )

            # noise average and threshold
            noise_average    = mean_noise_floor[i]
            noise_average_dB = 10 * np.log10(noise_average)
            noise_threshold  = noise_average + np.abs(noise_average * factor)

            llim   = int( ( len(df.iloc[i])/2 ) - fixed_envelope_bins / 20 )            
            ulim   = np.argmax(df.iloc[i][ int( len(df.iloc[i])/2 )::] < noise_threshold)  + int( ( len(df.iloc[i])/2 ) )
            width  = ulim - llim
            
            width = np.diff([llim, ulim])

            # if upper limit is too far above bedrock pick
            if width > max_width:
                ulim = int( len(df.iloc[i])/2 + max_width )
                
            try:
                dB_integral   = 20 * np.log10( np.trapz(np.array(df.iloc[i][llim:ulim])) )

                # if window for integration is too small
                if width < min_width:
                    dB_integral = np.nan

            except:
                dB_integral   = np.nan
                        
            dB.append(dB_integral)
            dB_max_.append(dB_max)
            x_min_.append(llim)
            x_max_.append(ulim)
            noise_floor.append(noise_average)
            width_.append(width)


    dB          = np.array(dB)
    dB_max      = np.array(dB_max_)
    x_min       = np.array(x_min_)
    x_max       = np.array(x_max_)
    noise_floor = np.array(noise_floor)
    noise_mean  = np.nanmean(noise_floor)
    width       = np.array(width_)

    bed_twt   = np.array(bed_twt_list)
    surf_twt  = np.array(surf_twt_list)
    longitude = np.array(lon_list)
    latitude  = np.array(lat_list)

    # Peakiness
    peakyness   = dB_max / dB

    return df, df_no_average, dB, dB_max, dB_max_before, x_min, x_max, peakyness, bed_twt, surf_twt, longitude, latitude, bed_traces
