
def get_layer_backscatter(crd_object, layer_name="", envelope=200, fixed_envelope_bins=50, n_noise=200, percent=10000, min_width=10, max_width=75, mode='flexible_envelope', v_ice=1.69e8):

    import numpy as np
    import pandas as pd
    import copy


    # needs better way to determine noise average
    # only works for continuous bits
    # ==> needs to be divided into continuous sets

    #######################
    # get Layer 
    # index, trace, window 
    # and lat lon
    #######################
    data = np.transpose(crd_object.Data)
    lon  = crd_object.Longitude
    lat  = crd_object.Latitude

    bed_trace  = crd_object.Layer[layer_name]['trace'] - 1 #np.where(crd_object.Layer[layer_name]["value_idx"] != 0)[0] #
    bed_twt    = crd_object.Layer[layer_name]['value'][bed_trace]
    bed_index  = crd_object.Layer[layer_name]['value_idx'][bed_trace]

    surf_twt   = crd_object.Layer['Surface']['value'][bed_trace]

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
        
        if bed_index[i] == 0:
            pass
        else:
            trace       = bed_trace[i]           # trace number in bedpick
            bed_idx     = bed_index[i]       # time in 

            bed_twt_list.append(bed_twt[i])
            surf_twt_list.append(surf_twt[i])
            lon_list.append(lon[i])
            lat_list.append(lat[i])

            # get index where surface reflection is located
            # idx = (np.abs(twt - bed_twt)).argmin()

            # get single trace of radargram
            tr = data[trace]

            # get average value for noise based on lowest n values
            mean_noise = np.sort(tr)[0:n_noise].mean()

            # delete values around bedrock reflection
            lower_lim = int(bed_idx - envelope)
            upper_lim = int(bed_idx + envelope)

            # reduce data to bed section
            data_e1 = tr[lower_lim:upper_lim]                         #np.delete(data[:lower_lim], np.s_[0:upper_lim])

            dB_max_bf = data_e1.max()

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
    df            = df_near_bed#.rolling(50, center=True, win_type='hamming').mean().T
    df        = df.T


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

        # iterate over columns in dataframe (for real now)
        for col in df:

            # find maximum value in column
            dB_max = df[col].max()
            dB_max_.append(dB_max)
            db_max_idx = df[col].idxmax()

            llim = int( ( len(df[col])/2 ) - fixed_envelope_bins / 20 )
            ulim = int( ( len(df[col])/2 ) + fixed_envelope_bins )

            try:
                dB_integral   = np.trapz(np.array(df[col])) 
            except:
                dB_integral   = np.nan

            

            dB.append(dB_integral)
            x_min_.append(llim)
            x_max_.append(ulim)


    ##################################################
    #
    ##################################################

    # if mode == 'flexible_envelope':

    #     factor  = (percent/100)
        
    #     df = df.T

    #     # iterate over columns in dataframe (for real now)
    #     for i in range(len(df)):

    #         # find maximum value in column
    #         dB_max = 20 * np.log10( df.iloc[i].max() )

    #         # noise average and threshold
    #         noise_average    = mean_noise_floor[i]
    #         noise_average_dB = 10 * np.log10(noise_average)
    #         noise_threshold  = noise_average + np.abs(noise_average * factor)

    #         llim   = int( ( len(df.iloc[i])/2 ) - fixed_envelope_bins / 20 )            
    #         ulim   = np.argmax(df.iloc[i][ int( len(df.iloc[i])/2 )::] < noise_threshold)  + int( ( len(df.iloc[i])/2 ) )
    #         width  = ulim - llim
            
    #         width = np.diff([llim, ulim])

    #         # if upper limit is too far above bedrock pick
    #         if width > max_width:
    #             ulim = int( len(df.iloc[i])/2 + max_width )
                
    #         try:
    #             dB_integral   = 20 * np.log10( np.trapz(np.array(df.iloc[i][llim:ulim])) )

    #             # if window for integration is too small
    #             if width < min_width:
    #                 dB_integral = np.nan

    #         except:
    #             dB_integral   = np.nan
                        
    #         dB.append(dB_integral)
    #         dB_max_.append(dB_max)
    #         x_min_.append(llim)
    #         x_max_.append(ulim)
    #         noise_floor.append(noise_average)
    #         width_.append(width)


    # Peakiness
    peakyness   = dB / dB_max

    dB          = 20*np.log10(np.array(dB))
    dB_max      = 20*np.log10(np.array(dB_max_))
    x_min       = np.array(x_min_)
    x_max       = np.array(x_max_)
    noise_floor = np.array(noise_floor)
    noise_mean  = np.nanmean(noise_floor)
    width       = np.array(width_)

    bed_twt   = np.array(bed_twt_list)
    surf_twt  = np.array(surf_twt_list)
    longitude = np.array(lon_list)
    latitude  = np.array(lat_list)

    layer_thickness = ( (bed_twt - surf_twt) / 2 ) * v_ice




    ## Build Dictionary for return

    return_dict = {"DataFrame"        : df,
                    "DF_no_average"   : df_no_average,
                    "dB"              : dB,
                    "dB_max"          : dB_max,
                    "dB_max_before"   : dB_max_before,
                    "x_min"           : x_min,
                    "x_max"           : x_max,
                    "peakyness"       : peakyness,
                    "bed_twt"         : bed_twt,
                    "surf_twt"        : surf_twt,
                    "longitude"       : longitude,
                    "latitude"        : latitude,
                    "bed_traces"      : bed_traces,
                    "layer_thickness" : layer_thickness}

    # return [df, df_no_average, dB, dB_max, dB_max_before, x_min, x_max, peakyness, bed_twt, surf_twt, longitude, latitude, bed_traces]
    return return_dict