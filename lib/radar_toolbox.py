




def radar_twt2surface(data="", twt="", surf_idx=""):

    import numpy  as np
    import time


    # define some short variables
    data     = data                   # Data matrix
    twt      = twt                    # Time array
    surf_idx = surf_idx               # idx of surf. reflection

    # transpose into rows if necessary
    data    = np.transpose(data)

    # get dimension of output matrix
    x_dim   = len(surf_idx)
    y_dim   = len(twt) - surf_idx.min()

    # create empty output matrix 
    out_arr = np.zeros((x_dim, y_dim), dtype=float)

    #print("out_arr shape = {}".format(out_arr.shape))

    # Log every n lines.
    LOG_EVERY_N = 1000

    start = time.time()

    #plt.figure(figsize=(12,7))
    for i in np.arange(0, len(surf_idx)):
        if (i % LOG_EVERY_N) == 0:
            end = time.time()
            print('... processed  {}  of  {}  traces in {:.2f} s'.format(i + 1, len(surf_idx), end - start))

        # get single trace of radargram
        single_trace = data[i][surf_idx[i]::]

        # append as many NaN's as were cut off
        single_trace = np.append(single_trace, np.zeros(surf_idx[i]) + np.nan)

        # cut to minimum extent
        clip     = len(twt) - surf_idx.min()

        # insert single trace in empty out array
        out_arr[i,:] = single_trace[0:clip]

    data_new = np.transpose(out_arr)
    twt_new  = twt[0:clip]

    return data_new, twt_new





def radar_twt2elevation(data='',
                  twt='',
                  surf_idx='',
                  speed_of_ice=1.689e8,
                  dem_surf='',
                  sample_int=""
                  ):

    import time
    import numpy as np
    from scipy import interpolate
    

    speed_of_light = 2.99792458e8
    speed_of_ice   = speed_of_ice

    sample_int = sample_int

    # define some short variables
    data     = data                   # Data matrix
    twt      = twt                    # Time array
    surf_idx = surf_idx               # idx of surf. reflection
    dem_surf = dem_surf

    # transpose into rows
    data    = np.transpose(data)

    # time and depth step
    time_step  = np.diff(twt).mean()
    depth_step = time_step * speed_of_ice / 2

    # maximum height of air and ice to estimate output matrix dims
    max_air_height = int((twt[surf_idx.max()] * speed_of_light) + 1)
    max_ice_height = int(((twt[-1] - twt[surf_idx.min()]) * speed_of_ice / 2) + 1)
    total_height   = max_air_height + max_ice_height + ( dem_surf.max() - dem_surf.min() )

    # get dimension of output matrix
    x_dim   = len(surf_idx)
    y_dim   = int(total_height / sample_int)

    # create empty output matrix 
    out_arr = np.zeros((x_dim, y_dim), dtype=float)

    # elevation axis
    elev_axis = np.round(np.repeat(dem_surf.max(), y_dim) - (np.cumsum(np.repeat(sample_int, y_dim)) - sample_int),1)

    # Log every n lines.
    LOG_EVERY_N = 1000

    start = time.time()

    air_ice_list = []

    #plt.figure(figsize=(12,7))
    for i in np.arange(0, len(surf_idx)):
        if (i % LOG_EVERY_N) == 0:
            end = time.time()
            print('... processed  {}  of  {}  traces in {:.2f} s'.format(i + 1, len(surf_idx), end - start))

        # get single trace of radargram starting at surface reflection
        trace  = data[i][surf_idx[i]::]

        # get index in elevation axis where to put the trace
        tr_idx = (np.abs(elev_axis - dem_surf[i])).argmin()

        # twt array extension factor and length
        fc      = depth_step / sample_int
        new_len = int(fc * len(trace))

        # create ice dB array according to the new matrix dimensions
        x      = np.arange(0, len(trace))
        y      = trace
        f      = interpolate.interp1d(x, y, kind='linear', fill_value="extrapolate")
        x_new  = np.linspace(0, len(trace), new_len)
        y_new  = f(x_new)
        ice_dB = y_new

        # create air dB array
        air_dB = np.repeat(np.nan, tr_idx)

        # create leftover array for below ice
        # print("y_dim = {}, len(ice_dB) = {}, len(air_dB) = {}".format(y_dim, len(ice_dB), len(air_dB)))
        num_nans = y_dim - len(ice_dB) - len(air_dB)
        nan_dB   = np.repeat(np.nan, num_nans)

        # concatenate air, ice and leftover
        full_trace = np.concatenate([air_dB, ice_dB, nan_dB])

        out_arr[i,:] = full_trace

        air_ice_list.append(len(air_dB) + len(ice_dB))

    clip      = max(air_ice_list)
    data_new  = np.transpose(out_arr)[:][0:clip]
    elev_axis = elev_axis[0:clip]

    return data_new, elev_axis, tr_idx







######
def radar_pull2bed(data='', elevation_array='', bed_elevation='', range_resolution_m=''):



    import pandas as pd
    import numpy  as np
    import time


    data = np.array(data.T)

    # create empty numpy arrays for Elevation, dB and Trace number
    all_elevation = np.array([])
    all_dB        = np.array([])
    all_tracenum  = np.array([])

    # first bed elevation, will be our reference
    bed_first      = bed_elevation[0]            
    bed_first_idx  = (np.abs(np.array(elevation_array) - np.array(bed_elevation[0]))).argmin()

    # Log every 500 lines.
    LOG_EVERY_N = 1000

    start = time.time()

    for i in np.arange(0, len(bed_elevation)):

        try:
            if (i % LOG_EVERY_N) == 0:
                end = time.time()
                print('... processed  {}  of  {}  traces in {:.2f} s'.format(i + 1, len(bed_elevation), end - start))


            bed_idx             = (np.abs(np.array(elevation_array) - np.array(bed_elevation)[i])).argmin()
            difference_idx      = bed_first_idx - bed_idx
            new_elevation_array = elevation_array - difference_idx

            # get single trace of radargram
            single_trace = data[i]

            trace_elev = new_elevation_array
            trace_dB   = single_trace
            trace_num  = np.ones(int(single_trace.size)) * i

            all_elevation = np.append(all_elevation, trace_elev, axis=0)
            all_dB        = np.append(all_dB, trace_dB, axis=0)
            all_tracenum  = np.append(all_tracenum, trace_num, axis=0)
            
        except:
            break

        # del single_trace, surf_idx, trace_twt, trace_dB, trace_num, T, T_new

    df_comb             = pd.DataFrame(all_elevation)
    df_comb['dB']       = pd.DataFrame(all_dB)
    df_comb['trace']    = pd.DataFrame(all_tracenum)
    df_comb.columns     = ['elevation', 'dB', 'trace']

    # drop nan's
    df_comb = df_comb.dropna()

    ## create pivot table for heatmap
    df = df_comb.pivot('elevation', 'trace', 'dB')
    df = df.interpolate()
    df = df.iloc[::-1]

    depth_array = np.cumsum(np.repeat(range_resolution_m, len(df.index.values))) - range_resolution_m

    print('==> Done ...')

    return df, depth_array





def add_range_gain(data='', gain_type='', b=2, n=2, f=2):

    import numpy as np
    import pandas as pd

    data         = data
    gain_type    = gain_type
    b            = b
    n            = n
    f            = f
    xlen         = len(data.T)
    ylen         = len(data)
    
    data_gain_matrix = []

    if gain_type == 'linear':
        slope = np.linspace(1, f, ylen)
        
    elif gain_type == 'exponential':
        slope = np.geomspace(1, b**n, ylen)

    for i in np.arange(0, xlen):
        trace_gain  = np.array(data[i]) * slope
        data_gain_matrix.append(trace_gain)
        
    new_data = pd.DataFrame(np.array(data_gain_matrix)).T
        
    return new_data




def automatic_gain_control(data, window=50):

    '''

    '''

    import numpy as np
    import pandas as pd
    # from scipy import ndimage


    # k        = np.array([[a,a,a],[a,a,a],[a,a,a]])
    # new_data = ndimage.convolve(data, k, mode='constant', cval=1.0)

    win  = window
    data = data

    # transpose into rows if necessary
    data    = np.transpose(data)

    new_matrix = []

    LOG_EVERY_N = 1000

    for i in np.arange(0, data.shape[0]):
        if (i % LOG_EVERY_N) == 0:
            print('... processed  {}  of  {}  traces'.format(i + 1, data.shape[0]))
        
        trace = np.array(data[i])
        n     =  len(trace) // win
        rest  = len(trace) - (win * n)

        normalized_trace = []

        
        # for j in range(len(trace) - rest):
        #     ulim, llim = j, j+win
        #     str        = trace[ulim:llim]
        #     norm       = (str-np.min(str))/(np.max(str)-np.min(str))
        #     normalized_trace.append([norm[0]])

        # # last bit
        # ulim, llim = len(trace) - rest, len(trace)
        # subtrace   = trace[ulim:llim]
        # norm       = (str-np.min(str))/(np.max(str)-np.min(str))
        # normalized_trace.append(norm)

        # normalized_trace = np.concatenate(np.array(normalized_trace))

        for j in range(n):
            ulim, llim = j*win, (j+1)*win
            str        = trace[ulim:llim]
            norm       = (str-np.min(str))/(np.max(str)-np.min(str))
            normalized_trace.append(norm)

        # last bit
        ulim, llim = n*win, (n*win) + rest
        str        = trace[ulim:llim]
        norm       = (str-np.min(str))/(np.max(str)-np.min(str))
        normalized_trace.append(norm)

        normalized_trace = np.concatenate(np.array(normalized_trace))

        new_matrix.append(normalized_trace)  
        
    #new_data = pd.DataFrame(new_matrix).T
    new_data = np.transpose(new_matrix)
    print('...  done.')
    
    return new_data



def automatic_gain_control2(data, window=50):

    '''

    '''

    import numpy as np
    import pandas as pd
    # from scipy import ndimage


    # k        = np.array([[a,a,a],[a,a,a],[a,a,a]])
    # new_data = ndimage.convolve(data, k, mode='constant', cval=1.0)

    win  = window
    data = data

    # transpose into rows if necessary
    data    = np.transpose(data)

    # get dimension of output matrix
    x_dim   = data.shape[0]
    y_dim   = data.shape[1]

    # create empty output matrix 
    out_arr = np.zeros((x_dim, y_dim), dtype=float)

    nans       = np.repeat(np.nan, window*2)
    LOG_EVERY_N = 1000

    for i in np.arange(0, data.shape[0]):
        if (i % LOG_EVERY_N) == 0:
            print('... processed  {}  of  {}  traces'.format(i + 1, data.shape[0]))
        
        trace     = np.array(data[i])
        #print(trace.shape)
        new_trace = []
        for j in range(len(trace) - (2*window)):
            k = j + window
            section     = trace[k-window:k+window]
            vmin        = section.min()
            if vmin == np.nan:
                new_value = np.nan
            else:
                vmax        = section.max()
                diff        = vmax - vmin
                value       = vmax - trace[k+window]
                new_value   = (value * 100/diff)*-1

            new_trace.append(new_value)
        new_trace = np.array(new_trace)
        new_trace = np.concatenate([nans, new_trace])
        # new_matrix.append(new_trace) 
        out_arr[i,:] = new_trace 
        
    #new_data = pd.DataFrame(new_matrix).T
    new_data = np.transpose(out_arr)
    print('...  done.')
    
    return new_data
    



########################
########################

def rangegain(self, slope):
    """Apply a range gain.
    Parameters
    ----------
    slope: float
        The slope of the linear range gain to be applied. Maybe try 1.0e-2?
    """
    if isinstance(self.trig, (float, int, np.float, np.int64)):
        gain = self.travel_time[int(self.trig) + 1:] * slope
        self.data[int(self.trig + 1):, :] *= np.atleast_2d(gain).transpose()
    else:
        for i, trig in enumerate(self.trig):
            gain = self.travel_time[int(trig) + 1:] * slope
            self.data[int(trig) + 1:, i] *= gain
    self.flags.rgain = True


def agc(self, window=50, scaling_factor=50):
    """Try to do some automatic gain control
    This is from StoDeep--I'm not sure it is useful but it was easy to roll over so
    I'm going to keep it. I think you should have most of this gone with a bandpass,
    but whatever.
    Parameters
    ----------
    window: int, optional
        The size of window we use in number of samples (default 50)
    scaling_factor: int, optional
        The scaling factor. This gets divided by the max amplitude when we rescale the input.
        Default 50.
    """
    maxamp = np.zeros((self.snum,))
    # In the for loop, old code indexed used range(window // 2). This did not make sense to me.
    for i in range(self.snum):
        maxamp[i] = np.max(np.abs(self.data[max(0, i - window // 2):
                                            min(i + window // 2, self.snum), :]))
    maxamp[maxamp == 0] = 1.0e-6
    self.data *= (scaling_factor / np.atleast_2d(maxamp).transpose()).astype(self.data.dtype)
    self.flags.agc = True


#########################
#########################


def correct4attenuation(data, twt, surf_idx, v_ice=1.68914e8, mode=0, loss_factor=0):

    '''
    modes:  0 = geometric spreading 
            1 = constant ice thickness loss in dB
            2 = both, 1 & 2

    factor: ice thickness dependend loss in dB
    
    '''

    import numpy as np
    import pandas as pd

    data      = np.transpose(data)
    twt       = twt
    surf_idx  = surf_idx
    mode      = mode

    if mode == 0:
        print('==> Correcting for geometrical spreading')

    elif mode == 1:
        print('==> Correcting for geometrical spreading and ice thickness ({} dB/km)'.format(loss_factor))

    if mode == 1:

        # transform dB to amplitude (inverse of: 20 * log10 ?)
        #loss_factor == 10**loss_factor / 20

        # loss factor given in dB/km --> db/m
        loss_factor = loss_factor / 1000

    v_ice     = v_ice / 2
    v_air     = (2.99792458e8 / 2)

    new_matrix = []

    # a little tweak to avoid too small amplitudes in the first array
    twt[0] = twt[1]

    # for every trace
    for trace in range(data.shape[0]):
        
        air_col    = range(int(surf_idx[trace]))
        ice_col    = range(int(surf_idx[trace]), len(data[trace]))
        power      = np.array(data[trace])
        geom_range = []
        att_range  = []
        
        # for every pixel in the air col.
        for pixel in air_col:
            air_m = twt[pixel] * v_air
            geom_range.append(air_m)                    
            if mode == 1 or mode == 2:
                att_range.append(1)
            #print(air_m)
            
        for pixel in ice_col:
            ice_m = twt[pixel] * v_ice
            geom_range.append(ice_m)

            if mode == 1:
                atrange = twt[pixel] * v_ice - twt[int(surf_idx[trace])]
                #print(atrange)
                att_range.append(atrange)
            
        geom_range = np.array(geom_range)
        new_power  = power * (geom_range ** 2)
        #np.savetxt("new_power1.csv", new_power, delimiter=",")
        new_power  = 20 * np.log10(new_power)
        #np.savetxt("geom_range.csv", geom_range, delimiter=",")
        #np.savetxt("att_range.csv", att_range, delimiter=",")
        #np.savetxt("new_power2.csv", new_power, delimiter=",")
        
        if mode == 1:

            att_range = np.array(att_range)
            the_loss  = att_range * loss_factor * 2 
            new_power = new_power + the_loss
            #np.savetxt("the_loss.csv", the_loss, delimiter=",")
            # print(att_range)
            # print(the_loss)
            # print(new_power)
        
        new_matrix.append(new_power)
    
    new_data = np.transpose(np.array(new_matrix))

    return new_data


