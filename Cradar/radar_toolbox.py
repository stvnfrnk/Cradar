







def twt2elevation(data='',
                  twt='',
                  twt_surface='',
                  aircraft_elevation='',
                  speed_of_ice=1.689e8,
                  reference='GPS',
                  DEM_surface='',
                  setting='narrowband',
                  overlap=False,
                  overlap_traces=0,
                  decimate=[True, 2]
                  ):


    import pandas as pd
    import numpy  as np
    import time


    speed_of_light = 2.99792458e8
    speed_of_ice   = speed_of_ice

    data = np.array(data.T)

    # define data frames
    # df       = data
    # df       = df.apply(pd.to_numeric).astype(float)        
    # df       = df.reset_index(drop=True) # reset index
    # df_comb  = pd.DataFrame(columns = ['Elevation', 'dB', 'Trace'])


    # create empty numpy arrays for Elevation, dB and Trace number
    all_elevation = []
    all_dB        = []
    all_tracenum  = []

    # define some short variables
    twt      = twt                 # Time array
    elev     = aircraft_elevation  # Aircraft Elevation
    twt_surf = twt_surface         # twt of surf. reflection


    
    # Log every 500 lines.
    LOG_EVERY_N = 1000

    # create surface and bottom array (m)
    surface_m = []

    start = time.time()

    for i in np.arange(0, elev.size):
        if (i % LOG_EVERY_N) == 0:
            end = time.time()
            print('... processed  {}  of  {}  traces in {:.2f} s'.format(i + 1, elev.size, end - start))

        # get index where surface reflection is located
        if setting == 'emr':
            surf_idx = idx[i]
        else:
            surf_idx = (np.abs(np.array(twt) - np.array(twt_surf)[i])).argmin()

        # get single trace of radargram
        single_trace = data[i]

        # delete values until surface reflection
        single_trace = np.delete(single_trace,np.s_[0:surf_idx])

        #
        T                   = np.delete(twt,np.s_[0:surf_idx],axis=0) # delete traces abofe surface reflection
        T_                  = T - twt[surf_idx] # set start of new twt array to zero
        Depth               = T_ * (speed_of_ice / 2)
        Air_Column          = (twt_surf[i] * speed_of_light) / 2
        Airplane_Elevation  = elev[i]

        if reference == 'GPS':
            surface   = Airplane_Elevation - Air_Column
            Elevation = Airplane_Elevation - Air_Column - Depth
        elif reference == 'DEM':
            surface   = DEM_surface.flatten()
            Elevation = surface[i] - Depth

        surface_m.append(surface)

        trace_elevation     = Elevation
        trace_dB            = single_trace
        trace_num           = np.ones(int(single_trace.size)) * i

        all_elevation.append(trace_elevation)
        all_dB.append(trace_dB)
        all_tracenum.append(trace_num)

        del single_trace, surf_idx, surface, trace_elevation, trace_dB, trace_num, Elevation

    all_elevation = np.concatenate(all_elevation)
    all_dB        = np.concatenate(all_dB)
    all_tracenum  = np.concatenate(all_tracenum)
    
    df_comb             = pd.DataFrame(all_elevation)
    df_comb['dB']       = pd.DataFrame(all_dB)
    df_comb['Trace']    = pd.DataFrame(all_tracenum)
    df_comb.columns     = ['Elevation', 'dB', 'Trace']

    # drop nan's
    df_comb = df_comb.dropna()

    # create final data frame
    # the range interval depends on the frequency range
    if setting == 'wideband': 
        df_comb         = df_comb.round({'Elevation': 1})

        ## create pivot table for heatmap
        df = df_comb.pivot('Elevation', 'Trace', 'dB')
        df = df.interpolate()
        df = df.iloc[::-1]

        # special section for SEGY conversion:
        # traces with a length of more than 32767 samples are not supported
        if len(df) > 32767 :
            print('...(problem) ===> Number of samples: {}'.format(len(df)))
            print('...(problem) ===> Exceeding maximum number of samples for SEGY conversion ( >32767 samples )')
            print('...(action)  ===> Deleting every >>{}<< row'.format(int(decimate[1])))

            # decimate if it is set TRUE
            if decimate[0] == False:
                pass
            elif decimate[0] == True:
                decimate_factor = int(decimate[1])
                df = df.iloc[::decimate_factor, :]

    if setting == 'narrowband' or setting == 'emr':
        df_comb = df_comb.round({'Elevation': 0})

        ## create pivot table for heatmap
        df       = df_comb.pivot('Elevation', 'Trace', 'dB')
        df       = df.interpolate()
        df.index = df.index.astype(int)
        df       = df.iloc[::-1]

    if setting == 'snow':
        df_comb         = df_comb.round({'Elevation': 3})

        ## create pivot table for heatmap
        df = df_comb.pivot('Elevation', 'Trace', 'dB')
        df = df.interpolate()
        df.index = np.around(df.index.values, decimals=2)
        df = df.loc[~df.index.duplicated(keep='first')]
        df = df.iloc[::-1]

    # get Elevation (Z) array
    Z = df.index.values

    if reference == 'DEM':
        surface_m = DEM_surface

    # get the index value of surface reflection
    surface_m = np.array(surface_m).flatten()
    surf_m_idx = np.array([])

    for i in range(len(surface_m)):
        s_idx = (np.abs(df.index.values - np.array(surface_m)[i])).argmin()
        surf_m_idx = np.append(surf_m_idx, s_idx)

    surf_m_idx = surf_m_idx.astype(int)

    # delete crappy traces above surface reflection
    data     = np.array(df.T)
    new_data = []

    for i in range(len(surf_m_idx)):
        trace         = data[i]
        index         = surf_m_idx[i]
        trace[:index] = np.nan
        new_data.append(trace)

    new_data = np.array(new_data).T
    df = pd.DataFrame(new_data)

    if overlap == True:
        df.drop(df.columns[-65:], axis=1, inplace=True)
        df_meta.drop(df_meta.index[-65:], axis=0, inplace=True)

    print('==> Done ...')
    
    return df, Z, surface_m, surf_m_idx







def radar_pull2surface(data='', twt='', twt_surface='', setting=''):

    import pandas as pd
    import numpy  as np
    import time

    data = np.array(data.T)

    # create empty numpy arrays for Elevation, dB and Trace number
    all_twt      = []
    all_dB       = []
    all_tracenum = []

    # define some short variables
    twt      = twt                 # Time array
    twt_surf = twt_surface         # twt of surf. reflection

    # Log every 500 lines.
    LOG_EVERY_N = 1000

    start = time.time()

    for i in np.arange(0, len(twt_surf)):
        if (i % LOG_EVERY_N) == 0:
            end = time.time()
            print('... processed  {}  of  {}  traces in {:.2f} s'.format(i + 1, len(twt_surf), end - start))

        # get index where surface reflection is located
        if setting == 'emr':
            surf_idx = idx[i]
        else:
            surf_idx = (np.abs(np.array(twt) - np.array(twt_surf)[i])).argmin()

        # get single trace of radargram
        single_trace = data[i]

        # delete values until surface reflection
        single_trace = np.delete(single_trace,np.s_[0:surf_idx])

        #
        T                   = np.delete(twt,np.s_[0:surf_idx],axis=0) # delete traces abofe surface reflection
        T_new               = T - twt[surf_idx] # set start of new twt array to zero

        trace_twt           = T_new
        trace_dB            = single_trace
        trace_num           = np.ones(int(single_trace.size)) * i


        all_twt.append(trace_twt)
        all_dB.append(trace_dB)
        all_tracenum.append(trace_num)

    all_twt       = np.concatenate(all_twt)
    all_dB        = np.concatenate(all_dB)
    all_tracenum  = np.concatenate(all_tracenum)

    df_comb             = pd.DataFrame(all_twt)
    df_comb['dB']       = pd.DataFrame(all_dB)
    df_comb['trace']    = pd.DataFrame(all_tracenum)
    df_comb.columns     = ['twt', 'dB', 'trace']

    # drop nan's
    df_comb = df_comb.dropna()

    #df_comb         = df_comb.round({'twt': 9})

    ## create pivot table for heatmap
    df = df_comb.pivot('twt', 'trace', 'dB')
    df = df.interpolate()
    #df = df.iloc[::-1]

    new_twt_array = df.index.values

    print('==> Done ...')
    
    return df, new_twt_array










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






def automatic_gain_control(data, window=''):

    '''

    '''

    import numpy as np
    import pandas as pd

    window = window
    data   = data

    new_matrix = []
    nans       = np.repeat(np.nan, window*2)

    LOG_EVERY_N = 1000

    for i in np.arange(0, len(data.T)):
        if (i % LOG_EVERY_N) == 0:
            print('... processed  {}  of  {}  traces'.format(i + 1, len(data.T)))
        
        trace     = np.array(data[i])
        new_trace = []
        for i in range(len(data) - (2*window)):
            i = i+window
            section     = trace[i-window:i+window]
            vmin        = section.min()
            vmax        = section.max()
            diff        = vmax - vmin
            value       = vmax - trace[i+window]
            new_value   = (value * 100/diff)*-1

            new_trace.append(new_value)
        new_trace = np.array(new_trace)
        new_trace = np.concatenate([nans, new_trace])
        new_matrix.append(new_trace)  
        
    new_data = pd.DataFrame(new_matrix).T
    print('... magig done.')
    
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


    data      = data
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
    for trace in range(data.shape[1]):
        
        air_col    = range(int(surf_idx[trace]))
        ice_col    = range(int(surf_idx[trace]), len(data[trace]))
        power      = np.array(data[trace])
        geom_range = []
        att_range  = []
        
        # for every pixel in the air col.
        for pixel in air_col:
            air_m = twt[pixel] * v_air
            geom_range.append(air_m)                    
            if mode == 1:
                att_range.append(1)
            
        for pixel in ice_col:
            ice_m = twt[pixel] * v_ice
            geom_range.append(ice_m)

            if mode == 1:
                atrange = twt[pixel] * v_ice - twt[int(surf_idx[trace])]
                att_range.append(atrange)
            
        geom_range = np.array(geom_range)
        new_power  = power * (geom_range ** 2)
        new_power  = 20 * np.log10(new_power)
        
        if mode == 1:

            att_range = np.array(att_range)
            the_loss  = att_range * loss_factor * 2 
            new_power = new_power + the_loss
        
        new_matrix.append(new_power)
    
    new_data = pd.DataFrame(np.array(new_matrix)).T

    return new_data






























































































################################
# Pull to surface reflection
# from a cresis radar .mat file
################################

def pull2surface(in_path='', out_path='', file='', region='', overlap=False, setting='narrowband'):

    '''
        in_path         = path with segment folders with twt matfiles
        out_path        = path with segment folders where to save elevation matfiles
        segment         = cresis segment
        
        region          = Antarctica or Greenland, this is relevant for the reprojection
                          of the coordinates (lon, lat to X,Y)
                          EPSG: 3413 for Greenland
                          EPSG: 3031 for Antarctica
                      
        setting         = 'narrowband' or 'wideband'
    
    ''' 


    import h5py                 # for newer .mat (usually cresis matfiles)
    import scipy.io             # for older .mat (usually self-output matfiles)
    import numpy as np
    import pandas as pd
    import glob, os
    import pyproj
    import geopy.distance


    speed_of_light = 2.99792458e8

    #print('===> Processing File: __ {}'.format(file))
    print('')

    if file == '':
        print('Specify a file')
        exit()
    else:
        print('Filename: {}'.format(file))

    if region == '':
        print('Specify a region (Greenland or Antarctica)')
        exit()
    else:
        print('Region: {}'.format(region))

    if in_path == '':
        in_path = os.getcwd()
        print('No in_path defined, setting cwd as in_path')

    if out_path == '':
        out_path = os.getcwd()
        print('No out_path defined, setting cwd as out_path')


    if not os.path.exists(out_path):
            os.mkdir(out_path)

    #for file in sorted(glob.glob(input_file)):

    # don't process Data_img_... files      
    if 'img' not in file:
        # process only not already converted files
        if not file.endswith('elevation.mat'):
            print('')
            print('Calculating true elevation for: {}'.format(file))


            # depending on the .mat file version
            # either scipy.io (older versions) or h5py (newer versions)
            # will be used to load the file        
            try:
                mat     = scipy.io.loadmat(file)
                reader  = 'scipy'
            except NotImplementedError:
                mat     = h5py.File(file, mode='r')
                reader  = 'h5py'


            # loading with scipy.io
            if reader == 'scipy':   
                print('Using scipy.io to load matfile')

                df_meta               = pd.DataFrame(np.array(mat['GPS_time'])).T                 
                df_meta['Longitude']  = pd.DataFrame(np.array(mat['Longitude'])).T   
                df_meta['Latitude']   = pd.DataFrame(np.array(mat['Latitude'])).T     
                df_meta['Elevation']  = pd.DataFrame(np.array(mat['Elevation'])).T   
                df_meta['Filename']   = file
                df_meta.index.name    = 'index'
                df_meta.columns       = ['GPS_time', 'Longitude', 'Latitude', \
                                         'Aircraft_Elevation', 'Filename']
            
                df       = df.reset_index(drop=True) # reset index
                surf    = np.array(mat['Surface']).T # array with the twt of the surface reflection

                if mat['Time'].shape[0] == 1:
                    twt     = np.array(mat['Time']).T # array with the twt (y-axis)

                if mat['Time'].shape[0] > 1:
                    twt     = np.array(mat['Time']) # array with the twt (y-axis)

                else:
                    print('problem...')
                    pass

                elev    = np.array(mat['Elevation']).T # array with the aircraft elevation
                #bott    = np.array(mat['Bottom']).T

            #loading with h5py
            if reader == 'h5py':
                print('Using h5py to load matfile')

                df_meta = pd.DataFrame(np.array(mat['GPS_time']))                
                df_meta['Longitude'] = pd.DataFrame(np.array(mat['Longitude']))   
                df_meta['Latitude'] = pd.DataFrame(np.array(mat['Latitude']))     
                df_meta['Elevation'] = pd.DataFrame(np.array(mat['Elevation']))  
                df_meta['Filename'] = file

                df_meta.index.name = 'index'
                df_meta.columns = ['GPS_time', 'Longitude', 'Latitude', \
                                   'Aircraft_Elevation', 'Filename']
                df      = pd.DataFrame(np.log10(np.array(mat['Data']))).T
                surf    = np.array(mat['Surface'])
                twt     = np.array(mat['Time']).T
                elev    = np.array(mat['Elevation'])
                #bott    = np.array(mat['Bottom'])

                
            df       = df.apply(pd.to_numeric).astype(float)        
            df       = df.reset_index(drop=True) # reset index
            df_comb  = pd.DataFrame(columns = ['twt', 'dB', 'Trace'])

            # Log every 500 lines.
            LOG_EVERY_N = 500

            for i in np.arange(0, elev.size):
                if (i % LOG_EVERY_N) == 0:
                    print('===> Processed  {}  of  {}  Traces'.format(i + 1, elev.size))

                # get index where surface reflection is located
                idx = (np.abs(np.array(twt) - np.array(surf)[i])).argmin()

                # get single trace of radargram
                data = np.array(df[i])

                # delete values until surface reflection
                data = np.delete(data,np.s_[1:idx])

                # delete twt values until surface reflection
                TWT = np.delete(twt,np.s_[0:idx],axis=0)
                
                # subtract twt value from surface reflection from all values
                # to start with twt zero in the array
                TWT = TWT - twt[idx]
                
                trace = np.ones(data.size) * i

                df_trace            = pd.DataFrame(TWT)
                df_trace['dB']      = pd.DataFrame(data)
                df_trace['Trace']   = pd.DataFrame(trace)
                df_trace.columns    = ['twt', 'dB', 'Trace']
                df_comb             = df_comb.append(df_trace)
            
            df = df_comb.pivot('twt', 'Trace', 'dB')
            df = df.interpolate()
                
                
            if overlap == True:
                df.drop(df.columns[-65:], axis=1, inplace=True)
                df_meta.drop(df_meta.index[-65:], axis=0, inplace=True)


            ## Get real distance for traces 
            spacing = np.array([])
            for i in range(1, len(df_meta)-1):
                    coord_1    = (df_meta['Latitude'][i], df_meta['Longitude'][i])
                    coord_2    = (df_meta['Latitude'][i + 1], df_meta['Longitude'][i + 1])
                    f          = geopy.distance.geodesic(coord_1, coord_2).meters
                    spacing    = np.append(spacing, f) 

            distance       = np.cumsum(spacing[0:-1])
            distance       = np.insert(distance, 0, 0)
            distance       = np.insert(distance, len(distance), distance[-1] + spacing.mean())
            distance       = np.insert(distance, len(distance), distance[-1] + spacing.mean())
            spacing        = np.insert(spacing, 0, spacing.mean())
            spacing        = np.insert(spacing, len(spacing), spacing.mean())


            if region == 'Greenland':
                ## Project Lat, Lon to X, Y in EPSG:3413
                EPSG=pyproj.Proj("EPSG:3413")
                df_meta['X'], df_meta['Y'] = EPSG(np.array(df_meta['Longitude']), \
                                            np.array(df_meta['Latitude']))

            if region == 'Antarctica':       
                ## Project Lat, Lon to X, Y in EPSG:3413
                EPSG=pyproj.Proj("EPSG:3031")
                df_meta['X'], df_meta['Y'] = EPSG(np.array(df_meta['Longitude']), \
                                            np.array(df_meta['Latitude']))

            if not os.path.exists(out_path):
                os.makedirs(out_path)

            print('===> Done Calculating...')

            #####################
            # Save as .mat file
            #####################

            # test if the following variables exist
            # older files do not contain this information

            try:
                pitch = mat['Pitch'][0]
            except:
                pitch = 'empty'
                pass

            try:
                roll = mat['Roll'][0]
            except:
                roll = 'empty'
                pass

            try:
                heading = mat['Heading'][0]
            except:
                heading = 'empty'
                pass                


            full_dict = {'Data'                 : df.values,
                         'Time'                 : df.index.values,
                         'GPS_time'             : df_meta['GPS_time'].values,
                         'Latitude'             : df_meta['Latitude'].values, 
                         'Longitude'            : df_meta['Longitude'].values, 
                         'Aircraft_Elevation'   : df_meta['Aircraft_Elevation'].values,
                         'Spacing'              : spacing,
                         'Distance'             : distance,
                         'X'                    : df_meta['X'].values,
                         'Y'                    : df_meta['Y'].values,
                         'Pitch'                : pitch,
                         'Roll'                 : roll,
                         'Heading'              : heading
                         #'Bottom'               : bottom_m,
                         #'Surface'              : surface_m
                         }

            out_filename = file.split('/')[-1].split('.mat')[0] + '_surface.mat'

            scipy.io.savemat(out_path + '/' + out_filename, full_dict)
            print('===> Saved Frame as: {}'.format(out_filename))

    print('===> DONE !!')



################################################################################################
################################################################################################
################################################################################################


################################
# Calculate True Elevation
# from a cresis radar .mat file
################################

def calc_elevation(in_path='', out_path='', file='', region='', speed_of_ice=1.689e8, overlap=False, setting='', reference='DEM', geotif = '', number_of_gaps=100, decimate=[False, 2]):

    '''
        in_path         = path with segment folders with twt matfiles
        out_path        = path with segment folders where to save elevation matfiles
        segment         = cresis segment
        
        region          = Antarctica or Greenland, this is relevant for the reprojection
                          of the coordinates (lon, lat to X,Y)
                          EPSG: 3413 for Greenland
                          EPSG: 3031 for Antarctica
                      
        speed_of_ice    = Usually 1.689 for e=3.15, but it can be changed

        setting         = 'narrowband', 'wideband' for rds or 'snow' for uwbm-snowradar

        reference       = 'reflection' or 'Laserscanner'
    
    ''' 
                
    import h5py                 # for newer .mat (usually cresis matfiles)
    import scipy.io             # for older .mat (usually self-output matfiles)
    import numpy as np
    import pandas as pd
    import glob, os
    import pyproj
    import geopy.distance
    
    
    speed_of_light = 2.99792458e8

    #print('===> Processing File: __ {}'.format(file))
    print('')
    
    if file == '':
        print('Specify a file')
        exit()
    else:
        print('Filename: {}'.format(file))
    
    if region == '':
        print('Specify a region (Greenland or Antarctica)')
        exit()
    else:
        print('Region: {}'.format(region))
    
    if in_path == '':
        in_path = os.getcwd()
        print('No in_path defined, setting cwd as in_path')
    
    if out_path == '':
        out_path = os.getcwd()
        print('No out_path defined, setting cwd as out_path')
    
    
    if not os.path.exists(out_path):
            os.mkdir(out_path)
    
    
    ####################################################      
    # check number of DEM gaps (ALS_nans)
    
    try:
        mat     = scipy.io.loadmat(file)
    except NotImplementedError:
        mat     = h5py.File(file, mode='r')
    
    if reference == 'DEM': 

        from geo_toolbox import extract_geotif_values

        if region == 'Greenland':
            EPSG = 3413
        if region == 'Antarctica':       
            EPSG = 3031

        try:
            mat     = scipy.io.loadmat(file)
        except NotImplementedError:
            mat     = h5py.File(file, mode='r')

        df_tmp             = pd.DataFrame(np.array(mat['Longitude']))#.T   
        df_tmp['Latitude'] = pd.DataFrame(np.array(mat['Latitude']))#.T
        df_tmp.columns     = ['Longitude', 'Latitude']

        df_tmp['DEM'] = extract_geotif_values(geotif, df_tmp, EPSG=EPSG)

        try:
            DEM_nans     = int(df_tmp['DEM'].isna().sum())
        except:
            pass
    
        if DEM_nans > number_of_gaps:
            reference = 'reflection'
            print('===> More than {} gaps in DEM Data...'.format(number_of_gaps))
            print('===> Using surface reflection instead.')
        else: 
            pass
    else:
        pass
    
    ####################################################
    
    #for file in sorted(glob.glob(input_file)):
    
    # don't process Data_img_... files      
    if 'img' not in file:
        # process only not already converted files
        if not file.endswith('elevation.mat'):
            print('')
            print('Calculating true elevation for: {}'.format(file))
    
    
            # depending on the .mat file version
            # either scipy.io (older versions) or h5py (newer versions)
            # will be used to load the file        
            try:
                mat     = scipy.io.loadmat(file)
                reader  = 'scipy'
            except NotImplementedError:
                mat     = h5py.File(file, mode='r')
                reader  = 'h5py'
    
    
            # loading with scipy.io
            if reader == 'scipy':   
                print('Using scipy.io to load matfile')
                ######################### EMR ####################################
                if setting == 'emr':
                    print('Dealing with AWI EMR data...')
                    
                    df_meta              = pd.DataFrame(np.array(mat['Longitude'])).T   
                    df_meta['Latitude']  = pd.DataFrame(np.array(mat['Latitude'])).T     
                    df_meta['Elevation'] = pd.DataFrame(np.array(mat['Altitude'])).T  
                    df_meta['Filename']  = file
                    df_meta['GPS_time']  = np.ones(len(df_meta['Latitude'])) * -9999
                    df_meta['DEM']       = df_tmp['DEM']#pd.DataFrame(np.array(mat['DEM'])).T
    
                    df_meta.index.name = 'index'
                    df_meta.columns = ['Longitude', 'Latitude', 'Aircraft_Elevation', 'Filename', 'GPS_time', 'DEM']
                    df      = pd.DataFrame(np.log10(np.array(mat['Data'])))
    
                    twt     = np.array(mat['Time'][0]).T
                    elev    = np.array(mat['Altitude'][0])
                    surf    = np.ones(len(elev)) * -9999
                    bott    = np.ones(len(elev)) * -9999
                    idx     = np.array(mat['Surface_idx'][0])
                    
                
                else:
    
                    df_meta               = pd.DataFrame(np.array(mat['GPS_time'])).T                 
                    df_meta['Longitude']  = pd.DataFrame(np.array(mat['Longitude'])).T   
                    df_meta['Latitude']   = pd.DataFrame(np.array(mat['Latitude'])).T     
                    df_meta['Elevation']  = pd.DataFrame(np.array(mat['Elevation'])).T   
                    df_meta['Filename']   = file
                    df_meta.index.name    = 'index'
                    df_meta.columns       = ['GPS_time', 'Longitude', 'Latitude', \
                                             'Aircraft_Elevation', 'Filename']
    
                    df      = pd.DataFrame(np.log10(np.array(mat['Data']))) # radar matrix
    
                    # chech how the ice surface boundary is defined
                    # surface reflection or Laserscanner data?
                    if reference == 'reflection':
                        surf    = np.array(mat['Surface']).T
                        print('==> Using surface reflection as ice-surface boundary.')
                    elif reference == 'DEM':
                        #print('==> Trying to use DEM data as ice-surface boundary.')
                        try:
                            df_meta['DEM'] = df_tmp['DEM']
                        except:
                            print('No DEM data found')
                        #surf    = (np.array(df_meta['Aircraft_Elevation']) - np.array(df_meta['DEM'])) / 2.99792458e8 * 2
                        surf    = np.array(mat['Surface']).T
                        print('==> Using DEM data as ice-surface boundary.')
    
    
                    if mat['Time'].shape[0] == 1:
                        twt     = np.array(mat['Time']).T # array with the twt (y-axis)
    
                    elif mat['Time'].shape[0] > 1:
                        twt     = np.array(mat['Time']) # array with the twt (y-axis)
    
                    else:
                        print('...(problem) ==> Time array shape smaller than 1')
                        pass
    
                    elev    = np.array(mat['Elevation']).T # array with the aircraft elevation
                    bott    = np.array(mat['Bottom']).T
    
            #loading with h5py
            elif reader == 'h5py':
                print('Using h5py to load matfile')
    
                df_meta = pd.DataFrame(np.array(mat['GPS_time']))                
                df_meta['Longitude'] = pd.DataFrame(np.array(mat['Longitude']))   
                df_meta['Latitude'] = pd.DataFrame(np.array(mat['Latitude']))     
                df_meta['Elevation'] = pd.DataFrame(np.array(mat['Elevation']))  
                df_meta['Filename'] = file
    
                df_meta.index.name = 'index'
                df_meta.columns = ['GPS_time', 'Longitude', 'Latitude', \
                                   'Aircraft_Elevation', 'Filename']
                df      = pd.DataFrame(np.log10(np.array(mat['Data']))).T
    
    
                # chech how the ice surface boundary is defined
                # surface reflection or Laserscanner data?
                if reference == 'reflection':
                    surf    = np.array(mat['Surface'])
                    print('==> Using surface reflection as ice-surface boundary.')
                elif reference == 'DEM':
                    try:
                        df_meta['DEM'] = df_tmp['DEM']
                    except:
                        print('No DEM data found')
    
                    #surf    = (df_meta['Aircraft_Elevation'] - df_meta['ALS']) / 2.99792458e8 * 2
                    surf    = np.array(mat['Surface'])
                    print('==> Using DEM data as ice-surface boundary.')
    
    
                twt     = np.array(mat['Time']).T
                elev    = np.array(mat['Elevation'])
                bott    = np.array(mat['Bottom'])
    
    
            
    
    
    
            df       = df.apply(pd.to_numeric).astype(float)        
            df       = df.reset_index(drop=True) # reset index
            df_comb  = pd.DataFrame(columns = ['Elevation', 'dB', 'Trace'])
    
            # Log every 500 lines.
            LOG_EVERY_N = 500
    
            # create surface and bottom array (m)
            surface_m = []
            bottom_m  = []
    
            for i in np.arange(0, elev.size):
                if (i % LOG_EVERY_N) == 0:
                    print('===> Processed  {}  of  {}  Traces'.format(i + 1, elev.size))
    
    
                # get index where surface reflection is located
                if setting == 'emr':
                    surf_idx = idx[i]
                else:
                    surf_idx = (np.abs(np.array(twt) - np.array(surf)[i])).argmin()
    
                # get single trace of radargram
                data = np.array(df[i])
    
                # delete values until surface reflection
                data = np.delete(data,np.s_[0:surf_idx])
    
                #
                T  = np.delete(twt,np.s_[0:surf_idx],axis=0) # delete traces abofe surface reflection
                T_ = T - twt[surf_idx] # set start of new twt array to zero
    
                Depth               = T_ * (speed_of_ice / 2)
    
                Air_Column          = (surf[i] * speed_of_light) / 2
                Airplane_Elevation  = elev[i]
    
                if reference == 'reflection':
                    surface   = Airplane_Elevation - Air_Column
                    Elevation = Airplane_Elevation - Air_Column - Depth
                elif reference == 'DEM':
                    surface   = df_meta['DEM'][i]
                    Elevation = df_meta['DEM'][i] - Depth
    
                surface_m.append(surface)
    
                bottom = Airplane_Elevation - (bott[i] * speed_of_ice / 2)
                bottom_m.append(bottom)
    
    
    
                trace = np.ones(data.size) * i
    
                df_trace            = pd.DataFrame(Elevation)
                df_trace['dB']      = pd.DataFrame(data)
                df_trace['Trace']   = pd.DataFrame(trace)
                df_trace.columns    = ['Elevation', 'dB', 'Trace']
    
    
    
                df_comb             = df_comb.append(df_trace)
                #df_comb             = df_comb.round({'ElevationWGS84': 0})
    
            df_comb = df_comb.dropna()
    
            if setting == 'wideband': 
                df_comb         = df_comb.round({'Elevation': 1})
                
                ## create pivot table for heatmap
                df = df_comb.pivot('Elevation', 'Trace', 'dB')
                df = df.interpolate()
                df = df.iloc[::-1]
    
    
                if len(df) > 32767 :
                    print('...(problem) ===> Number of samples: {}'.format(len(df)))
                    print('...(problem) ===> Exceeding maximum number of samples for SEGY conversion ( >32767 samples )')
                    print('...(action)  ===> Deleting every >>{}<< row'.format(int(decimate[1])))
    
                    if decimate[0] == False:
                        pass
                    elif decimate[0] == True:
    
                        decimate_factor = int(decimate[1])
                        df = df.iloc[::decimate_factor, :]
    
            if setting == 'narrowband' or setting == 'emr':
                df_comb = df_comb.dropna()
                df_comb = df_comb.round({'Elevation': 0})
    
                ## create pivot table for heatmap
                df = df_comb.pivot('Elevation', 'Trace', 'dB')
                df = df.interpolate()
                df.index = df.index.astype(int)
                df = df.iloc[::-1]
    
    
            if setting == 'snow':
                df_comb         = df_comb.round({'Elevation': 3})
                
                ## create pivot table for heatmap
                df = df_comb.pivot('Elevation', 'Trace', 'dB')
                df = df.interpolate()
                df.index = np.around(df.index.values, decimals=2)
                df = df.loc[~df.index.duplicated(keep='first')]
                df = df.iloc[::-1]
                #df.index = (df.index * 10).astype(int)
    
    
            surface_m = np.array(surface_m)
            bottom_m  = np.array(bottom_m)
    
            if overlap == True:
                df.drop(df.columns[-65:], axis=1, inplace=True)
                df_meta.drop(df_meta.index[-65:], axis=0, inplace=True)
    
    
            ## Get real distance for traces 
            spacing = np.array([])
            for i in range(1, len(df_meta)-1):
                    coord_1    = (df_meta['Latitude'][i], df_meta['Longitude'][i])
                    coord_2    = (df_meta['Latitude'][i + 1], df_meta['Longitude'][i + 1])
                    f          = geopy.distance.geodesic(coord_1, coord_2).meters
                    spacing    = np.append(spacing, f) 
    
            distance       = np.cumsum(spacing[0:-1])
            distance       = np.insert(distance, 0, 0)
            distance       = np.insert(distance, len(distance), distance[-1] + spacing.mean())
            distance       = np.insert(distance, len(distance), distance[-1] + spacing.mean())
            spacing        = np.insert(spacing, 0, spacing.mean())
            spacing        = np.insert(spacing, len(spacing), spacing.mean())
    
    
            if region == 'Greenland':
                ## Project Lat, Lon to X, Y in EPSG:3413
                EPSG=pyproj.Proj("EPSG:3413")
                df_meta['X'], df_meta['Y'] = EPSG(np.array(df_meta['Longitude']), \
                                            np.array(df_meta['Latitude']))
    
            if region == 'Antarctica':       
                ## Project Lat, Lon to X, Y in EPSG:3413
                EPSG=pyproj.Proj("EPSG:3031")
                df_meta['X'], df_meta['Y'] = EPSG(np.array(df_meta['Longitude']), \
                                            np.array(df_meta['Latitude']))
    
            if not os.path.exists(out_path):
                os.makedirs(out_path)
    
            print('===> Done Calculating...')
    
            #####################
            # Save as .mat file
            #####################
    
            # test if the following variables exist
            # older files do not contain this information
    
            try:
                pitch = mat['Pitch']
            except:
                pitch = 'empty'
                pass
    
            try:
                roll = mat['Roll']
            except:
                roll = 'empty'
                pass
    
            try:
                heading = mat['Heading']
            except:
                heading = 'empty'
                pass   
    
    
    
    
    
            if reference == 'DEM':
                full_dict = {'Data'                 : df.values,
                             'Elevation'            : df.index.values,
                             'GPS_time'             : df_meta['GPS_time'].values,
                             'Latitude'             : df_meta['Latitude'].values, 
                             'Longitude'            : df_meta['Longitude'].values, 
                             #'Aircraft_Elevation'   : df_meta['Aircraft_Elevation'].values,
                             'Spacing'              : spacing,
                             'Distance'             : distance,
                             'X'                    : df_meta['X'].values,
                             'Y'                    : df_meta['Y'].values,
                             'Pitch'                : pitch,
                             'Roll'                 : roll,
                             'Heading'              : heading,
                             'Bottom'               : bottom_m,
                             'Surface'              : surface_m,
                             'DEM'                  : df_meta['DEM'].values,
                             'DEM_nans'             : DEM_nans
                             }
           
            else:
                full_dict = {'Data'                 : df.values,
                             'Elevation'            : df.index.values,
                             'GPS_time'             : df_meta['GPS_time'].values,
                             'Latitude'             : df_meta['Latitude'].values, 
                             'Longitude'            : df_meta['Longitude'].values, 
                             #'Aircraft_Elevation'   : df_meta['Aircraft_Elevation'].values,
                             'Spacing'              : spacing,
                             'Distance'             : distance,
                             'X'                    : df_meta['X'].values,
                             'Y'                    : df_meta['Y'].values,
                             'Pitch'                : pitch,
                             'Roll'                 : roll,
                             'Heading'              : heading,
                             'Bottom'               : bottom_m,
                             'Surface'              : surface_m,
                             }
                
            if setting == 'emr':
                full_dict = {'Data'                 : df.values,
                             'Elevation'            : df.index.values,
                             'Latitude'             : df_meta['Latitude'].values, 
                             'Longitude'            : df_meta['Longitude'].values,
                             'X'                    : df_meta['X'].values,
                             'Y'                    : df_meta['Y'].values,
                             'Spacing'              : spacing,
                             'Distance'             : distance,
                             'DEM'                  : df_meta['DEM'].values,
                              }
    
    
            suffix = ''
            if reference == 'DEM':
                suffix = 'DEM'
            elif reference == 'reflection':
                suffix = 'reflection'
    
    
            out_filename = file.split('/')[-1].split('.mat')[0] + '_elevation_' + suffix + '.mat'
    
            scipy.io.savemat(out_path + '/' + out_filename, full_dict)
            print('===> Saved Frame as: {}'.format(out_filename))
    
    print('===> DONE !!')




################################
# Connect single frames to a
# larger piece
################################


def combine_frames(frame_list='', output_filename='', z_mode='elevation', overlap=False, overlap_traces=0, Time_Array=False):

    '''
    Reads a list of frames and connects them
    The list should have the following format:

        frames = ['Data_...._001_.mat', 'Data_...._002_.mat', 'Data_...._003_.mat', etc...]

    '''
        
    import scipy.io
    import h5py
    import pandas as pd
    import numpy as np
    import sys
    import glob
    import os


    # CHECKS

    if output_filename == '':
        print('you did not provide a output filename')
        print('setting output filename to output.mat')
        output_filename = 'output.mat'

    ##################################################
    ##################################################
    # ELEVATION
    ##################################################
    ##################################################

    if z_mode == 'elevation':

        Data_               = []
        Elevation_          = []
        GPS_time_           = []
        Latitude_           = []
        Longitude_          = []
        X_                  = []
        Y_                  = []
        Aircraft_Elevation_ = []
        Spacing_            = []
        Pitch_              = []
        Roll_               = []
        Heading_            = []
        Bottom_             = []
        Surface_            = []


        for file in frame_list:

            mat                = scipy.io.loadmat(file)

            # for older frames with overlap
            # cuts off 66 traces
            if overlap == True:

                Data               = pd.DataFrame(np.array(mat['Data']))
                Elevation          = np.array(mat['Elevation'])[0]

                ol                 = Data.shape[1] - overlap_traces
                
                #print(ol)
                #print(len(np.array(mat['Latitude'])))

                GPS_time           = np.array(mat['GPS_time'])[0][0:ol]
                Latitude           = np.array(mat['Latitude'])[0][0:ol]
                Longitude          = np.array(mat['Longitude'])[0][0:ol]
                X                  = np.array(mat['X'])[0][0:ol]
                Y                  = np.array(mat['Y'])[0][0:ol]
                #Aircraft_Elevation = np.array(mat['Aircraft_Elevation'])[0][0:ol]
                Spacing            = np.array(mat['Spacing'])[0][0:ol]
                Pitch              = np.array(mat['Pitch'])[0][0:ol]
                Roll               = np.array(mat['Roll'])[0][0:ol]
                Heading            = np.array(mat['Heading'])[0][0:ol]
                Bottom             = np.array(mat['Bottom'])[0][0:ol]
                Surface            = np.array(mat['Surface'])[0][0:ol]

                Data.index         = Elevation
                Data               = Data[Data.columns[0:ol]]

            if overlap == False:

                Data               = pd.DataFrame(np.array(mat['Data']))
                Elevation          = np.array(mat['Elevation'])[0]
                GPS_time           = np.array(mat['GPS_time'])[0]
                Latitude           = np.array(mat['Latitude'])[0]
                Longitude          = np.array(mat['Longitude'])[0]
                X                  = np.array(mat['X'])[0]
                Y                  = np.array(mat['Y'])[0]
               # Aircraft_Elevation = np.array(mat['Aircraft_Elevation'])[0]
                Spacing            = np.array(mat['Spacing'])[0]
                Pitch              = np.array(mat['Pitch'])[0]
                Roll               = np.array(mat['Roll'])[0]
                Heading            = np.array(mat['Heading'])[0]
                try:
                    Bottom             = np.array(mat['Bottom'])[0]
                    Surface            = np.array(mat['Surface'])[0]
                except:
                    Bottom = 'empty'
                    Surface = 'empty'

                Data.index         = Elevation

            # check if navigation data is available
            # it might not be for files older than 2012
            if Pitch == '':
                navigation = 'off'
            else:
                navigation = 'off'

            Data_.append(Data)
            Elevation_.append(Elevation)
            GPS_time_.append(GPS_time)
            Latitude_.append(Latitude)
            Longitude_.append(Longitude)
            X_.append(X)
            Y_.append(Y)
            #Aircraft_Elevation_.append(Aircraft_Elevation)
            Spacing_.append(Spacing)
            Pitch_.append(Pitch)
            Roll_.append(Roll)
            Heading_.append(Heading)
            Bottom_.append(Bottom)
            Surface_.append(Surface)


        Data               = pd.concat(Data_, axis=1, ignore_index=True)
        Data               = Data[::-1]
        Elevation          = np.array(Data.index.astype(int))#np.concatenate(Elevation_WGS84_)
        GPS_time           = np.concatenate(GPS_time_)
        Latitude           = np.concatenate(Latitude_)
        Longitude          = np.concatenate(Longitude_)
        X                  = np.concatenate(X_)
        Y                  = np.concatenate(Y_)
        #Aircraft_Elevation = np.concatenate(Aircraft_Elevation_)
        Spacing            = np.concatenate(Spacing_)
        try:
            Bottom             = np.concatenate(Bottom_)
            Surface            = np.concatenate(Surface_)
        except:
            Bottom = 'empty'
            Surface = 'empty'


        # check if navigation data is available
        if navigation == 'off':

            Pitch   = 'empty'
            Roll    = 'empty'
            Heading = 'empty'

        else:

            Pitch              = np.concatenate(Pitch_)
            Roll               = np.concatenate(Roll_)
            Heading            = np.concatenate(Heading_)



        Distance = np.cumsum(Spacing)


        full_dict = {'Data'                 : Data.values,
                     'Elevation'            : Elevation,
                     'GPS_time'             : GPS_time,
                     'Latitude'             : Latitude, 
                     'Longitude'            : Longitude,
                     'X'                    : X,
                     'Y'                    : Y,
                     #'Aircraft_Elevation'   : Aircraft_Elevation,
                     'Spacing'              : Spacing,
                     'Distance'             : Distance,
                     'Pitch'                : Pitch,
                     'Roll'                 : Roll,
                     'Heading'              : Heading,
                     'Bottom'               : Bottom,
                     'Surface'              : Surface
                     }

        scipy.io.savemat(output_filename, full_dict)
        print('===> Saved Frame as: {}'.format(output_filename))

    
    ##################################################
    ##################################################
    # TWT
    ##################################################
    ##################################################

    navigation = 'off'

    if z_mode == 'twt':

        Data_               = []
        GPS_time_           = []
        Latitude_           = []
        Longitude_          = []
        Elevation_          = [] # aircraft elevation
        Pitch_              = []
        Roll_               = []
        Heading_            = []
        Bottom_             = []
        Surface_            = []
        Shot_ID_            = []

        
        for file in frame_list:
    
            mat     = h5py.File(file, mode='r')

            # for older frames with overlap
            # cuts off 66 traces
            if overlap == True:

                Data               = pd.DataFrame(np.array(mat['Data'])).T
                Time               = np.array(mat['Time'])[0]

                ol                 = Data.shape[1] - overlap_traces
                
                Elevation          = np.array(mat['Elevation'])[0][0:ol] # aircraft elevation

                GPS_time           = np.array(mat['GPS_time'])[0][0:ol]
                Latitude           = np.array(mat['Latitude'])[0][0:ol]
                Longitude          = np.array(mat['Longitude'])[0][0:ol]
                Pitch              = np.array(mat['Pitch'])[0][0:ol]
                Roll               = np.array(mat['Roll'])[0][0:ol]
                Heading            = np.array(mat['Heading'])[0][0:ol]
                try:
                    Bottom             = np.array(mat['Bottom'])[0][0:ol]
                    Surface            = np.array(mat['Surface'])[0][0:ol]
                except:
                    Bottom = 'empty'
                    Surface = 'empty'

                Data.index         = Time
                Data               = Data[Data.columns[0:ol]]

            if overlap == False:

                Data               = pd.DataFrame(np.array(mat['Data'])).T
                Time               = np.array(mat['Time'])[0]

                ol                 = Data.shape[1] - overlap_traces
                
                Elevation          = np.array(mat['Elevation']) # aircraft elevation

                GPS_time           = np.array(mat['GPS_time'])
                Latitude           = np.array(mat['Latitude'])
                Longitude          = np.array(mat['Longitude'])
                Pitch              = np.array(mat['Pitch'])
                Roll               = np.array(mat['Roll'])
                Heading            = np.array(mat['Heading'])
                try:
                    Bottom             = np.array(mat['Bottom'])
                    Surface            = np.array(mat['Surface'])
                except:
                    Bottom = 'empty'
                    Surface = 'empty'

                Data.index         = Time


            # check if navigation data is available
            # it might not be for files older than 2012
            if Pitch == '':
                navigation = 'off'
            else:
                navigation = 'off'


            # create shot ID
            frame   = file.split('.')[0].split('Data_')[1]
            shots   = np.arange(1, len(Longitude) + 1)
            shot_id = np.array([frame + '_' + str(shot) for shot in shots])

            Data_.append(Data)
            Elevation_.append(Elevation)
            GPS_time_.append(GPS_time)
            Latitude_.append(Latitude)
            Longitude_.append(Longitude)
            Pitch_.append(Pitch)
            Roll_.append(Roll)
            Heading_.append(Heading)
            Bottom_.append(Bottom)
            Surface_.append(Surface)
            Shot_ID_.append(shot_id)


        Data               = pd.concat(Data_, axis=1, ignore_index=True)
        #Data               = Data[::-1]
        Time               = np.array(Data.index)
        GPS_time           = np.concatenate(GPS_time_)
        Elevation          = np.concatenate(Elevation_)
        Latitude           = np.concatenate(Latitude_)
        Longitude          = np.concatenate(Longitude_)
        Shot_ID            = np.concatenate(Shot_ID_)

        
        try:
            Bottom             = np.concatenate(Bottom_)
            Surface            = np.concatenate(Surface_)
        except:
            Bottom = 'empty'
            Surface = 'empty'


        # check if navigation data is available
        if navigation == 'off':
            Pitch   = 'empty'
            Roll    = 'empty'
            Heading = 'empty'

        else:
            Pitch              = np.concatenate(Pitch_)
            Roll               = np.concatenate(Roll_)
            Heading            = np.concatenate(Heading_)




        full_dict = {'Data'                 : Data.values,
                     'Time'                 : Time,
                     'Elevation'            : Elevation.T,
                     'GPS_time'             : GPS_time.T,
                     'Latitude'             : Latitude.T, 
                     'Longitude'            : Longitude.T,
                     'Pitch'                : Pitch,
                     'Roll'                 : Roll,
                     'Heading'              : Heading,
                     'Bottom'               : Bottom.T,
                     'Surface'              : Surface.T,
                     }
        scipy.io.savemat(output_filename, full_dict)
        print('===> Saved Frame as: {}'.format(output_filename))

    

    if z_mode == 'elevation':

        for file in frame_list:

            mat     = scipy.io.loadmat(file)


            if overlap == False:

                Data               = pd.DataFrame(np.array(mat['Data']))

                if Time_Array == False:
                    Time               = np.array(mat['Time'].T[0])
                elif Time_Array == True:
                    Time               = np.cumsum(np.repeat(1.666666666666657e-08, len(Data)))

                ol                 = Data.shape[1] - overlap_traces

                Elevation          = np.array(mat['Elevation'][0]) # aircraft elevation

                GPS_time           = np.array(mat['GPS_time'][0])
                Latitude           = np.array(mat['Latitude'][0])
                Longitude          = np.array(mat['Longitude'][0])
                Pitch              = np.array(mat['Pitch'][0])
                Roll               = np.array(mat['Roll'][0])
                Heading            = np.array(mat['Heading'][0])
                try:
                    Bottom             = np.array(mat['Bottom'][0])
                    Surface            = np.array(mat['Surface'][0])
                except:
                    Bottom = 'empty'
                    Surface = 'empty'


                Data.index         = Time


            frame   = file.split('.')[0].split('Data_')[1]
            shots   = np.arange(1, len(Longitude) + 1)
            shot_id = np.array([frame + '_' + str(shot) for shot in shots])

            Data_.append(Data)
            Elevation_.append(Elevation)
            GPS_time_.append(GPS_time)
            Latitude_.append(Latitude)
            Longitude_.append(Longitude)
            Pitch_.append(Pitch)
            Roll_.append(Roll)
            Heading_.append(Heading)
            Bottom_.append(Bottom)
            Surface_.append(Surface)
            Shot_ID_.append(shot_id)


        Data               = pd.concat(Data_, axis=1, ignore_index=True)
        #Data               = Data[::-1]
        Time               = np.array(Data.index)
        GPS_time           = np.concatenate(GPS_time_)
        Elevation          = np.concatenate(Elevation_)
        Latitude           = np.concatenate(Latitude_)
        Longitude          = np.concatenate(Longitude_)
        Shot_ID            = np.concatenate(Shot_ID_)


        try:
            Bottom             = np.concatenate(Bottom_)
            Surface            = np.concatenate(Surface_)
        except:
            Bottom = 'empty'
            Surface = 'empty'


        # check if navigation data is available
        if navigation == 'off':
            Pitch   = 'empty'
            Roll    = 'empty'
            Heading = 'empty'

        else:
            Pitch              = np.concatenate(Pitch_)
            Roll               = np.concatenate(Roll_)
            Heading            = np.concatenate(Heading_)




        full_dict = {'Data'                 : Data.values,
                     'Time'                 : Time,
                     'Elevation'            : Elevation.T,
                     'GPS_time'             : GPS_time.T,
                     'Latitude'             : Latitude.T, 
                     'Longitude'            : Longitude.T,
                     'Pitch'                : Pitch,
                     'Roll'                 : Roll,
                     'Heading'              : Heading,
                     'Bottom'               : Bottom.T,
                     'Surface'              : Surface.T,
                     'Shot_ID'              : Shot_ID.T
                     }

        scipy.io.savemat(output_filename, full_dict)
        print('===> Saved Frame as: {}'.format(output_filename))


##################################
# Flip single frame left to right
##################################


def flip_frame(file='', output_filename='', z_mode='elevation'):

    '''
    Reads a list of frames and connects them
    The list should have the following format:

        frames = ['Data_...._001_.mat', 'Data_...._002_.mat', 'Data_...._003_.mat', etc...]

    '''
        
    import scipy.io
    import pandas as pd
    import numpy as np
    import sys
    import glob
    import os


    # CHECKS

    if output_filename == '':
        print('you did not provide a output filename')
        print('setting output filename to output.mat')
        output_filename = 'output.mat'




    mat                = scipy.io.loadmat(file)


    navigation = 'off'

    Data               = pd.DataFrame(np.array(mat['Data'])).T[::-1].T
    Elevation          = np.array(mat['Elevation'])[0]
    GPS_time           = np.array(mat['GPS_time'])[0][::-1]
    Latitude           = np.array(mat['Latitude'])[0][::-1]
    Longitude          = np.array(mat['Longitude'])[0][::-1]
    X                  = np.array(mat['X'])[0][::-1]
    Y                  = np.array(mat['Y'])[0][::-1]
   # Aircraft_Elevation = np.array(mat['Aircraft_Elevation'])[0]
    Spacing            = np.array(mat['Spacing'])[0][::-1]
    Pitch              = np.array(mat['Pitch'])[0][::-1]
    Roll               = np.array(mat['Roll'])[0][::-1]
    Heading            = np.array(mat['Heading'])[0][::-1]
    try:
        Bottom             = np.array(mat['Bottom'])[0][::-1]
        Surface            = np.array(mat['Surface'])[0][::-1]
    except:
        Bottom = 'empty'
        Surface = 'empty'


    if navigation == 'off':

        Pitch   = 'empty'
        Roll    = 'empty'
        Heading = 'empty'


    Distance = np.cumsum(Spacing)


    full_dict = {'Data'                 : Data.values,
                 'Elevation'            : Elevation,
                 'GPS_time'             : GPS_time,
                 'Latitude'             : Latitude, 
                 'Longitude'            : Longitude,
                 'X'                    : X,
                 'Y'                    : Y,
                 #'Aircraft_Elevation'   : Aircraft_Elevation,
                 'Spacing'              : Spacing,
                 'Distance'             : Distance,
                 'Pitch'                : Pitch,
                 'Roll'                 : Roll,
                 'Heading'              : Heading,
                 'Bottom'               : Bottom,
                 'Surface'              : Surface
                 }

    scipy.io.savemat(output_filename, full_dict)
    print('===> Saved Frame as: {}'.format(output_filename))





##################################
# Truncate a frame
##################################


def truncate_frame(file='', output_filename='', start=0, end=-1, z_mode='elevation'):

    '''
    Reads a list of frames and connects them
    The list should have the following format:

        frames = ['Data_...._001_.mat', 'Data_...._002_.mat', 'Data_...._003_.mat', etc...]

    '''
        
    import scipy.io
    import pandas as pd
    import numpy as np
    import sys
    import glob
    import os


    # CHECKS

    if output_filename == '':
        print('you did not provide a output filename')
        print('setting output filename to output.mat')
        output_filename = 'output.mat'




    mat                = scipy.io.loadmat(file)


    navigation = 'off'

    Data               = pd.DataFrame(np.array(mat['Data'])).T[start:end].T
    Elevation          = np.array(mat['Elevation'])
    GPS_time           = np.array(mat['GPS_time'])[0][start:end]
    Latitude           = np.array(mat['Latitude'])[0][start:end]
    Longitude          = np.array(mat['Longitude'])[0][start:end]
    X                  = np.array(mat['X'])[0][start:end]
    Y                  = np.array(mat['Y'])[0][start:end]
   # Aircraft_Elevation = np.array(mat['Aircraft_Elevation'])[0]
    Spacing            = np.array(mat['Spacing'])[0][start:end]
    try:
        Pitch              = np.array(mat['Pitch'])[0][start:end]
        Roll               = np.array(mat['Roll'])[0][start:end]
        Heading            = np.array(mat['Heading'])[0][start:end]
    except:
        pass
    try:
        Bottom             = np.array(mat['Bottom'])[0][start:end]
        Surface            = np.array(mat['Surface'])[0][start:end]
    except:
        pass

    if navigation == 'off':

        Pitch   = 'empty'
        Roll    = 'empty'
        Heading = 'empty'


    Distance = np.cumsum(Spacing)


    full_dict = {'Data'                 : Data.values,
                 'Elevation'            : Elevation,
                 'GPS_time'             : GPS_time,
                 'Latitude'             : Latitude, 
                 'Longitude'            : Longitude,
                 'X'                    : X,
                 'Y'                    : Y,
                 #'Aircraft_Elevation'   : Aircraft_Elevation,
                 'Spacing'              : Spacing,
                 'Distance'             : Distance,
                 'Pitch'                : Pitch,
                 'Roll'                 : Roll,
                 'Heading'              : Heading,
                 'Bottom'               : Bottom,
                 'Surface'              : Surface
                 }

    scipy.io.savemat(output_filename, full_dict)
    print('===> Saved Frame as: {}'.format(output_filename))


################################
# Quickplot a radar section
# from a cresis radar .mat file
################################

def plot_mat(file, in_path='', out_path='', z_type='elevation', scale=15, flip=False, cmap='bone_r', legend=False, dpi=300):

    import h5py
    import scipy.io 
    import numpy as np
    import pandas as pd
    from matplotlib import pyplot as plt
    import os, glob
    import geopy.distance

    #from matplotlib.ticker import FormatStrFormatter

    if in_path == '':
        in_path = os.getcwd()
        print('No in_path defined, setting cwd as in_path')

    if out_path == '':
        out_path = os.getcwd()
        print('No out_path defined, setting cwd as out_path')


    # depending on the .mat file version
    # either scipy.io (older versions) or h5py (newer versions)
    # will be used to load the file        
    try:
        mat     = scipy.io.loadmat(file)
        reader  = 'scipy'

    except NotImplementedError:
        mat     = h5py.File(file, mode='r')
        reader  = 'h5py'
        

    # loading with scipy.io
    if reader == 'scipy':   
        print('Using scipy.io to load matfile')

        if z_type == 'twt':
            df        = pd.DataFrame(np.log10(np.array(mat['Data']))) # radar matrix
            twt       = np.array(mat['Time']).T
            longitude = np.array(mat['Longitude'])
            latitude  = np.array(mat['Latitude'])

        elif z_type == 'elevation':
            df        = pd.DataFrame(np.array(mat['Data'])) # radar matrix
            distance  = np.array(mat['Distance'][0]) # radar matrix
            elevation = np.array(mat['Elevation'][0])

        else:
            print("provide a z_type (z_type='twt' or 'elevation')")
            exit()
        
    #loading with h5py
    if reader == 'h5py':
        print('Using h5py to load matfile')
        
        if z_type == 'twt':
            df        = pd.DataFrame(np.log10(np.array(mat['Data']))).T # radar matrix
            twt       = np.array(mat['Time']).T
            longitude = np.array(mat['Longitude'])
            latitude  = np.array(mat['Latitude'])

        elif z_type == 'elevation':
            df        = pd.DataFrame(np.array(mat['Data'])) # radar matrix
            distance  = np.array(mat['Distance'][0]) # radar matrix
            elevation = np.array(mat['Elevation'][0])
        else:
            print("provide a z_type (z_type='twt' or 'elevation')")
            exit()


    if z_type == 'twt':

        ## Get real distance for traces 
        spacing = np.array([])
        for i in range(1, len(latitude)-1):
                coord_1    = (latitude[i], longitude[i])
                coord_2    = (latitude[i + 1], longitude[i + 1])
                f          = geopy.distance.geodesic(coord_1, coord_2).meters
                spacing    = np.append(spacing, f) 
                
        distance       = np.cumsum(spacing[0:-1])
        distance       = np.insert(distance, 0, 0)
        distance       = np.insert(distance, len(distance), distance[-1] + spacing.mean())
        spacing        = np.insert(spacing, 0, spacing.mean())
        spacing        = np.insert(spacing, len(spacing), spacing.mean())

    # 
    df = df.reset_index(drop=True) # reset index

    if flip == True:
        df = df.T[::-1].T
    
    elif flip == False:
        pass

        
    height      = scale
    width       = df.shape[1] / 150

    ## Plot Radargram
    fig, ax = plt.subplots(figsize=(width, height))
    
    ims = ax.imshow(df, cmap=cmap, aspect="auto")
    
    if z_type == 'elevation':

        x_step      = 250
        y_spacing   = 250
        offset      = int(np.min([y for y in elevation[0::y_spacing] if y > 0]))

        plt.xticks(np.array(range(1, len(distance), 500)), \
                   np.round((distance[0::500] / 1000), 0))
        plt.yticks(df.index.values[0::500], elevation[0::500])
        plt.xlabel('Along-track distance (km)', fontsize='16')
        plt.ylabel('Elevation (m)', fontsize='16')
        plt.title(file + ' ' + z_type, fontsize = '20')

    if z_type == 'twt':

        x_step      = 250
        #y_spacing   = 250
        #offset      = int(np.min([y for y in elevation[0::y_spacing] if y > 0]))

        plt.xticks(np.array(range(1, len(distance), 500)), \
                   np.round((distance[0::500] / 1000), 0))
        plt.yticks(df.index.values[0::100], np.round(twt[0::100] * 1000 * 1000, 0))
        plt.xlabel('Along-track distance (km)', fontsize='16')
        plt.ylabel('TWT (ns)', fontsize='16')
        plt.title(file + ' ' + z_type, fontsize = '20')
        
    if legend == True:
        plt.colorbar(ims)
    elif legend == False:
        pass
    else:
        pass 
            
    fig.savefig(file.split('.')[0] + '.png', dpi=dpi, bbox_inches='tight')












def twt2elevation_slow(data='',
                  twt='',
                  twt_surface='',
                  aircraft_elevation='',
                  speed_of_ice=1.689e8,
                  reference='GPS',
                  DEM_surface='',
                  setting='narrowband',
                  overlap=False,
                  overlap_traces=0,
                  decimate=[False, 1]
                  ):


    import pandas as pd
    import numpy  as np


    speed_of_light = 2.99792458e8
    speed_of_ice   = speed_of_ice

    data = data

    # define data frames
    df       = data
    df       = df.apply(pd.to_numeric).astype(float)        
    df       = df.reset_index(drop=True) # reset index
    df_comb  = pd.DataFrame(columns = ['Elevation', 'dB', 'Trace'])

    # define some short variables
    twt      = twt                 # Time array
    elev     = aircraft_elevation  # Aircraft Elevation
    twt_surf = twt_surface         # twt of surf. reflection


    
    # Log every 500 lines.
    LOG_EVERY_N = 1000

    # create surface and bottom array (m)
    surface_m = []

    for i in np.arange(0, elev.size):
        if (i % LOG_EVERY_N) == 0:
            print('... processed  {}  of  {}  traces'.format(i + 1, elev.size))

        # get index where surface reflection is located
        if setting == 'emr':
            surf_idx = idx[i]
        else:
            surf_idx = (np.abs(np.array(twt) - np.array(twt_surf)[i])).argmin()

        # get single trace of radargram
        data = np.array(df[i])

        # delete values until surface reflection
        data = np.delete(data,np.s_[0:surf_idx])

        #
        T                   = np.delete(twt,np.s_[0:surf_idx],axis=0) # delete traces abofe surface reflection
        T_                  = T - twt[surf_idx] # set start of new twt array to zero
        Depth               = T_ * (speed_of_ice / 2)
        Air_Column          = (twt_surf[i] * speed_of_light) / 2
        Airplane_Elevation  = elev[i]

        if reference == 'GPS':
            surface   = Airplane_Elevation - Air_Column
            Elevation = Airplane_Elevation - Air_Column - Depth
        elif reference == 'DEM':
            surface   = DEM_surface.flatten()
            Elevation = surface[i] - Depth

        surface_m.append(surface)

        trace               = np.ones(data.size) * i
        df_trace            = pd.DataFrame(Elevation)
        df_trace['dB']      = pd.DataFrame(data)
        df_trace['Trace']   = pd.DataFrame(trace)
        df_trace.columns    = ['Elevation', 'dB', 'Trace']
        df_comb             = df_comb.append(df_trace)

    # drop nan's
    df_comb = df_comb.dropna()

    # create final data frame
    # the range interval depends on the frequency range
    if setting == 'wideband': 
        df_comb         = df_comb.round({'Elevation': 1})

        ## create pivot table for heatmap
        df = df_comb.pivot('Elevation', 'Trace', 'dB')
        df = df.interpolate()
        df = df.iloc[::-1]

        # special section for SEGY conversion:
        # traces with a length of more than 32767 samples are not supported
        if len(df) > 32767 :
            print('...(problem) ===> Number of samples: {}'.format(len(df)))
            print('...(problem) ===> Exceeding maximum number of samples for SEGY conversion ( >32767 samples )')
            print('...(action)  ===> Deleting every >>{}<< row'.format(int(decimate[1])))

            # decimate if it is set TRUE
            if decimate[0] == False:
                pass
            elif decimate[0] == True:
                decimate_factor = int(decimate[1])
                df = df.iloc[::decimate_factor, :]

    if setting == 'narrowband' or setting == 'emr':
        df_comb = df_comb.dropna()
        df_comb = df_comb.round({'Elevation': 0})

        ## create pivot table for heatmap
        df = df_comb.pivot('Elevation', 'Trace', 'dB')
        df = df.interpolate()
        df.index = df.index.astype(int)
        df = df.iloc[::-1]

    if setting == 'snow':
        df_comb         = df_comb.round({'Elevation': 3})

        ## create pivot table for heatmap
        df = df_comb.pivot('Elevation', 'Trace', 'dB')
        df = df.interpolate()
        df.index = np.around(df.index.values, decimals=2)
        df = df.loc[~df.index.duplicated(keep='first')]
        df = df.iloc[::-1]

    # get Elevation (Z) array
    Z = df.index.values

    if reference == 'DEM':
        surface_m = DEM_surface

    # get the index value of surface reflection
    surface_m = np.array(surface_m).flatten()
    surf_m_idx = np.array([])

    for i in range(len(surface_m)):
        s_idx = (np.abs(df.index.values - np.array(surface_m)[i])).argmin()
        surf_m_idx = np.append(surf_m_idx, s_idx)

    surf_m_idx = surf_m_idx.astype(int)

    # delete crappy traces above surface reflection
    data     = np.array(df.T)
    new_data = []

    for i in range(len(surf_m_idx)):
        trace         = data[i]
        index         = surf_m_idx[i]
        trace[:index] = np.nan
        new_data.append(trace)

    new_data = np.array(new_data).T
    df = pd.DataFrame(new_data)

    if overlap == True:
        df.drop(df.columns[-65:], axis=1, inplace=True)
        df_meta.drop(df_meta.index[-65:], axis=0, inplace=True)

    print('==> Done ...')
    
    return df, Z, surface_m, surf_m_idx