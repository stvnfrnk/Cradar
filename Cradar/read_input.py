
'''


'''


###################
# CReSIS Matfiles
###################

def read_cresis_mat(matfile):

    import h5py
    import scipy.io
    import numpy as np
    import pandas as pd


    ##############################
    # Try with scipy.io.loadmat()
    try:
        File      = scipy.io.loadmat(matfile)
        Frame     = matfile.split('.mat')[0]
        Reader    = 'scipy.io'
        Domain    = 'twt'

        Data      = pd.DataFrame(np.array(File['Data']))
        Time      = np.array(File['Time']).flatten()
        Longitude = np.array(File['Longitude']).flatten()
        Latitude  = np.array(File['Latitude']).flatten()
        Elevation = np.array(File['Elevation']).flatten()
        GPS_time  = np.array(File['GPS_time']).flatten()
        
        Layer = {}
        if 'Surface' in File.keys():
            surface    = {'trace'  : np.arange(len(np.array(File['Surface']).flatten())) + 1,
                            'value'  : np.array(File['Surface']).flatten(),
                            'color'  : 'blue'}
            Layer['Surface'] = surface

        if 'Bottom' in File.keys():
            bottom    = {'trace'  : np.arange(len(np.array(File['Bottom']).flatten())) + 1,
                         'value'  : np.array(File['Bottom']).flatten(),
                         'color'  : 'red'}
            Layer['Bottom'] = bottom

    ##############################
    # Try with h5py.File()
    except:
        File      = h5py.File(matfile, 'r')
        Frame     = matfile.split('.mat')[0]
        Reader    = 'h5py'
        Domain    = 'twt'

        Data      = pd.DataFrame(np.array(File['Data'])).T
        Time      = np.array(File['Time']).flatten()
        Longitude = np.array(File['Longitude']).flatten()
        Latitude  = np.array(File['Latitude']).flatten()
        Elevation = np.array(File['Elevation']).flatten()
        GPS_time  = np.array(File['GPS_time']).flatten()
        
        Layer = {}
        if 'Surface' in File.keys():
            surface    = {'trace'  : np.arange(len(np.array(File['Surface']).flatten())) + 1,
                            'value'  : np.array(File['Surface']).flatten(),
                            'color'  : 'blue'}
            Layer['Surface'] = surface

        if 'Bottom' in File.keys():
            bottom    = {'trace'  : np.arange(len(np.array(File['Bottom']).flatten())) + 1,
                         'value'  : np.array(File['Bottom']).flatten(),
                         'color'  : 'red'}
            Layer['Bottom'] = bottom
        

    print('Loaded {} with {}'.format(Frame, Reader))

    return Frame, Reader, Domain, Data, Time, Longitude, Latitude, Elevation, GPS_time, Layer

################################################################################################
################################################################################################

###################
# AWI SEGY Files
###################

def read_awi_segy(segy_file):

    # load_awi_segy(self, segy_file='', coordinate_file='', Longitude='', Latitude='', dB=False, correct_gps=True):

    '''

    '''

    from obspy.io.segy.segy import _read_segy
    import numpy as np
    import pandas as pd

    stream  = _read_segy(segy_file, headonly=True)
    data    = pd.DataFrame(np.array([t.data for t in list(stream.traces)]))

    header            = stream.binary_file_header.__dict__
    frame             = header['line_number']
    sample_interval   = header['sample_interval_in_microseconds']
    sample_interval   = sample_interval * 10e-13
    number_of_samples = header['number_of_samples_per_data_trace']

    # build time array
    time    = np.repeat(sample_interval, number_of_samples)
    time[0] = 0
    time    = np.cumsum(time)

    Data   = data.T
    Stream = stream
    Time   = time
    Frame  = str(frame)

    return Data, Time, Frame, Stream




def read_awi_nc(nc_file):

    # load_awi_segy(self, segy_file='', coordinate_file='', Longitude='', Latitude='', dB=False, correct_gps=True):

    '''

    '''

    import xarray as xr
    import numpy as np

    dx =  xr.load_dataset(nc_file)

    Data       = dx.variables['WAVEFORM'].values.T
    Longitude  = dx.variables['LONGITUDE'].values
    Latitude   = dx.variables['LATITUDE'].values

    Aircraft_altitude      = dx.variables['ALTITUDE'].values
    Ice_surface_elevation  = dx.variables['ELEVATION'].values
    Range                  = dx.variables['RANGE'].values

    sample_interval   = 1.3333 * 10e-9
    number_of_samples = dx.attrs['SAMPLES']

    # build time array
    time    = np.repeat(sample_interval, number_of_samples)
    time[0] = 0
    Time    = np.cumsum(time)

    surface_twt = ( (Aircraft_altitude - Ice_surface_elevation) / 299792458 ) * 2

    Layer = {}
    surface  = {'trace'  : np.arange(len(np.array(surface_twt).flatten())) + 1,
                'value'  : np.array(surface_twt).flatten(),
                'color'  : 'blue'}
    Layer['Surface'] = surface

    return Data, Time, Longitude, Latitude, Aircraft_altitude, Ice_surface_elevation, Layer
