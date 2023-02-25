


class Cradar:

    #_allObjects = []

    # initiate an instance and read all important parameters
    # from a CReSIS mcords matfile
    def __init__(self):
        #self._allObjects.append(self)
        pass

    
    # def create_traces(self):
    #     import numpy as np

    #     traces = np.array(range(len(self.Longitude))).astype(int)
    #     return traces

    #############################
    # Method: load_cresis_mat
    #############################

    def load_cresis_mat(self, filename, dB=False):

        import h5py
        import scipy.io
        import numpy as np
        import pandas as pd


        ##############################
        # Try with scipy.io.loadmat()
        try:
            self.File      = scipy.io.loadmat(filename)
            self.Frame     = filename.split('.mat')[0]
            try:
                self.Frame = self.Frame.split('/')[-1]
            except:
                pass

            self.Reader    = 'scipy'

            # Iterate over almost all items in HDF5 File
            for k, v in self.File.items():
                if 'Time' not in k:
                    try:
                        setattr(self, k, np.array(v).flatten())
                    except:
                        setattr(self, k, v)

                if 'Time' in k:
                    self.Time = np.array(self.File['Time']).flatten()

                if 'Time' in k and 'Z' not in k:
                    self.Domain = 'twt'

                if 'Z' in k:
                    self.Domain = 'Z'

            self.Data = pd.DataFrame(np.array(self.File['Data']))

            if dB == False:
                self.dB = False
            elif dB == True:
                self.dB = True
            else:
                print('==> dB True or False? Is set to False.')




        ######################
        # Try with h5py.File()
        except:
            self.File      = h5py.File(filename, 'r')
            self.Frame     = filename.split('.mat')[0]
            self.Reader    = 'h5py'

            # Iterate over almost all items in HEF5 File
            for k, v in self.File.items():
                if 'Time' not in k:
                    if '#' not in k:
                        setattr(self, k, np.array(v).flatten())

                if 'Time' in k:
                    self.Time = np.array(self.File['Time']).flatten()

                if 'Time' in k and 'Z' not in k:
                    self.Domain    = 'twt'

            self.Data = pd.DataFrame(np.array(self.File['Data'])).T

            if dB == False:
                self.dB = False
            elif dB == True:
                self.dB = True
            else:
                print('==> dB True or False? Is set to False.')

        # Create Traces
        self.Trace = np.array(range(len(self.Longitude))).astype(int) + 1

        self.Layer = {}

        # Delete the HDF5 file
        del self.File
        # print('')
        print('==> Loaded {}'.format(self.Frame))
        return self

    ########## END of load_cresis_mat() ###########





    #############################
    # Method: load_h5
    #############################

    def load_h5(self, filename, dB=False):

        import h5py
        import numpy as np
        import pandas as pd

        self.File      = h5py.File(filename, 'r')

        # Iterate over almost all items in HEF5 File
        for k, v in self.File.items():
            if 'Time' not in k:
                setattr(self, k, np.array(v).flatten())

            if 'Time' in k:
                self.Time = np.array(self.File['Time']).flatten()

            if 'Time' in k and 'Z' not in k:
                self.Domain    = 'twt'

        self.Data = pd.DataFrame(np.array(self.File['Data']))
        self.Frame     = filename.split('.h5')[0]
        self.Reader    = 'h5py'

        if dB == False:
            self.dB = False
        elif dB == True:
            self.dB = True
        else:
            print('==> dB True or False? Is set to False.')

        # Delete the HDF5 file
        del self.File
        # print('')
        print('==> Loaded {}'.format(self.Frame))

        return self

    ########## END of load_h5() ###########




    #############################
    # Method: load_awi_segy
    #############################

    def load_awi_segy(self, segy_file='', coordinate_file='', dB=False, correct_gps=True):

        '''


        '''


        from obspy.io.segy.segy import _read_segy
        import numpy as np
        import pandas as pd

        segy_file       = segy_file
        coordinate_file = coordinate_file

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


        self.Data   = data.T
        self.Stream = stream
        self.Time   = time
        self.Frame  = str(frame)
        self.Domain = 'twt'

        coords         = pd.read_csv(coordinate_file, delim_whitespace=True)
        self.Longitude = coords['GPSLon'].values
        self.Latitude  = coords['GPSLat'].values
        self.Elevation = coords['GPSAlt'].values
        self.GPS_time  = coords['Time'].values

        # if correct_gps == True:
        #     try:
        # 		self.Latitude  = np.array( pd.DataFrame(self.Latitude).mask(pd.DataFrame(self.Latitude).duplicated(keep='first'), np.nan).interpolate() ).T[0]
        #         self.Longitude = np.array( pd.DataFrame(self.Longitude).mask(pd.DataFrame(self.Longitude).duplicated(keep='first'), np.nan).interpolate() ).T[0]
        # 	except:
        # 		print('... could not correct gps positions.')

        # if dB == False:
        #         self.dB = False
        # elif dB == True:
        #     self.dB = True
        # else:
        #     print('==> dB True or False? Is set to False.')

        # print('')
        print('==> Loaded {}'.format(self.Frame))

        del stream, data, header, time, coords

        return self

    ########## END of load_awi_segy() ###########





    def load_awi_nc(self, nc_file='', dB=True):

        '''


        '''


        import numpy as np
        import pandas as pd
        import xarray as xr

        
        dx =  xr.load_dataset(nc_file)

        self.Data   = dx.variables['WAVEFORM'].values[::-1]
        self.Time   = dx.variables['TWT'].values
        # self.Frame  = str(frame)
        self.Domain = 'twt'

        self.Longitude = dx.variables['LONGITUDE'].values
        self.Latitude  = dx.variables['LATITUDE'].values
        
        print('==> Loaded {}'.format(nc_file))

        return self

    ########## END of load_awi_nc() ###########




    #############################
    # Method: load_bas_nc
    #############################

    def load_bas_nc(self, filename, data_type='pulse'):

        import xarray as xr
        import numpy as np
        import pandas as pd
        import os

        if data_type == 'pulse':

            ds = xr.open_dataset(filename, decode_times=False)

            self.Frame       = os.path.split(filename)[1].split('.nc')[0]
            self.Time        = ds['fast_time'].values.astype(float) * 10e-7
            self.Data        = pd.DataFrame(ds.variables['pulse_data'].values)
            self.Domain      = 'twt'
            self.dB          = False

            self.Longitude   = ds.variables['longitude_layerData'].values
            self.Latitude    = ds.variables['latitude_layerData'].values
            self.Elevation   = ds.variables['aircraft_altitude_layerData'].values
            self.GPS_time    = ds.variables['UTC_time_layerData'].values

            self.Surface_idx = ds.variables['surface_pick_layerData']
            self.Surface_m   = ds.variables['surface_altitude_layerData'].values
            self.Bed_idx     = ds.variables['bed_pick_layerData']
            self.Bed_m       = ds.variables['bed_altitude_layerData'].values

            print('')
            print('==> Loaded {} >> pulse data <<'.format(self.Frame))

            del ds

            return self

        if data_type == 'chirp':

            from bas_io import correct_chirp_data

            data_chirp, time, longitude, latitude, elevation, gps_time, surface_idx, surface_m, bed_idx, bed_m = correct_chirp_data(filename)

            self.Frame       = os.path.split(filename)[1].split('.nc')[0]
            self.Time        = time
            self.Data        = data_chirp
            self.Domain      = 'twt'
            self.dB          = False

            self.Longitude   = longitude
            self.Latitude    = latitude
            self.Elevation   = elevation
            self.GPS_time    = gps_time

            self.Surface_idx = surface_idx
            self.Surface_m   = surface_m
            self.Bed_idx     = bed_idx
            self.Bed_m       = bed_m

            print('')
            print('==> Loaded {} >> chirp data <<'.format(self.Frame))

            del data_chirp, time, longitude, latitude, elevation, gps_time, surface_idx, surface_m, bed_idx, bed_m

            return self

    ########## END of load_bas_nc() ###########
        


    #############################
    # Method: load_awi_segy
    #############################

    def load_custom(self,
                    data='',
                    time='',
                    frame='',
                    domain='',
                    longitude='',
                    latitude='',
                    dB=False):

        '''


        '''

        import numpy as np


        data[np.isnan(data)] = 0
        print(data.min().min())
        # try:
        #     data[np.isnan(data)] = 0
        # except:
        #     print('error')

        print(data.min().min())

        self.Data      = data
        self.Time      = time
        self.Frame     = frame
        self.Domain    = domain
        self.Longitude = longitude
        self.Latitude  = latitude
        self.dB        = dB

        return self


    #############################
    # Method: LAyERS
    #############################

    def add_layer_by_frame_trace(self, layer_name, traces, values):

        layer = {'trace'      : traces,
                 'value_twt'  : values}
        
        try:
            self.Layer
        except:
            self.Layer             = {}
        
        self.Layer[layer_name] = layer
        print('==> added layer: {}'.format(layer_name))


    def add_layer_by_coords(self, layer_name, layer_lon, layer_lat, layer_val):

        import numpy as np
        from scipy import spatial

        crd_locs    = list(zip(self.Longitude, self.Latitude))
        tree        = spatial.KDTree(crd_locs)
        npt_traces  = []
        npt_values  = []

        for i in range(len(layer_lon)):
            pick_loc  = (layer_lon[i], layer_lat[i])
            val       = layer_val[i]
            npt_trace = tree.query([pick_loc])[1][0]

            npt_traces.append(npt_trace)
            npt_values.append(val)

        traces = np.array(npt_traces)
        values = np.array(npt_values)

        # delete duplicates
        idx    = np.unique(traces, return_index=True)[1]
        traces = traces[idx]
        values = values[idx]

        layer = {'trace'      : traces,
                 'value_twt'  : values}

        self.Layer[layer_name] = layer
        print('==> added layer: {}'.format(layer_name))


    def get_layer_idx(self):

        import numpy as np

        layer_list = list(self.Layer.keys())
        time        = self.Time

        for lr in layer_list:
            values_idx = np.array([])
            v_twt      = self.Layer[lr]['value_twt']

            for i in range(len(v_twt)):
                v_idx      = (np.abs(np.array(time) - np.array(v_twt)[i])).argmin()
                values_idx = np.append(values_idx, v_idx)
                #print(surf_idx)
            
            values_idx[values_idx == 0] = np.nan

            self.Layer[lr]['value_twt_idx'] = values_idx
            print('... getting layer idx for {}'.format(lr))

    
    
            
            
    def track_surface(self, skip=100, gauss_factor=1, use_gradient=False, offset=2):

        '''

        '''

        import numpy as np
        import pandas as pd
        from scipy.ndimage import gaussian_filter

        data   = self.Data
        offset = offset
        #sample_interval = float(int(sample_interval) / 1000)

        ### Find Surface Reflection id's
        factor = 10

        # filter to avoid large jumps
        data = pd.DataFrame(gaussian_filter(data, sigma=gauss_factor))
        if use_gradient == False:
            data_filtered  = data.rolling(factor, center=True, win_type='hamming').mean()
        else:
            data_filtered  = data.diff().rolling(factor, center=True, win_type='hamming').mean()
        surf_idx       = ( data_filtered[skip::].idxmax(axis=0, skipna=True) ).astype(int)

        # get surface values
        srf = []
        
        for i in range(len(surf_idx)):
            index = surf_idx[i]
            val   = self.Time[index]
            srf.append(val)
        self.Surface = np.array(srf)


        # construct time (twt) array

        #twt_ns = np.ones(data.shape[0]) * sample_interval # in nanoseconds
        #twt_   = twt_ns / 10e8 # in seconds
        #twt_s  = np.cumsum(twt_)

        #self.Time2        = twt_s
        self.Surface_idx = np.array(surf_idx)



        del data, data_filtered, surf_idx#, twt_ns, twt_

    def emr_preprocess(self, skip=100, gauss_factor=1):

        '''

        '''

        import numpy as np
        import pandas as pd
        from scipy.ndimage import gaussian_filter

        stream = self.Stream
        data   = self.Data

        sample_interval = str(stream.binary_file_header).split('sample_interval_in_microseconds: ')[1].split('sample_interval_in_microseconds_of_original_field_recording')[0].split('\n')[0]
        sample_interval = float(int(sample_interval) / 1000)

        ### Find Surface Reflection id's
        factor = 10

        # filter to avoid large jumps


        data = pd.DataFrame(gaussian_filter(data.values, sigma=gauss_factor))
        data_filtered  = data.diff().rolling(factor, center=True, win_type='hamming').mean()
        surf_idx       = ( data_filtered[skip::].idxmax(axis=0, skipna=True) - (factor/2) ).astype(int)

        # get surface values
        srf = []
        for i in range(len(surf_idx)):
            index = surf_idx[i]
            val   = self.Time[index]
            srf.append(val)
        self.Surface = np.array(srf)


        # construct time (twt) array

        #twt_ns = np.ones(data.shape[0]) * sample_interval # in nanoseconds
        #twt_   = twt_ns / 10e8 # in seconds
        #twt_s  = np.cumsum(twt_)

        #self.Time2        = twt_s
        self.Surface_idx = np.array(surf_idx)



        del stream, data, data_filtered, surf_idx#, twt_ns, twt_


    #############################
    # Method: add_raster_values
    #############################

    def grdtrack(self, geotif='', geotif_name='DEM_surface', geotif_epsg=3413):


        geotif      = geotif
        geotif_name = geotif_name
        Longitude   = self.Longitude
        Latitude    = self.Latitude

        #print(geotif_epsg)


        print('==> Applying gridtrack method 1 ...')
        from geo_toolbox import gridtrack
        raster_vals = gridtrack(Longitude=Longitude, Latitude=Latitude, geotif=geotif, geotif_name=geotif_name, geotif_epsg=geotif_epsg)

        #print('==> Failed, trying gridtrack method 2 ...')
        #from geo_toolbox import gridtrack2
        #raster_vals = gridtrack(geotif, Latitude, Longitude, EPSG_xy=4326, EPSG_raster=geotif_epsg)

        setattr(self, geotif_name, raster_vals)
        print('==> Added {} to the data'.format(geotif_name))

        del raster_vals, Longitude, Latitude




    ########## END of rename() ###########










    #############################
    # Method: get_surf_idx()
    #############################

    def get_surf_idx(self):

        import numpy as np

        twt      = self.Time
        twt_surf = self.Surface
        #traces   = self.Data.shape[1]

        surf_idx = np.array([])

        for i in range(len(twt_surf)):
            s_idx = (np.abs(np.array(twt) - np.array(twt_surf)[i])).argmin()
            surf_idx = np.append(surf_idx, s_idx)
            #print(surf_idx)

        self.Surface_idx = surf_idx.astype(int)
        print('==> Added pixel index of surface reflection')
        del twt, twt_surf, surf_idx


    ########## END of get_surf_idx() ###########




    #############################
    # Method: get_surf_m_idx()
    #############################

    def get_surf_m_idx(self):

        import numpy as np

        Z          = self.Z
        Z_surf     = self.Surface_m
        surf_m_idx = np.array([])

        for i in range(len(Z_surf)):
            idx        = len(Z) - ( (np.abs(np.array(Z) - np.array(Z_surf)[i])).argmin() )
            surf_m_idx = np.append(surf_m_idx, idx)

        self.Surface_m_idx = surf_m_idx
        print('==> Added pixel index of surface elevation')
        del Z, Z_surf, surf_m_idx


    ########## END of get_surf_idx() ###########


    #############################
    # Method: get_bed_m_idx()
    #############################

    def get_bed_m_idx(self):

        import numpy as np

        Z         = self.Z
        Z_bed     = self.Bed_m
        bed_m_idx = np.array([])

        for i in range(len(Z_bed)):
            idx        = ( (np.abs(np.array(Z) - np.array(Z_bed)[i])).argmin() )
            bed_m_idx  = np.append(bed_m_idx, idx)

        self.Bed_m_idx = bed_m_idx
        print('==> Added pixel index of bed elevation')
        del Z, Z_bed, bed_m_idx


    ########## END of get_bed_idx() ###########



    #############################
    # Method: retrack_surf()
    #############################

    def retrack_surf(self, roll_factor=10, sigma=2, gatesize=50, offset=2):

        import copy
        import numpy as np
        import pandas as pd
        from scipy.ndimage import gaussian_filter

        roll_factor = roll_factor
        sigma       = sigma
        gatesize    = gatesize
        offset      = offset

        self.Surface_idx_old = copy.deepcopy(self.Surface_idx)
        self.Surface_old     = copy.deepcopy(self.Surface)

        data_filtered = pd.DataFrame(gaussian_filter(self.Data.values, sigma=2))
        data_filtered = data_filtered.diff().rolling(roll_factor, center=True, win_type='hamming').mean()

        # get surface values
        srf      = []
        srf_idx  = []

        for i in range(len(self.Longitude)):
            trx         = np.array(data_filtered[i])[self.Surface_idx[i] - gatesize:self.Surface_idx[i] + gatesize]
            offset_gate = np.argmax(trx) - gatesize
            index       = self.Surface_idx[i] + offset_gate + offset
            val         = self.Time[index]
            srf.append(val)
            srf_idx.append(index)

        self.Surface     = np.array(srf)
        self.Surface_idx = np.array(srf_idx)

        print('==> Retracked ice surface')

        #return self.Surface_idx_old








    #############################
    # Method: correct_geom_spreading
    #############################

    def correct4attenuation(raw_object, mode=0, loss_factor=0):

        from radar_toolbox import correct4attenuation
        import copy

        if raw_object.dB == False:
            pass
        else:
            print('==> Data should not be in dB on input.')

        mode        = mode
        loss_factor = loss_factor

        try:
            geom_obj = copy.deepcopy(raw_object)
        except:
            geom_obj = copy.copy(raw_object)

        geom_obj.get_surf_idx()

        data     = geom_obj.Data
        twt      = geom_obj.Time
        surf_idx = geom_obj.Surface_idx

        data_new = correct4attenuation(data, twt, surf_idx, v_ice=1.68914e8, mode=mode, loss_factor=loss_factor)

        geom_obj.Data = data_new
        geom_obj.dB   = True

        return geom_obj

        del geom_obj, data, twt, surf_idx, data_new


    ########## END of get_surf_idx() ###########





    #############################
    # Method: add_distance
    #############################


    def add_distance(self):

        from geo_toolbox import coords2distance

        X = self.Longitude
        Y = self.Latitude

        Spacing, Distance = coords2distance(X, Y, EPSG=4326)

        self.Spacing  = Spacing
        self.Distance = Distance

        del X, Y, Spacing, Distance
        print('==> Added Spacing and Distance.')


    ########## END of add_distance() ###########




    #############################
    # Method: calc_elevation
    #############################

    def twt2elevation(twt_object,
                       reference='',
                       setting='',
                       speed_of_ice=1.689e8,
                       overlap=False,
                       number_of_gaps=100,
                       decimate=[True, 2]):

        '''
        ==> Takes a Cradar object and transforms from twt domain
            to the elvation domain.

        ==> The air-ice interface can either be set by the twt from the
            aircraft to the surface reflection or by a ice surface DEM

        reference       = 'GPS' or 'DEM'
        speed_of_ice    = Usually 1.689 for e=3.15, but it can be changed
        setting         = 'narrowband', 'wideband' for rds or 'snow' for uwbm-snowradar
        overlap         = some older CReSIS .mat files used to be processed with an overlap
                          Hence, if you for instance concatenate several frames, you will
                          end up with multiple segments...
        number_of_gaps  = refers to the (by the user) defined number of allowed gaps in a DEM
                          for the twt to elevation conversion. If the number is exceeded, it will
                          switch automatically to reference='reflection'
        '''

        import numpy as np
        import pandas as pd
        from radar_toolbox import twt2elevation
        import copy

        # makes a copy of the first object (serves as a blue print)
        try:
            elev_obj = copy.deepcopy(twt_object)
        except:
            elev_obj = copy.copy(twt_object)

        print('==> Now: twt2elevation...')

        if reference == 'GPS':
            print('... Using aircraft GPS and radar surface reflection to derive elevation')

        elif reference == 'DEM':
            print('... Using a ice surface DEM to derive elevation')

        elif reference == '':
            reference = 'GPS'
            print("... !! You didn't define a reference, it is now automatically set to 'GPS'")
            print('... Using aircraft GPS and radar surface reflection to derive elevation')

        # Get ice surface elevation values from DEM
        if reference == 'DEM':
            twt_surface = elev_obj.Surface
            DEM_surface = elev_obj.DEM_surface
            aircraft_elevation = np.ones(len(elev_obj.Longitude))

        elif reference == 'GPS':
            twt_surface = elev_obj.Surface
            DEM_surface = ''
            aircraft_elevation = elev_obj.Elevation

        # imput variables from instance
        data               = elev_obj.Data
        twt                = elev_obj.Time
        

        # input variables defined in prior steps
        twt_surface   = twt_surface

        # input variables from input of method above
        reference = reference
        setting   = setting
        overlap   = overlap

        df, Z, surf_m, surf_m_idx =  twt2elevation(data=data,
                                                   twt=twt,
                                                   twt_surface=twt_surface,
                                                   aircraft_elevation=aircraft_elevation,
                                                   speed_of_ice=1.689e8,
                                                   reference=reference,
                                                   DEM_surface=DEM_surface,
                                                   setting=setting,
                                                   overlap=overlap,
                                                   overlap_traces=0,
                                                   decimate=[True, 2]
                                                   )


        # re-define instance atributes
        elev_obj.Z             = Z
        elev_obj.Data          = pd.DataFrame(np.array(df))
        elev_obj.Surface_m     = surf_m
        elev_obj.Surface_m_idx = surf_m_idx
        elev_obj.Domain        = 'Z'

        # define range resolution depending on the setting
        if setting == 'narrowband':
            elev_obj.Range_Resolution = '1 m'

        if setting == 'wideband':
            elev_obj.Range_Resolution = '0.1 m'

        if setting == 'snow':
            elev_obj.Range_Resolution = '0.001 m'


        return elev_obj

        del df


    ########## END of twt2elevation() ###########







#############################
    # Method: pull2surface
#############################

    def pull2surface(twt_object, setting=''):

        '''
        
        '''

        import numpy as np
        import pandas as pd
        from radar_toolbox import radar_pull2surface
        import copy

        # makes a copy of the first object (serves as a blue print)
        try:
            p2s_obj = copy.deepcopy(twt_object)
        except:
            p2s_obj = copy.copy(twt_object)

        print('==> Now: pull2elevation...')

        # imput variables from instance
        data        = p2s_obj.Data
        twt         = p2s_obj.Time
        twt_surface = p2s_obj.Surface

        # input variables defined in prior steps
        twt_surface   = twt_surface

        df, new_time_array =  radar_pull2surface(data=data,
                                                 twt=twt,
                                                 twt_surface=twt_surface,
                                                 setting=setting
                                                 )


        # re-define instance atributes
        p2s_obj.Time   = new_time_array
        p2s_obj.Data   = pd.DataFrame(np.array(df))
        p2s_obj.Domain = 'twt'

        del df

        return p2s_obj

        


    ########## END of pull2surface() ###########




    #############################
    # Method: pull2bed
    #############################

    def pull2bed(Z_object):

        '''
        
        '''

        import numpy as np
        import pandas as pd
        from radar_toolbox import radar_pull2bed
        import copy

        # makes a copy of the first object (serves as a blue print)
        try:
            p2b_obj = copy.deepcopy(Z_object)
        except:
            p2b_obj = copy.copy(Z_object)

        print('==> Now: pull2bed...')

        # imput variables from instance
        data             = p2b_obj.Data
        elevation_array  = p2b_obj.Z
        bed_elevation    = p2b_obj.Bed_m

        # define range resolution depending on the setting
        if p2b_obj.Range_Resolution == '1 m':
            range_resolution_m = 1

        if p2b_obj.Range_Resolution == '0.1 m':
            range_resolution_m = 0.1

        if p2b_obj.Range_Resolution == '0.0001 m':
            range_resolution_m = 0.001


        df, depth_array = radar_pull2bed(data=data, 
                                         elevation_array=elevation_array, 
                                         bed_elevation=bed_elevation, 
                                         range_resolution_m=range_resolution_m
                                         )

        # re-define instance atributes
        p2b_obj.Depth  = depth_array
        p2b_obj.Data   = pd.DataFrame(np.array(df))
        p2b_obj.Domain = 'depth'

        return p2b_obj

        del df


    ########## END of pull2surface() ###########



    #############################
    # Method: flip left-right
    #############################

    def flip_lr(self):

        '''
        Flips the whole radar matrix and all its along-track attributes
        '''

        self.Data      = self.Data.T[::-1].T
        self.Longitude = self.Longitude[::-1]
        self.Latitude  = self.Latitude[::-1]
        self.GPS_time  = self.GPS_time[::-1]
        self.Elevation = self.Elevation[::-1]

        if self.Domain == 'twt':
            self.Surface     = self.Surface[::-1]
            try:
                self.Surface_idx = self.Surface_idx[::-1]
            except:
                pass

        elif self.Domain == 'Z':
            self.Surface_m     = self.Surface_m[::-1]
            try:
                self.Surface_m_idx = self.Surface_m_idx[::-1]
            except:
                pass


        # optional
        try:
            self.Heading = self.Heading[::-1]
        except:
            pass
        try:
            self.Pitch   = self.Pitch[::-1]
        except:
            pass
        try:
            self.Roll    = self.Roll[::-1]
        except:
            pass


    ########## END of flip_lr() ###########




    ##################################
    # Method: clip data along-track
    ##################################

    def clip_along(self, start=0, end=-1):

        '''
        '''

        import pandas as pd
        import numpy as np

        self.Data      = self.Data.T[start:end].T

        # tweak to make col. names start with zero
        # important later for indexing
        self.Data      = pd.DataFrame(np.array(self.Data))

        self.Longitude = self.Longitude[start:end]
        self.Latitude  = self.Latitude[start:end]
        

        # optional
        try:
            self.Elevation = self.Elevation[start:end]
        except:
            pass
        try:
            self.GPS_time  = self.GPS_time[start:end]
        except:
            pass
        try:
            self.Surface   = self.Surface[start:end]
        except:
            pass
        try:
            self.Surface_idx = self.Surface_idx[start:end]
        except:
            pass
        try:
            self.Bed   = self.Bed[start:end]
        except:
            pass
        try:
            self.Bed_idx = self.Bed_idx[start:end]
        except:
            pass
        try:
            self.Heading = self.Heading[start:end]
        except:
            pass
        try:
            self.Pitch   = self.Pitch[start:end]
        except:
            pass
        try:
            self.Roll    = self.Roll[start:end]
        except:
            pass

        print('==> Clipped along-track: traces {}--{}'.format(start, end))


    ########## END of clip_along() ###########





    ##################################
    # Method: clip data in range
    ##################################

    def clip_range(self, start, end):

        '''
        '''

        self.Data  = self.Data[start:end]
        domain     = self.Domain

        if domain == 'twt':
            self.Time = self.Time[start:end]
            # if not starting at 0, ggf. noch was tun?

        if domain == 'Z':
            self.Z = self.Z[start:end]
        try:
            print('adsfasdf')
        except:
            pass

        try:
            self.Surface_idx = self.Surface_idx - start
        except:
            pass

        try:
            self.Surface_idx_old = self.Surface_idx_old - start
        except:
            pass

        print('==> Clipped in range: bins {}--{}'.format(start, end))


    ########## END of clip_range() ###########




    ##################################
    # Method: concatenate frames
    ##################################

    def concat_frames(added_objects):

        '''

        '''

        import pandas as pd
        import numpy as np
        import copy

        # makes a copy of the first object (serves as a blue print)
        new_obj  = copy.deepcopy(added_objects[0])
        namelist = []

        # default attributes
        Data      = []
        Longitude = []
        Latitude  = []
        Elevation = []
        GPS_time  = []
        Heading   = []
        Roll      = []
        Pitch     = []
        Spacing   = []
        Frames    = []

        # additional important attributes
        if added_objects[0].Domain == 'twt':
            Surface        = [] 
        elif added_objects[0].Domain == 'Z':
            Surface_m      = []


        for obj in added_objects:

            namelist.append(obj.Frame)

            if obj.Domain == 'Z':
                obj.Data.index = obj.Z
                Surface_m.append(obj.Surface_m)
            elif obj.Domain == 'twt':
                obj.Data.index = obj.Time
                Surface.append(obj.Surface)

            Data.append(obj.Data)
            Longitude.append(obj.Longitude)
            Latitude.append(obj.Latitude)
            Elevation.append(obj.Elevation)
            GPS_time.append(obj.GPS_time)
            
            #Bottom.append(obj.Bottom)
            try:
                Heading.append(obj.Heading)
                Roll.append(obj.Roll)
                Pitch.append(obj.Pitch)
            except:
                pass
            try:
                Frames.append(obj.Frame)
            except:
                pass
            try:
                Spacing.append()
            except:
                pass

        new_obj.Reader    = added_objects[0].Reader

        # concatenate object attributes
        new_obj.Data      = pd.concat(Data, axis=1, ignore_index=True)

        if added_objects[0].Domain == 'Z':
            new_obj.Z = new_obj.Data.index
        elif added_objects[0].Domain == 'twt':
            new_obj.Time = new_obj.Data.index

        # delete Z or Time from index
        new_obj.Data.reset_index(inplace=True, drop=True)

        new_obj.Longitude = np.concatenate(Longitude)
        new_obj.Latitude  = np.concatenate(Latitude)
        new_obj.Elevation = np.concatenate(Elevation)
        new_obj.GPS_time  = np.concatenate(GPS_time)

        if obj.Domain == 'twt':
            new_obj.Surface   = np.concatenate(Surface)
            new_obj.get_surf_idx()

        if obj.Domain == 'Z':
            new_obj.Surface_m = np.concatenate(Surface_m)
            new_obj.Data      = new_obj.Data[::-1]
            new_obj.get_surf_m_idx()
            new_obj.Z         = new_obj.Z[::-1]
            


        
        try:
            new_obj.Heading   = np.concatenate(Heading)
            new_obj.Roll      = np.concatenate(Roll)
            new_obj.Pitch     = np.concatenate(Pitch)
        except:
            pass

        
        new_obj.Frames    = Frames
        new_obj.Frame     = added_objects[0].Frame + '_concat'

        print('==> Concatenated {}'.format(namelist))

        return new_obj



    ########## END of concat_frames() ###########





    ##################
    # Method: rename
    ##################

    def rename_frame(self, new_framename):
        self.Frame = new_framename


    ########## END of rename() ###########



    ##################
    # Method: to_dB
    ##################

    def to_dB(self):

        import numpy as np
        import pandas as pd

        if self.dB == False:
            # check if 0 in data because log10 will then return -inf
            # replacing 0 with 1 and log10 will return 0.0
            if 0 in self.Data:
                self.Data = pd.DataFrame(np.where(self.Data==0, 1, self.Data))
                print('... zeroes [0] in self.Data --> replacing with ones [1] before log10.')

            self.Data = 20 * np.log10(self.Data)
            self.dB   = True
            print('==> Converted to dB (20*log10).')

        elif self.dB == True:
            print('... already in dB.')


    ########## END of to_dB() ###########



    #####################
    # Method: inverse_dB
    #####################

    def inverse_dB(self):

        import numpy as np

        if self.dB == True:
            self.Data = 10**(self.Data/20)
            self.dB   = False
            print('==> Converted to 10**(Data/20).')

        elif self.dB == False:
            print('... already NOT in dB.')


    ########## END of inverse_db() ###########




    #############################
    # Method: range gain
    #############################

    '''
    Adds a linear or exponential gain to the radar data
    The function should be uses >!not< in db

    b = base
    n = exponent
    f = (linear) factor

    gain_type = 'linear' (f) or 'exponential' (b) and (n)

    self.Data will be overwritten

    '''

    def range_gain(self, gain_type='', b=2, n=2, f=2):

        from radar_toolbox import add_range_gain

        #if self.dB == True:
        #    self.inverse_dB()
        #    self.dB == False
        #else:
        #    pass
        print('==> adding range gain ({})'.format(gain_type))
        new_data  = add_range_gain(self.Data, gain_type=gain_type, b=b, n=n, f=f)
        self.Data = new_data

        del new_data


    ########## END of range_gain() ###########



    #################################
    # Method: automatic gain control
    #################################

    '''

    '''

    def agc(self, window=50):

        from radar_toolbox import automatic_gain_control

        window = window

        print('==> applying automatic gain control for layer sharpening')
        new_data  = automatic_gain_control(self.Data, window=window)
        self.Data = new_data

        del new_data


    ########## END of magic_gain() ###########



    #############################
    # Method: write shape
    #############################

    def write_shape(self,
                    geometry='Point',
                    step=1,
                    out_filename='',
                    out_folder='',
                    out_format='shapefile',
                    attributes=[]):

        import numpy as np
        import geopandas
        from geo_toolbox import coords2shape
        import copy
        import os

        out_object = copy.deepcopy(self)

        X      = out_object.Longitude
        Y      = out_object.Latitude
        Frame  = out_object.Frame
        Season = 'dummy'

        geometry     = geometry
        step         = step
        out_filename = out_filename
        out_format   = out_format

        if out_filename == '':
            shape_filename = out_object.Frame# + '_' + out_object.Domain

        else:
            shape_filename = out_filename

        if out_format == '':
            out_format = 'shapefile'


        out = coords2shape(X,
                           Y,
                           Frame,
                           Season,
                           EPSG_in=4326,
                           EPSG_out=4326,
                           geometry=geometry,
                           step=step,
                           attributes=attributes)

        if out_folder == '':
            if out_format == 'shapefile':
                if not os.path.exists('shapes'):
                    os.makedirs('shapes')

                out.to_file('shapes/' + shape_filename + '.shp')
                print('==> Written: shapes/{}.shp'.format(shape_filename))

            if out_format == 'geojson':
                if not os.path.exists('geojson'):
                    os.makedirs('geojson')

                out.to_file('geojson/' + shape_filename + '.geojson', driver='GeoJSON')
                print('==> Written: geojson/{}.geojson'.format(shape_filename))

            if out_format == 'kml':
                import fiona
                fiona.supported_drivers['KML'] = 'rw'
                fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'

                if not os.path.exists('geojson'):
                    os.makedirs('geojson')

                out.to_file('kml/' + shape_filename + '.kml', driver='KML')
                print('==> Written: kml/{}.kml'.format(shape_filename))

        else:
            if out_format == 'shapefile':
                out.to_file(out_folder + '/' + shape_filename + '.shp')
                print('==> Written: {}/{}.shp'.format(out_folder, shape_filename))

            if out_format == 'geojson':
                out.to_file(out_folder + '/' + shape_filename + '.geojson', driver='GeoJSON')
                print('==> Written: {}/{}.geojson'.format(out_folder, shape_filename))


            if out_format == 'kml':
                import fiona
                fiona.supported_drivers['KML'] = 'rw'
                fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'

                out.to_file(out_folder + '/' + shape_filename + '.kml', driver='KML')
                print('==> Written: {}/{}.kml'.format(out_folder, shape_filename))

        del X, Y, out

    ########## END of write_mat() ###########





    #############################
    # Method: save matfile
    #############################

    def save_mat(self, out_filename=''):

        import scipy.io
        import pandas as pd
        import numpy as np
        import copy

        try:
            out_object = copy.deepcopy(self)
        except:
            out_object = copy.copy(self)

        try:
            del out_object.Stream
        except:
            pass

        out_object.Data    = out_object.Data.values
        full_dict          = out_object.__dict__
        mat_filename       = out_object.Frame + '_' + out_object.Domain + '.mat'

        if out_filename=='':
            pass
        else:
            mat_filename = out_filename

        try:
            scipy.io.savemat(mat_filename, full_dict)
        except:
            try:
                del out_object.__header__, out_object.__version__, out_object.__globals__
            except:
                pass
            try:
                del out_object.array_param, out_object.param_combine_wf_chan, out_object.param_records
            except:
                pass
            try:
                del out_object.param_csarp, out_object.param_radar
            except:
                pass

            full_dict          = out_object.__dict__
            scipy.io.savemat(mat_filename, full_dict)

        print('==> Written: {}'.format(mat_filename))

        del out_object
        del full_dict
        del mat_filename


    ########## END of write_mat() ###########



    #############################
    # Method: save h5
    #############################


    def save_h5(self, out_filename=''):

        import h5py
        import copy

        try:
            out_object = copy.deepcopy(self)
        except:
            out_object = copy.copy(self)

        try:
            del out_object.Stream
        except:
            pass

        out_object.Data  = out_object.Data.values
        full_dict        = out_object.__dict__
        h5_filename      = out_object.Frame + '_' + out_object.Domain + '.h5'

        if out_filename=='':
            pass
        else:
            h5_filename = out_filename

        
        hf = h5py.File(h5_filename, 'w')


        for key, value in full_dict.items():
                
            try:
                hf.create_dataset(key, data=value)
            except:
                pass

        print('==> Written: {}'.format(h5_filename))
        hf.close()


        del out_object, hf, full_dict, h5_filename

    



    #####################################
    # Converts CReSIS format .mat files
    # to SEGY format
    #####################################

    def to_segy(self, region='', out_filename='', differenciate=False, step=1, save_segy=True, to_dB=False):

        '''
        ==>    Writes Cradar Object as SEGY-Format File.
        Usage: write_segy(self, region='', differenciate=False, step=1)

        - Parameters are:
            1) Cradar Object
            2) region='Greenland' or 'Antarctica'
            3) differenciate=True/False (False is default)
            4) step

                1) Cradar Object, either in 'twt' or 'Z' Domain
                2) if the parameter 'region' is set empty, you already
                   should have 'X' and 'Y' coordinates available for your region
                   --> if you set it to 'Greenland'  they will be projected to EPSG:3413
                   --> if you set it to 'Antarctica' they will be projected to EPSG:3031
                3) choose weather to take the data as it is or to differenciate
                   --> set differenciate=True (default is differenciate=False)
                4) step=N, every N'th trace will be considered (to reduce data size if needed)
        '''

        import numpy as np
        from segy_toolbox import radar2segy

        from obspy import Trace, Stream
        #from obspy.core import AttribDict
        #from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader

        from pyproj import Transformer
        import pdb
        #import time
        #import sys

        step = step

        print('==> Processing Frame: {} located in {}'.format(self.Frame, region))
        print('... This file is in >> {} << domain'.format(self.Domain))

        # check if to_dB is set True, but data is already in dB
        if self.dB and to_dB == True:
            print('... to_dB is set True, but the data is already in dB!')


        # Check on TWT or Elevation Data
        if self.Domain == 'twt':
            domain              = self.Time
            receiver_elevation  = 0
            num_of_samples      = len(self.Time)
            segy_filename       = self.Frame + '_twt.segy'
        elif self.Domain == 'Z':
            domain              = self.Z
            receiver_elevation  = self.Z.max()
            num_of_samples      = len(self.Z)
            segy_filename       = self.Frame + '_Z.segy'

        # the data
        data = self.Data

        # the time
        try:
            gps_time = self.GPS_time
        except:
            gps_time = np.ones(len(self.Longitude))

        # re-project lon, lat to X, Y depending on EPSG
        if region == 'Greenland':
            EPSG = 3413
        elif region == 'Antarctica':
            EPSG = 3031
        else:
            print('define region....')

        transformer = Transformer.from_crs(4326, EPSG, always_xy=True)
        Lon, Lat    = self.Longitude, self.Latitude
        X, Y        = transformer.transform(Lon, Lat)



        # Figure out TWT sampling interval
        diff_domain = np.array([])

        for i in range(len(domain) - 1):
            diff_domain = np.append(diff_domain, domain[i] - domain[i+1])

        # get sample interval | doesn't matter if twt or Z
        sample_interval = np.abs(diff_domain.mean()) * 1000000

        if self.Domain == 'Z':
            sample_interval = 0.01


        # apply radar2segy method
        stream = radar2segy(data=data,
                            receiver_elevation=receiver_elevation,
                            num_of_samples=num_of_samples,
                            sample_interval=sample_interval,
                            X=X,
                            Y=Y,
                            step=step,
                            time_mode='gmtime',
                            gps_time=gps_time,
                            differenciate=differenciate,
                            to_dB=to_dB
                            )

        if out_filename=='':
            pass
        else:
            segy_filename = out_filename

        self.Stream             = stream
        self.Receiver_Elevation = receiver_elevation

        #if to_dB == True:
        #    self.dB = True

        if save_segy == True:
            stream.write(segy_filename, format='SEGY', data_encoding=5, byteorder='>',textual_file_encoding='ASCII')
            print('==> Written: {}'.format(segy_filename))

        elif save_segy == False:
            print('==> Returning: SEGY Stream Object')

        del Lon, Lat, X, Y
        del sample_interval, gps_time
        del stream, data, domain
        del receiver_elevation, num_of_samples


    ########## END of write_mat() ###########


    #####################################
    # Filter
    #
    #####################################

    def filter(raw_object, filter_type='', freq='', freq_min='', freq_max='', corners=2, zerophase=False):

        '''
        This method will apply a filter on an existing Stream object.

            type      =  lowpass
                         highpass
                         bandpass
            freq      =  a frequency as a limit for lowpass or highpass filter
            freq_min  =  lower frequency for bandpass filter
            freq_max  =  upper frequency for bandpass filter
            corners   =  corners
            zerophase =  zerophase
        '''

        import pandas as pd
        import numpy as np
        import copy

        filtered_object = copy.deepcopy(raw_object)

        if filtered_object.Stream:
            pass
        else:
            print('... No Stream Object found.')
            print('==> Consider running to_segy() to create a Stream Object.')

        stream = filtered_object.Stream

        if filter_type == 'highpass':
            print('==> Applying a highpass filter at {} MHz'.format(freq))
            stream_filtered = stream.filter('highpass', freq=freq, corners=corners, zerophase=zerophase)

        elif filter_type == 'lowpass':
            print('==> Applying a lowpass filter at {} MHz'.format(freq))
            stream_filtered = stream.filter('lowpass', freq=freq, corners=corners, zerophase=zerophase)

        elif filter_type == 'bandpass':
            print('==> Applying a bandpass filter between {} - {} MHz'.format(freq_min, freq_max))
            stream_filtered = stream.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=corners, zerophase=zerophase)

        data_filtered = np.array([t.data for t in list(stream_filtered.traces)]).T

        filtered_object.Stream = stream_filtered
        filtered_object.Data   = pd.DataFrame(data_filtered)

        filtered_object.dB = True

        del stream, data_filtered

        return filtered_object






    #####################################
    # Plotting Echograms
    #####################################


    def plot_echogram(self, 
                      figsize_x=10,
                      figsize_y=6,
                      range_mode='twt',
                      every_km_dist=10,
                      every_m_elev=1000,
                      every_twt_ms=10,
                      plot_surface=True,
                      plot_old_surface=False,
                      plot_bed=False,
                      plot_layers=False,
                      show_legend=True,
                      xlabels_as_int=True,
                      ylabels_as_int=True,
                      fontsize=12,
                      show_figure=True, 
                      show_cbar=False,
                      cmap='binary',
                      save_svg=False, 
                      save_png=False, 
                      suffix='',
                      out_folder='',
                      interactive=False,
                      dpi=200):

        '''
        
            range_mode     : 'twt' or 'elevation'
            xlabels_as_int : 
            ylabels_as_int : 
            every_km_dist  :
            every_m_elev   : 
            every_twt_ms   : 


        '''

        import matplotlib.pyplot as plt
        import numpy as np 
        import os

        if interactive == True:
            pass
            #%matplotlib qt

        

        # before plotting check and run these
        try:
            if hasattr(self, 'Surface_idx'):
                pass
            else:
                self.get_surf_idx()
        except:
            pass

        #self.to_dB()
        self.add_distance()


        # Build xticks every n kilometers
        distance_km = self.Distance / 1000
        num_xticks  = int(distance_km[-1]/every_km_dist)

        xticks_km     = np.linspace(0, num_xticks*every_km_dist, num_xticks + 1)
        xticks_km_idx = []

        for i in range(len(xticks_km)):
            idx = np.abs(distance_km - xticks_km[i]).argmin()
            xticks_km_idx.append(idx)
            
        xticks_km_idx = np.array(xticks_km_idx)

        xticks        = xticks_km_idx
        xtick_labels  = xticks_km

        if xlabels_as_int == True:
            xtick_labels  = xticks_km.astype(int)
            
        xaxis_label   = 'Distance (km)'


        # Build yticks every n µs
        if range_mode == 'twt':
            
            twt_ms     = self.Time * 10e5
            num_yticks = int(twt_ms[-1]/every_twt_ms)

            yticks_ms     = np.linspace(0, num_yticks*every_twt_ms, num_yticks + 1)
            yticks_ms_idx = []

            for i in range(len(yticks_ms)):
                idx = np.abs(twt_ms - yticks_ms[i]).argmin()
                yticks_ms_idx.append(idx)

            yticks_ms_idx = np.array(yticks_ms_idx)

            yticks        = yticks_ms_idx
            ytick_labels  = yticks_ms
            
            if ylabels_as_int == True:
                ytick_labels  = yticks_ms.astype(int)
                
            yaxis_label   = 'TWT (µs)'
            

        # Build yticks every n meters
        if range_mode == 'elevation':
            
            elevation_m = self.Z
            below_zero  = int(elevation_m.min() / every_m_elev)
            above_zero  = int(elevation_m.max() / every_m_elev)

            below_zero_steps = np.linspace(below_zero*every_m_elev, 0, np.abs(below_zero) + 1)
            above_zero_steps = np.linspace(every_m_elev, above_zero*every_m_elev, np.abs(above_zero))

            yticks_m     = np.concatenate((below_zero_steps, above_zero_steps))
            yticks_m_idx = []

            for i in range(len(yticks_m)):
                idx = np.abs(elevation_m - yticks_m[i]).argmin()
                yticks_m_idx.append(idx)

            yticks_m_idx = np.array(yticks_m_idx)

            yticks        = yticks_m_idx
            ytick_labels  = yticks_m
            
            if ylabels_as_int == True:
                ytick_labels  = yticks_m.astype(int)
                
            yaxis_label   = 'Elevation (m)'


        # Build yticks every n meters
        if range_mode == 'depth':
            
            depth_m     = self.Depth
            num_yticks = int(depth_m[-1]/every_m_elev)

            yticks_m     = np.linspace(0, num_yticks*every_m_elev, num_yticks + 1)
            yticks_m_idx = []

            for i in range(len(yticks_m)):
                idx = np.abs(depth_m - yticks_m[i]).argmin()
                yticks_m_idx.append(idx)

            yticks_m_idx = np.array(yticks_m_idx)

            yticks        = yticks_m_idx
            ytick_labels  = yticks_m
            
            if ylabels_as_int == True:
                ytick_labels  = yticks_m.astype(int)
                
            yaxis_label   = 'Depth (m)'
            
            

        plt.figure(figsize=(figsize_x,figsize_y))
        
        # plot echogram
        img = plt.imshow(self.Data, aspect='auto', cmap=cmap, alpha=0.8)

        # plot surface ?
        if plot_surface == True:
            if range_mode == 'twt':
                plt.plot(self.Surface_idx, label='surface')
            if range_mode == 'elevation':
                plt.plot(self.Surface_m_idx, label='surface')

        # plot old surface ?
        if plot_old_surface == True:
            if range_mode == 'twt':
                plt.plot(self.Surface_idx_old, color='red', alpha=0.4, label='old surface')
            if range_mode == 'elevation':
                plt.plot(self.Surface_m_idx_old, label='old surface')

        # plot bed ?
        if plot_bed == True:
            if range_mode == 'twt':
                plt.plot(self.Bed_idx, color='red', linewidth=1, linestyle='dashed', label='bed')
            if range_mode == 'elevation':
                plt.plot(self.Bed_m_idx, color='red', linewidth=0.5, linestyle='dashed', label='bed')

        # plot layers ?

        if plot_layers == True:

            layer_list = list(self.Layer.keys())
            n          = len(layer_list)
            colors     = plt.cm.gist_rainbow(np.linspace(0,1,n))
            c          = 0

            if range_mode == 'twt':
                for lr in layer_list:
                    if lr == 'surface':
                        plt.plot(self.Layer[lr]['trace'], self.Layer[lr]['value_twt_idx'], 
                            color='white', linewidth=1, linestyle='dashed')

                    if lr == 'bed':
                        plt.plot(self.Layer[lr]['trace'], self.Layer[lr]['value_twt_idx'], 
                            color='brown', linewidth=0.5, linestyle='dashed')

                    if lr != 'surface':
                        if lr != 'bed':
                            plt.plot(self.Layer[lr]['trace'], self.Layer[lr]['value_twt_idx'], 
                                    color=colors[c], linewidth=0.5)
                            c = c + 1


        plt.xticks(xticks, xtick_labels, fontsize=fontsize)
        plt.xlabel(xaxis_label, fontsize=fontsize)
        plt.yticks(yticks, ytick_labels, fontsize=fontsize)
        plt.ylabel(yaxis_label, fontsize=fontsize)
        plt.title(self.Frame, fontsize=fontsize)

        if show_cbar == True:
            cbr = plt.colorbar(img)
            cbr.set_label('dB')

        if show_legend == True:
            plt.legend()

        if save_png == True:
            if out_folder == '':
                if not os.path.exists('figures'):
                    os.makedirs('figures')

                figname = str(self.Frame) + suffix + '.png'
                plt.savefig('figures/' + figname, dpi=dpi, bbox_inches='tight')
                print('==> Written: figures/{}'.format(figname))
            
            else:
                if not os.path.exists(out_folder):
                    os.makedirs(out_folder)

                figname = str(self.Frame) + suffix + '.png'
                plt.savefig(out_folder + '/' + figname, dpi=dpi, bbox_inches='tight')
                print('==> Written: {}/{}'.format(out_folder, figname))


        if save_svg == True:
            if out_folder == '':
                if not os.path.exists('figures'):
                    os.makedirs('figures')

                figname = str(self.Frame) + suffix + '.svg'
                plt.savefig('figures/' + figname)
                print('==> Written: figures/{}'.format(figname))
            
            else:
                if not os.path.exists(out_folder):
                    os.makedirs(out_folder)

                figname = str(self.Frame) + suffix + '.svg'
                plt.savefig(out_folder + '/' + figname)
                print('==> Written: {}/{}'.format(out_folder, figname))


        if show_figure == True:
            plt.show()
        else:
            plt.clf()
            plt.close('all')

        #del xticks, xlabels, yticks, ylabels



















    def plot_overview(self, flight_lines, save_png=True, dpi=100, out_folder='', cmap='binary', divergent=False, show=True, domain='twt'):

        import matplotlib.pyplot as plt
        from matplotlib import gridspec
        import numpy as np
        import pandas as pd
        import geopandas as gpd
        import os

        # remove old crap
        plt.clf()
        plt.close('all')

        cmap = cmap

        flight_lines = flight_lines
        Lon          = self.Longitude
        Lat          = self.Latitude

        survey_lines = gpd.read_file(flight_lines)

        # generate data frames for points
        df          = pd.DataFrame(Lon)
        df['Lat']   = pd.DataFrame(Lat)
        df.columns  = ['Lon', 'Lat']
        df_first    = df[0:1]

        # create geopandas data frames
        frame  = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Lon, df.Lat))
        first  = gpd.GeoDataFrame(df_first, geometry=gpd.points_from_xy(df_first.Lon, df_first.Lat))

        # set crs to EPSG:4326
        frame        = frame.set_crs(epsg=4326)
        first        = first.set_crs(epsg=4326)
        survey_lines = survey_lines.set_crs(epsg=4326)

        # get distance and spacing
        self.add_distance()

        # define x and y ticks + labels
        xticks  = np.array(range(0, len(self.Distance)))[0::334]
        xlabels = (self.Distance[0::334]/1000).astype(int)

        if domain == 'twt':
            yticks      = np.array(range(0, len(self.Time)))[0::152]
            ylabels     = (self.Time[0::152]*1000000).astype(int)
            yaxis_label = 'TWT (µs)'

        elif domain == 'elevation':
            arr         = self.Z[::500]
            arr         = arr[arr > 0]
            idx         = arr[np.abs(arr).argmin()]
            yticks      = np.array(range(0, len(self.Z)))[idx::500]
            ylabels     = (self.Z[idx::500]).astype(int)
            yaxis_label = 'elevation (m)'


        # Plotting
        fig, ax = plt.subplots(figsize=(30,10))
        gs      = gridspec.GridSpec(1, 2,
                                    width_ratios=[1.2, 3],
                                    height_ratios=[1]
                                    )

        ax0 = plt.subplot(gs[0])
        survey_lines.plot(ax=ax0, color='black', linewidth=0.5, zorder=1)
        frame.plot(ax=ax0, color='blue', markersize=15, zorder=2)
        first.plot(ax=ax0, color='red', markersize=35, zorder=3)
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')

        ax1 = plt.subplot(gs[1])
        if self.dB == False:
            img = plt.imshow(20 * np.log10(self.Data), aspect='auto', cmap=cmap)
        elif self.dB == True:
            if divergent == True:
                data = self.Data
                ende = data.min().min()
                vmin = np.abs(ende) * -1
                vmax = np.abs(ende)
                img  = plt.imshow(data, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
            else:
                img = plt.imshow(self.Data, aspect='auto', cmap=cmap)

        plt.axvline(x=0, color='red', linewidth=5)
        plt.xticks(xticks, xlabels, fontsize=16)
        plt.yticks(yticks, ylabels, fontsize=16)
        plt.xlabel('along-track distance (m)', fontsize=16)
        plt.ylabel(yaxis_label, fontsize=16)
        plt.title(self.Frame, fontsize=16)

        cbr = plt.colorbar(img)
        cbr.set_label('dB', fontsize='16')

        if save_png == True:
            if out_folder == '':
                if not os.path.exists('figures'):
                    os.makedirs('figures')

                figname = str(self.Frame) + '.png'

                plt.savefig('figures/' + figname, dpi=dpi, bbox_inches='tight')
                print('==> Written: figures/{}'.format(figname))
            else:

                figname = str(self.Frame) + '.png'

                plt.savefig(out_folder + figname, dpi=dpi, bbox_inches='tight')
                print('==> Written: {}/{}'.format(out_folder, figname))

        if show == True:
            plt.show()
        else:
            plt.clf()
            plt.close('all')

        del xticks, xlabels, yticks, ylabels
        del frame, first, survey_lines
        del df, flight_lines, Lon, Lat