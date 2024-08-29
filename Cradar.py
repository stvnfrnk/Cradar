

class Cradar:

    #_allObjects = []

    # initiate an instance and read all important parameters
    # from a CReSIS mcords matfile
    def __init__(self):
        #self._allObjects.append(self)
        pass


    #############################
    # Method: load_cresis_mat
    #############################

    def load_cresis_mat(self, filename, dB=False):

        import numpy as np
        from lib.read_input import read_cresis_mat

        Frame, Reader, Domain, Data, Time, Longitude, Latitude, Elevation, GPS_time, Layer = read_cresis_mat(filename)

        # check if Data is transposed
        # if Data.shape[0] == len(Time):
        #     pass
        # else:
        #     Data = np.transpose(Data)

        self.Frame      = Frame
        self.Reader     = Reader
        self.Domain     = Domain
        self.Data       = Data
        self.Time       = Time
        self.Longitude  = Longitude
        self.Latitude   = Latitude
        self.Elevation  = Elevation
        self.GPS_time   = GPS_time
        self.Layer      = Layer
        self.dB         = dB

        return self


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
            
            if 'Surface' in k:
                    self.Layer = {}
                    surface    = {'trace'    : np.arange(len(np.array(self.File['Surface']))).flatten() + 1,
                                'value'  : np.array(self.File['Surface']).flatten()}
                    self.Layer['Surface'] = surface
        

        self.Data      = pd.DataFrame(np.array(self.File['Data']))
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
        print("")
        print('==> Loaded {}'.format(self.Frame))

        return self

    ########## END of load_h5() ###########


    #############################
    # Method: load_awi_segy
    #############################

    def load_awi_segy(self, segy_file='', Longitude='', Latitude='', dB=False, correct_gps=True):

        from lib.read_input import read_awi_segy

        Longitude = Longitude
        Latitude  = Latitude
        Reader    = 'obspy'
        Domain    = 'twt'
        dB        = dB

        # try:
        #     self.Latitude  = np.array( pd.DataFrame(self.Latitude).mask(pd.DataFrame(self.Latitude).duplicated(keep='first'), np.nan).interpolate() ).T[0]
        #     self.Longitude = np.array( pd.DataFrame(self.Longitude).mask(pd.DataFrame(self.Longitude).duplicated(keep='first'), np.nan).interpolate() ).T[0]
        # except:
        #     print('... could not correct gps positions.')

        Data, Time, Frame, stream = read_awi_segy(segy_file)

        self.Frame      = Frame
        self.Reader     = Reader
        self.Domain     = Domain
        self.Data       = Data
        self.Time       = Time
        self.Longitude  = Longitude
        self.Latitude   = Latitude
        # self.Elevation  = Elevation
        # self.GPS_time   = GPS_time
        # self.Layer      = Layer
        self.dB         = dB
        self.Stream     = stream

        return self


    ########## END of load_awi_segy() ###########


    #############################
    # Method: load_awi_nc
    #############################


    def load_awi_nc(self, nc_file='', read_agc=False):

        '''


        '''

        from lib.read_input import read_awi_nc

        
        Data, Time, Longitude, Latitude, Aircraft_altitude, Ice_surface_elevation, Layer = read_awi_nc(nc_file, read_agc=read_agc)

        #self.Frame      = Frame
        self.Reader     = 'xarray'
        self.Domain     = 'twt'
        self.Data       = Data
        self.Time       = Time
        self.Longitude  = Longitude
        self.Latitude   = Latitude

        # self.Aircraft_altitude  = Aircraft_altitude
        # self.GPS_time   = GPS_time
        self.Layer      = Layer
        self.dB         = True

        print("")
        print('==> Loaded {}'.format(nc_file))

        return self

    ########## END of load_awi_nc() ###########




    #############################
    # Method: load_awi_nc
    #############################


    def load_awi_nc_old(self, nc_file='', read_agc=False):

        '''


        '''

        from lib.read_input import read_awi_nc_old

        
        Data, Time, Longitude, Latitude, Aircraft_altitude, Ice_surface_elevation, Layer = read_awi_nc_old(nc_file)

        #self.Frame      = Frame
        self.Reader     = 'xarray'
        self.Domain     = 'twt'
        self.Data       = Data
        self.Time       = Time
        self.Longitude  = Longitude
        self.Latitude   = Latitude

        # self.Aircraft_altitude  = Aircraft_altitude
        # self.GPS_time   = GPS_time
        self.Layer      = Layer
        self.dB         = True

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
            try:
                self.Data        = pd.DataFrame(ds.variables['pulse_data'].values)
            except:
                pass
            try:
                self.Data        = pd.DataFrame(ds.variables['polarised_pulse_SPHV_data'].values)
            except:
                pass
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
            print('==> Loaded {} >> {} data <<'.format(self.Frame, data_type))

            del ds

            return self

        if data_type == 'chirp':

            ds = xr.open_dataset(filename, decode_times=False)

            self.Frame       = os.path.split(filename)[1].split('.nc')[0]
            self.Time        = ds['fast_time'].values.astype(float) * 10e-7
            try:
                self.Data        = pd.DataFrame(ds.variables['chirp_data'].values)
            except:
                pass
            try:
                self.Data        = pd.DataFrame(ds.variables['polarised_chirp_SSHH_data'].values)
            except:
                pass
            try:
                self.Data        = pd.DataFrame(ds.variables['chirp_cHG_data'].values)
            except:
                pass

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


            dist            = self.Elevation - self.Surface_m
            dist[dist == 0] = np.nan
            dist            = np.array(pd.DataFrame(dist).interpolate(limit_direction='both'))
            twt_s           = dist / (299792458) * 4

            print('')
            print('==> Loaded {} >> {} data <<'.format(self.Frame, data_type))

            del ds

            

            ### load Layers

            self.Layer = {}
            
            surface    = {'trace'  : np.arange(len(np.array(twt_s).flatten())) + 1,
                        'value'    : np.array(twt_s).flatten(),
                        'color'    : [255, 255, 9]}
            self.Layer['Surface'] = surface

            # bottom    = {'trace'  : np.arange(len(np.array(self.Bed_idx).flatten())) + 1,
            #                 'value'  : np.array(self.Time[self.Surface_idx.astype(int)]).flatten(),
            #                 'color'  : 'red'}
            # self.Layer['Bottom'] = bottom

            return self

            # from lib.bas_io import correct_chirp_data

            # data_chirp, time, longitude, latitude, elevation, gps_time, surface_idx, surface_m, bed_idx, bed_m = correct_chirp_data(filename)

            # self.Frame       = os.path.split(filename)[1].split('.nc')[0]
            # self.Time        = time
            # self.Data        = data_chirp
            # self.Domain      = 'twt'
            # self.dB          = False

            # self.Longitude   = longitude
            # self.Latitude    = latitude
            # self.Elevation   = elevation
            # self.GPS_time    = gps_time

            # self.Surface_idx = surface_idx
            # self.Surface_m   = surface_m
            # self.Bed_idx     = bed_idx
            # self.Bed_m       = bed_m

            # print('')
            # print('==> Loaded {} >> chirp data <<'.format(self.Frame))

            # del data_chirp, time, longitude, latitude, elevation, gps_time, surface_idx, surface_m, bed_idx, bed_m

            # return self

    ########## END of load_bas_nc() ###########
        


    #############################
    # Method: load_custom
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

    def add_layer_by_frame_trace(self, layer_name, traces, values, color=''):

        import numpy as np
        import pandas as pd

        # if layer_name == 'Surface':
        #     color = [255, 255, 9]

        # if layer_name == 'Bed':
        #     color = [227, 26, 28]

        color = tuple(np.round(np.array( color ) / 255, 3))

        df_idx       = pd.DataFrame(np.arange(len(self.Longitude)) + 1)
        df_idx.index = df_idx.index + 1

        df_lyr = pd.DataFrame(traces)
        df_lyr["value"] = values

        df_lyr.columns = ["trace", "value"]
        df_lyr.index   = df_lyr["trace"]
        del df_lyr["trace"]

        df     = df_idx.join(df_lyr)
        values = df["value"].values
        traces = df.index.values

        layer = {'trace'  : traces,
                 'value'  : values,
                 'color'  : color}
        
        try:
            self.Layer
        except:
            self.Layer             = {}
        
        self.Layer[layer_name] = layer
        #self.reshape_layer_values(layer_name)

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

        layer = {'trace'  : traces,
                 'value'  : values}

        self.Layer[layer_name] = layer
        print('==> added layer: {}'.format(layer_name))



    def fix_missing_surface_picks(self):

        import numpy as np

        missing = len(self.Longitude) - len(self.Layer['Surface']['trace'])
        if missing == 0:
            pass
        elif missing > 0:
            for i in range(missing):
                self.Layer['Surface']['trace'] = np.append(self.Layer['Surface']['trace'], self.Layer['Surface']['trace'][-1] + 1)
                self.Layer['Surface']['value'] = np.append(self.Layer['Surface']['value'], self.Layer['Surface']['value'][-1])
        elif missing < 0:
            self.Layer['Surface']['trace'] = self.Layer['Surface']['trace'][:missing]
            self.Layer['Surface']['value'] = self.Layer['Surface']['value'][:missing]


    def get_layer_idx(self):

        import numpy as np

        layer_list = list(self.Layer.keys())
        time        = self.Time

        for lr in layer_list:
            values_idx = np.array([])
            v_twt      = self.Layer[lr]['value']

            for i in range(len(v_twt)):
                if v_twt[i] == np.nan:
                    v_idx = np.nan
                else:
                    v_idx      = (np.abs(np.array(time) - np.array(v_twt)[i])).argmin()
                    values_idx = np.append(values_idx, v_idx)
            
            #values_idx[values_idx == 0] = np.nan

            self.Layer[lr]['value_idx'] = values_idx.astype(int)
            print('... getting layer idx for {}'.format(lr))


    def reshape_layer_values(self, layer_name):

        import numpy as np

        new_values = []
        for x in range(len(self.Longitude)):
            try:
                if x + 1 in self.Layer[layer_name]['trace']:
                    new_values.append(self.Layer[layer_name]['value_idx'][x])
                else:
                    new_values.append(np.nan)
            except:
                new_values.append(np.nan)

        self.Layer[layer_name]['trace_']     = np.arange(len(self.Longitude)) + 1
        self.Layer[layer_name]['value_idx_'] = np.array(new_values).flatten()


    def get_layer_m_idx(self):

        import numpy as np

        layer_list = list(self.Layer.keys())
        Z          = self.Z

        for lr in layer_list:
            if '_m' in lr:
                values_idx = np.array([])
                vz         = self.Layer[lr]['value']
                for i in range(len(vz)):
                    if vz[i] == np.nan:
                        v_idx = np.nan
                    else:
                        try:
                            v_idx      = (np.abs(np.array(Z) - int(np.array(vz)[i]))).argmin()
                        except:
                            v_idx      = np.nan

                        values_idx = np.append(values_idx, v_idx)
                
                #values_idx[values_idx == 0] = np.nan

                self.Layer[lr]['value_idx'] = values_idx.astype(int)
                print('... getting layer idx for {}'.format(lr))

    
    def layer2surface(self):

        import numpy as np

        layer_list = list(self.Layer.keys())

        for lr in layer_list:
            if '_m' not in lr:
                if 'Surface' not in lr:
                    common_traces = np.intersect1d(self.Layer[lr]['trace'], self.Layer['Surface']['trace'])
                    new_vals = []
                    color = self.Layer[lr]['color']
                    for tr in common_traces:
                        idx_hr   = np.where(self.Layer[lr]['trace'] == tr)[0][0]
                        idx_sf   = np.where(self.Layer['Surface']['trace'] == tr)[0][0]
                        val      = self.Layer[lr]['value'][idx_hr] - self.Layer['Surface']['value'][idx_sf]
                        new_vals.append(val)

                    self.Layer[lr]['trace'] = np.array(common_traces)
                    self.Layer[lr]['value'] = np.array(new_vals)
                    
                    
                    try:
                        common_traces = np.intersect1d(self.Layer[lr]['trace'], self.Layer['Surface']['trace'])
                        new_vals_idx = []
                        color = self.Layer[lr]['color']
                        for tr in common_traces:
                            idx_hr   = np.where(self.Layer[lr]['trace'] == tr)[0][0]
                            idx_sf   = np.where(self.Layer['Surface']['trace'] == tr)[0][0]
                            val_idx  = self.Layer[lr]['value_idx'][idx_hr] - self.Layer['Surface']['value_idx'][idx_sf]
                            new_vals_idx.append(val_idx)

                        self.Layer[lr]['value_idx'] = np.array(new_vals_idx)
                    except:
                        pass

                    print('... layer2surface for {}'.format(lr))




    def layer2elevation(self, speed_of_ice=168900000.0):

        import numpy as np

        speed_of_ice = speed_of_ice

        layer_list = list(self.Layer.keys())

        for lr in layer_list:
            if '_m' not in lr:
                if 'Surface' not in lr:
                    common_traces = np.intersect1d(self.Layer[lr]['trace'], self.Layer['Surface']['trace'])
                    new_vals = []
                    color = self.Layer[lr]['color']
                    for tr in common_traces:
                        idx_hr   = np.where(self.Layer[lr]['trace'] == tr)[0][0]
                        idx_sf   = np.where(self.Layer['Surface']['trace'] == tr)[0][0]
                        val_m     = self.Layer['Surface_m']['value'][idx_sf] - ( (self.Layer[lr]['value'][idx_hr] - self.Layer['Surface']['value'][idx_sf]) * 168900000.0 / 2 )
                        new_vals.append(val_m)

                    new_vals     = np.array(new_vals)
                    new_lyr_name = lr + '_m' 
                    new_lyr      = layer = {'trace'  : common_traces,
                                            'value'  : new_vals,
                                            'color'  : color}
                    self.Layer[new_lyr_name] = new_lyr

                    print('... layer2elevation for {}'.format(lr))


    # def get_layer_m_idx():

        # import numpy as np

        # Z          = self.Z
        # layer_list = list(self.Layer.keys())
        
        # for lr in layer_list:
        #     Z_layer




        # Z_surf     = self.Surface_m
        # surf_m_idx = np.array([])

        # for i in range(len(Z_surf)):
        #     idx        = len(Z) - ( (np.abs(np.array(Z) - np.array(Z_surf)[i])).argmin() )
        #     surf_m_idx = np.append(surf_m_idx, idx)

        # self.Surface_m_idx = surf_m_idx
        # print('==> Added pixel index of surface elevation')
        # del Z, Z_surf, surf_m_idx


    def track_surface(self, skip=100, llim='', gauss_factor=1, use_gradient=False, offset=2):

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

        if llim == '':
            surf_idx       = ( data_filtered[skip::].idxmax(axis=0, skipna=True) ).astype(int)
        else:
            surf_idx       = ( data_filtered[skip:llim].idxmax(axis=0, skipna=True) ).astype(int)

        # get surface values
        srf = []
        
        for i in range(len(surf_idx)):
            index = surf_idx[i]
            val   = self.Time[index]
            srf.append(val)
        
        self.Layer = {}
        surface    = {'trace'   : np.arange(len(np.array(srf))) + 1,
                    'value'     : np.array(srf),
                    'color'     : 'yellow'}
        
        self.Layer['Surface'] = surface
        
        self.get_layer_idx()

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

        import numpy as np

        geotif      = geotif
        geotif_name = geotif_name
        Longitude   = self.Longitude
        Latitude    = self.Latitude

        print('==> Applying gridtrack method 1 ...')
        from lib.geo_toolbox import gridtrack
        raster_vals = gridtrack(Longitude=Longitude, Latitude=Latitude, geotif=geotif, geotif_name=geotif_name, geotif_epsg=geotif_epsg)

        #print('==> Failed, trying gridtrack method 2 ...')
        #from geo_toolbox import gridtrack2
        #raster_vals = gridtrack(geotif, Latitude, Longitude, EPSG_xy=4326, EPSG_raster=geotif_epsg)
        setattr(self, geotif_name, raster_vals)
        surface_m    = {'trace'  : np.arange(len(np.array(raster_vals))).flatten() + 1,
                        'value'  : np.array(raster_vals).flatten(),
                        'color'  : self.Layer['Surface']['color']}
        self.Layer['Surface_m'] = surface_m
        
        print('==> Added {} to the data'.format(geotif_name))

        del raster_vals, Longitude, Latitude




    ########## END of rename() ###########



    #############################
    # Method: retrack_surf()
    #############################

    def retrack_surface(self, roll_factor=10, sigma=2, gate=[0, 0], differenciate=False, offset=1):

        import copy
        import numpy as np
        import pandas as pd
        from scipy.ndimage import gaussian_filter

        roll_factor = roll_factor
        sigma       = sigma
        ulim        = gate[0]
        llim        = gate[1]
        offset      = offset

        data_filtered = pd.DataFrame(gaussian_filter(self.Data.values, sigma=2))
        if differenciate == True:
            data_filtered = data_filtered.diff().rolling(roll_factor, center=True, win_type='hamming').mean()
        else:
            data_filtered = data_filtered.rolling(roll_factor, center=True, win_type='hamming').mean()

        # get surface values
        srf         = []
        srf_idx     = []
        surface_trc = self.Layer['Surface']['trace']
        surface_val = self.Layer['Surface']['value']
        surface_idx = self.Layer['Surface']['value_idx']
        surface_col = self.Layer['Surface']['color']

        for i in range(len(surface_idx)):
            trx         = np.array(data_filtered[i])[int(surface_idx[i]) - int(ulim):int(surface_idx[i]) + int(llim)]
            offset_gate = np.argmax(trx) - int(((llim + ulim)/2))
            index       = surface_idx[i] + offset_gate + offset
            val         = self.Time[int(index)]
            srf.append(val)
            srf_idx.append(index)

        del self.Layer['Surface']

        self.Layer['Surface']              = {}
        self.Layer['Surface']['trace']     = surface_trc
        self.Layer['Surface']['value']     = np.array(srf)
        self.Layer['Surface']['value_idx'] = np.array(srf_idx)
        self.Layer['Surface']['color']     = surface_col

        self.Layer['Surface_old']              = {}
        self.Layer['Surface_old']['trace']     = surface_trc
        self.Layer['Surface_old']['value']     = surface_val
        self.Layer['Surface_old']['value_idx'] = surface_idx
        self.Layer['Surface_old']['color']     = (0.95, 0.95, 0.95)
 
        print('==> Retracked ice Surface')

        #return self.Surface_idx_old




    #############################
    # Method: correct_geom_spreading
    #############################

    def correct4attenuation(raw_object, mode=0, loss_factor=0):

        from lib.radar_toolbox import correct4attenuation
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

        # geom_obj.get_surf_idx()

        data     = geom_obj.Data
        twt      = geom_obj.Time
        surf_idx = geom_obj.Layer['Surface']['value_idx']

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

        from lib.geo_toolbox import coords2distance

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

    def twt2elevation( twt_object,
                       reference='',
                       sample_int=1,
                       speed_of_ice=1.689e8):

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
        from lib.radar_toolbox import radar_twt2elevation
        import copy

        # makes a copy of the first object (serves as a blue print)
        try:
            elev_obj = copy.deepcopy(twt_object)
        except:
            elev_obj = copy.copy(twt_object)

        print('==> Now: twt2elevation...')

        if reference == 'GPS':
            print('... Using aircraft GPS and radar Surface reflection to derive elevation')

        elif reference == 'DEM':
            print('... Using a ice Surface DEM to derive elevation')

        elif reference == '':
            reference = 'GPS'
            print("... !! You didn't define a reference, it is now automatically set to 'GPS'")
            print('... Using aircraft GPS and radar Surface reflection to derive elevation')

        # Get ice Surface elevation values from DEM
        if reference == 'DEM':
            twt_surface = elev_obj.Layer['Surface']['value']
            dem_surf = elev_obj.DEM_surface
            aircraft_elevation = np.ones(len(elev_obj.Longitude))

        elif reference == 'GPS':
            twt_surface = elev_obj.Surface
            dem_surf = ''
            aircraft_elevation = elev_obj.Elevation

        # imput variables from instance
        data               = elev_obj.Data
        twt                = elev_obj.Time
        
        # input variables defined in prior steps
        twt_surface   = twt_surface

        # input variables from input of method above
        reference = reference
        sample_int   = sample_int

        surf_idx = twt_object.Layer["Surface"]["value_idx"]

        data_elev, Z, surf_m_idx =  radar_twt2elevation(data=data,
                                                        twt=twt,
                                                        surf_idx=surf_idx,
                                                        speed_of_ice=speed_of_ice,
                                                        dem_surf=dem_surf,
                                                        sample_int=sample_int
                                                        )


        # re-define instance atributes
        elev_obj.Z             = Z
        elev_obj.Data          = data_elev

        surface_m    = {'trace'     : elev_obj.Layer['Surface']['trace'],
                        'value'     : dem_surf,
                        'value_idx' : surf_m_idx,
                        'color'     : elev_obj.Layer['Surface']['color']}
        
        elev_obj.Layer['Surface_m'] = surface_m


        elev_obj.Domain          = 'Z'
        elev_obj.Sample_Interval = "{} m".format(sample_int)

        del data_elev
        
        return elev_obj

        


    ########## END of twt2elevation() ###########



#############################
    # Method: twt2surface
#############################

    def twt2surface(twt_object):

        '''
        
        '''

        import numpy as np
        import copy
        from lib.radar_toolbox import radar_twt2surface

        # makes a copy of the first object (serves as a blue print)
        try:
            p2s_obj = copy.deepcopy(twt_object)
        except:
            p2s_obj = copy.copy(twt_object)

        print('==> Now: pull2elevation...')

        # imput variables from instance
        data     = p2s_obj.Data
        twt      = p2s_obj.Time
        surf_idx = p2s_obj.Layer['Surface']['value_idx'] #p2s_obj.Surface

        data_new, time_new =  radar_twt2surface(data=data, twt=twt, surf_idx=surf_idx)


        # re-define instance atributes
        p2s_obj.Time   = time_new
        p2s_obj.Data   = data_new
        p2s_obj.Domain = 'twt'
        p2s_obj.Layer['Surface']['value'] = np.repeat(0, len(surf_idx))

        del data, twt, surf_idx

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
        from lib.radar_toolbox import radar_pull2bed
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
        try:
            self.GPS_time  = self.GPS_time[::-1]
            self.Elevation = self.Elevation[::-1]
        except:
            pass

        if self.Domain == 'twt':
            try:
                layer_list = list(self.Layer.keys())
                for lr in layer_list:
                    # self.Layer[lr]['trace']     = self.Layer[lr]['trace'][::-1]
                    self.Layer[lr]['value']     = self.Layer[lr]['value'][::-1]
                    try:
                        self.Layer[lr]['value_idx']     = self.Layer[lr]['value_idx'][::-1]
                    except:
                        pass
            except:
                pass

        elif self.Domain == 'Z':
            try:
                self.Surface_m     = self.Surface_m[::-1]
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

    def clip_along(self, start_val, end_val, mode='trace'):

        '''
        '''

        import pandas as pd
        import numpy as np

        if mode == 'trace':
            start = start_val 
            end   = end_val

        if mode == 'distance':
            start = (np.abs(self.Distance - start_val)).argmin()
            end   = (np.abs(self.Distance - end_val)).argmin()

        self.Data      = self.Data[:,start:end]
        self.Longitude = self.Longitude[start:end]
        self.Latitude  = self.Latitude[start:end]
        
        # optional
        #try:
        layer_list = list(self.Layer.keys())
        for lr in layer_list:
            self.Layer[lr]["value"] = self.Layer[lr]["value"][start:end]
            self.Layer[lr]['trace'] = np.arange(int(end - start)).flatten() + 1
                
                
                # layer_keys = list(self.Layer[lr].keys())
                # for key in layer_keys:
                #     if "color" in key:
                #         pass
                #     else:
                #         self.Layer[lr][key]     = self.Layer[lr][key][start:end]
                #         self.Layer[lr]['trace'] = np.arange(len(self.Layer[lr]['trace'])).flatten() + 1
        #except:
        #    pass

        try:
            self.Elevation = self.Elevation[start:end]
        except:
            pass
        try:
            self.GPS_time  = self.GPS_time[start:end]
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
        try:
            self.Spacing    = self.Spacing[start:end]
        except:
            pass
        try:
            self.Distance    = self.Distance[start:end]
        except:
            pass

        print('==> Clipped along-track: traces {}--{}'.format(start, end))


    ########## END of clip_along() ###########





    ##################################
    # Method: clip data in range
    ##################################

    def clip_range(self, start_val, end_val, mode='index'):

        import numpy as np

        '''

        mode = 'index'
        mode = 'elevation'
        mode = 'twt'
        '''

        if mode == 'index':
            start = start_val
            end   = end_val

        if mode == 'elevation':
            lower_lim = np.where(self.Z == start_val)[0][0]
            upper_lim = np.where(self.Z == end_val)[0][0]

            if lower_lim > upper_lim:
                start = upper_lim
                end   = lower_lim
            else:
                start = lower_lim
                end   = upper_lim

        if mode == 'twt':
            start_val = start_val / 10**9
            end_val   = end_val   / 10**9

            start = (np.abs(self.Time - start_val)).argmin()
            end   = (np.abs(self.Time - end_val)).argmin()

            # if lower_lim > upper_lim:
            #     start = upper_lim
            #     end   = lower_lim
            # else:
            #     start = lower_lim
            #     end   = upper_lim
        
        self.Data  = self.Data[start:end,:]
        domain     = self.Domain

        if domain == 'twt':
            t_start   = self.Time[start]
            self.Time = self.Time[start:end]
            self.Time = self.Time - t_start
            
            # if not starting at 0, ggf. noch was tun?

        if domain == 'Z':
            self.Z = self.Z[start:end]
        
        try:
            layer_list = list(self.Layer.keys())
            for lr in layer_list:
                self.Layer[lr]['value']     = self.Layer[lr]['value'] - t_start
                self.Layer[lr]['value_idx'] = self.Layer[lr]['value_idx'] - start
        except:
            pass
        # try:
        #     print('')
        # except:
        #     pass

        # try:
        #     self.Surface_idx = self.Surface_idx - start
        # except:
        #     pass

        # try:
        #     self.Surface_idx_old = self.Surface_idx_old - start
        # except:
        #     pass

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

            # if obj.Domain == 'Z':
            #     obj.Data.index = obj.Z
            #     Surface_m.append(obj.Surface_m)
            # elif obj.Domain == 'twt':
            #     obj.Data.index = obj.Time
            #     Surface.append(obj.Surface)

            Data.append(obj.Data)
            Longitude.append(obj.Longitude)
            Latitude.append(obj.Latitude)
            
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
            try:
                Elevation.append(obj.Elevation)
            except:
                pass
            try:
                GPS_time.append(obj.GPS_time)
            except:
                pass

        new_obj.Reader    = added_objects[0].Reader

        # concatenate object attributes
        # new_obj.Data      = pd.concat(Data, axis=1, ignore_index=True)
        new_obj.Data      = np.concatenate(Data, axis=1)

        if added_objects[0].Domain == 'Z':
            new_obj.Z = new_obj.Data.index
        elif added_objects[0].Domain == 'twt':
            new_obj.Time = new_obj.Time

        # delete Z or Time from index
        # new_obj.Data.reset_index(inplace=True, drop=True)

        new_obj.Longitude = np.concatenate(Longitude)
        new_obj.Latitude  = np.concatenate(Latitude)

        try:
            new_obj.Heading   = np.concatenate(Heading)
            new_obj.Roll      = np.concatenate(Roll)
            new_obj.Pitch     = np.concatenate(Pitch)
            new_obj.Elevation = np.concatenate(Elevation)
            new_obj.GPS_time  = np.concatenate(GPS_time)
        except:
            pass

        ###########################
        # concat layer

        crd_list   = added_objects
        
        lr_key_list = []
        for crd in crd_list:
            for ky in crd.Layer.keys():
                lr_key_list.append(ky)
            
        layer_list = np.unique(lr_key_list)
        
        # create empty Layer dicts
        for lr in layer_list:
            lr_dct   = {"trace"  : 0,
                        "value"  : 0,
                        "color"  : 0}
            new_obj.Layer[lr] = lr_dct

        for lr in layer_list:
            print("==> Concatenating Layer: {}".format(lr))
            
            # find first frame where layer appears
            nn = []
            n = 0
            for crd in crd_list:
                if lr in crd.Layer:
                    # print("SnAcc_02 at n={}".format(n))
                    nn.append(n)
                n = n + 1
            k = min(nn)
            
            if k == 0:
                n_traces      = [crd_list[k].Layer[lr]['trace']]
                n_values      = [crd_list[k].Layer[lr]['value']]
                color         = [crd_list[k].Layer[lr]['color']]
            else:
                pass
                
            concat_length = []

            for i in np.arange(len(crd_list) - 1) + 1:
                prev_frame_length = len(crd_list[i-1].Longitude)
                concat_length.append(prev_frame_length)
                
                if k == 0:
                    try:
                        n_traces.append(crd_list[i].Layer[lr]['trace'] + np.sum(concat_length))
                        n_values.append(crd_list[i].Layer[lr]['value'])
                    except:
                        pass
                else:
                    if i < k:
                        pass
                    
                    if i == k:
                        n_traces      = [crd_list[k].Layer[lr]['trace'] + np.sum(concat_length)]
                        n_values      = [crd_list[k].Layer[lr]['value']]
                        color         = [crd_list[k].Layer[lr]['color']]
                        
                    if i > k:
                        try:
                            n_traces.append(crd_list[i].Layer[lr]['trace'] + np.sum(concat_length))
                            n_values.append(crd_list[i].Layer[lr]['value'])
                        except:
                            pass
                
                    

            new_traces = np.concatenate(n_traces)
            new_values = np.concatenate(n_values)

            new_obj.Layer[lr]['trace'] = new_traces
            new_obj.Layer[lr]['value'] = new_values
            new_obj.Layer[lr]['color'] = color

                
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
                self.Data = np.where(self.Data==0, 1, self.Data)
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

        from lib.radar_toolbox import add_range_gain

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

        from lib.radar_toolbox import automatic_gain_control


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
        from lib.geo_toolbox import coords2shape
        import copy
        import os

        self.add_distance()

        out_object = copy.copy(self)

        X       = out_object.Longitude
        Y       = out_object.Latitude
        Frame   = out_object.Frame
        Segment = out_object.Segment
        Season  = out_object.Season

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
                           Segment,
                           Season,
                           EPSG_in=4326,
                           EPSG_out=4326,
                           geometry=geometry,
                           step=step,
                           attributes=attributes)
        
        if geometry == "Point":
            out['Distance_m'] = self.Distance.astype(int)
            out['Trace']      = np.arange(len(self.Longitude)) + 1
        else:
            pass

        if geometry == "Point":
            shape_filename = shape_filename + "_point"

        if geometry == "Line":
            shape_filename = shape_filename + "_line"

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
    # Method: write matfile
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

        out_object.Data    = out_object.Data
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

    def to_segy(self, epsg='', out_filename='', differenciate=False, step=1, save_segy=True, to_dB=False):

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
        from lib.segy_toolbox import radar2segy

        from obspy import Trace, Stream
        #from obspy.core import AttribDict
        #from obspy.io.segy.segy import SEGYTraceHeader, SEGYBinaryFileHeader

        from pyproj import Transformer
        import pdb
        #import time
        #import sys

        step = step
        epsg = epsg

        print('==> Processing Frame: {} located in EPSG{}'.format(self.Frame, epsg))
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
            segy_filename       = self.Frame + '_' + str(receiver_elevation) + '.segy'

        # the data
        data = self.Data

        # the time
        try:
            gps_time = self.GPS_time
        except:
            gps_time = np.ones(len(self.Longitude))


        transformer = Transformer.from_crs(4326, epsg, always_xy=True)
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
                      range_mode='',
                      every_km_dist=10,
                      every_m_elev=1000,
                      every_twt=['ms', 10],
                      plot_layers=False,
                      markersize=0.2,
                      show_legend=True,
                      xlabels_as_int=True,
                      ylabels_as_int=True,
                      vline='',
                      fontsize=12,
                      show_figure=True, 
                      show_cbar=True,
                      cmap='binary',
                      vmin='',
                      vmax='',
                      save_svg=False, 
                      save_raster=[False, "jpg"], 
                      suffix='',
                      out_folder='',
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

            if every_twt[0] == 'ns':
                twt        = self.Time * 10e8
            if every_twt[0] == 'ms':
                twt        = self.Time * 10e5


            num_yticks = int(twt[-1]/every_twt[1])

            yticks_ms     = np.linspace(0, num_yticks*every_twt[1], num_yticks + 1)
            yticks_ms_idx = []

            for i in range(len(yticks_ms)):
                idx = np.abs(twt - yticks_ms[i]).argmin()
                yticks_ms_idx.append(idx)

            yticks_ms_idx = np.array(yticks_ms_idx)

            yticks        = yticks_ms_idx
            ytick_labels  = yticks_ms
            
            if ylabels_as_int == True:
                ytick_labels  = yticks_ms.astype(int)
                
            if every_twt[0] == 'ns':
                yaxis_label   = 'TWT (ns)'
            if every_twt[0] == 'ms':
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
            

        if vmin == '':
            vmin = np.nanmin(self.Data)
        else:
            vmin = vmin

        if vmax == '':
            vmax = np.nanmax(self.Data)
        else:
            vmax = vmax    
            

        plt.figure(figsize=(figsize_x,figsize_y))
        
        # plot echogram
        img = plt.imshow(self.Data, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax, alpha=0.8)

        # plot layers
        if plot_layers == True:
            layer_list = list(self.Layer.keys())
            n          = len(layer_list) #- 2
            offset     = 0
            colors     = plt.cm.spring(np.linspace(0,1,n + offset))
            c = 0

            if range_mode == 'twt':
                for lr in layer_list:
                    self.reshape_layer_values(lr)

                for lr in layer_list:    
                    if lr == 'Surface':
                        plt.scatter(x=self.Layer[lr]['trace'], y=self.Layer[lr]['value_idx'], 
                            s=markersize, label=lr)

                    elif lr == 'Base':
                        plt.scatter(x=self.Layer[lr]['trace'], y=self.Layer[lr]['value_idx'], 
                            color=self.Layer[lr]['color'], s=markersize, label=lr)

                    elif lr != 'Surface' and lr != 'Base' and lr != 'Surface_m':
                        plt.scatter(x=self.Layer[lr]['trace'], y=self.Layer[lr]['value_idx'], 
                                color=self.Layer[lr]['color'], s=markersize, label=lr)

            elif range_mode == 'elevation':
                for lr in layer_list:
                    if '_m' not in lr:
                        pass
                    else:
                        if lr == 'Surface_m':
                            plt.scatter(x=self.Layer[lr]['trace'], y=self.Layer[lr]['value_idx'], 
                                color=self.Layer[lr]['color'], s=markersize, label=lr)

                        elif lr == 'Bed_m':
                            plt.scatter(x=self.Layer[lr]['trace'], y=self.Layer[lr]['value_idx'], 
                                color=self.Layer[lr]['color'], s=markersize, linestyle='dashed', label=lr)

                        elif lr != 'Surface_m':
                            if lr != 'Bed_m':
                                if "_m" in lr:
                                    plt.scatter(x=self.Layer[lr]['trace'], y=self.Layer[lr]['value_idx'], 
                                            color=self.Layer[lr]['color'], s=markersize, label=lr)

        if vline != '':
            print('vlines')
            plt.vlines(x=vline, ymin=0, ymax=len(self.Time), color='white')

        plt.xticks(xticks, xtick_labels, fontsize=fontsize)
        plt.xlabel(xaxis_label, fontsize=fontsize)
        plt.yticks(yticks, ytick_labels, fontsize=fontsize)
        plt.ylabel(yaxis_label, fontsize=fontsize)
        plt.xlim(0, len(self.Longitude))
        plt.ylim(self.Data.shape[0], 0)
        plt.title(self.Frame, fontsize=fontsize)



        if show_cbar == True:
            cbr = plt.colorbar(img)
            cbr.set_label('dB')

        if show_legend == True:
            plt.legend()

        if save_raster[0] == True:
            if save_raster[1] == "jpg":
                raster_type = ".jpg"
            if save_raster[1] == "png":
                raster_type = ".png" 

            if out_folder == '':
                if not os.path.exists('figures'):
                    os.makedirs('figures')

                figname = str(self.Frame) + suffix + raster_type
                plt.savefig('figures/' + figname, dpi=dpi, bbox_inches='tight')
                print('==> Written: figures/{}'.format(figname))
            
            else:
                if not os.path.exists(out_folder):
                    os.makedirs(out_folder)

                figname = str(self.Frame) + suffix + raster_type
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

                figname = str(self.Frame) + '.jpg'

                plt.savefig('figures/' + figname, dpi=dpi, bbox_inches='tight')
                print('==> Written: figures/{}'.format(figname))
            else:

                figname = str(self.Frame) + '.jpg'

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