

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

    def load_cresis_mat(self, filename, frame="", dB=False):

        import numpy as np
        from lib.read_input import read_cresis_mat

        Frame, Reader, Domain, Data, Time, Longitude, Latitude, Elevation, GPS_time, Layer = read_cresis_mat(filename)

        # check if Data is transposed
        # if Data.shape[0] == len(Time):
        #     pass
        # else:
        #     Data = np.transpose(Data)

        self.Frame      = frame
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
    # Method: load_segy
    #############################

    def load_segy(self, segy_file, Longitude='', Latitude='', dB=False):

        from lib.read_input import read_segy

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

        Data, Time, Frame, stream = read_segy(segy_file)

        self.Frame      = Frame
        self.Reader     = Reader
        self.Domain     = Domain
        self.Data       = Data
        self.Time       = Time
        self.Longitude  = Longitude
        self.Latitude   = Latitude
        self.dB         = dB
        self.Stream     = stream

        return self


    ########## END of load_segy() ###########


    #############################
    # Method: load_awi_nc
    #############################


    def load_awi_nc(self, nc_file, frame="", read_agc=False):

        '''


        '''

        from lib.read_input import read_awi_nc

        
        Data, Time, Longitude, Latitude, Aircraft_altitude, Ice_surface_elevation, Ice_thickness, Layer = read_awi_nc(nc_file, read_agc=read_agc)

        self.Frame      = frame
        self.Reader     = 'xarray'
        self.Domain     = 'twt'
        self.Data       = Data
        self.Time       = Time
        self.Longitude  = Longitude
        self.Latitude   = Latitude
        
        self.Ice_thickness = Ice_thickness

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


    #############################
    # Method: LAyERS
    #############################


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



    def check_surface_pick(self):

        import numpy as np

        try:
            missing = len(self.Longitude) - len(self.Layer['Surface']['trace'])
            if missing == 0:
                print("Surface pick OK.")
            elif missing > 0:
                for i in range(missing):
                    self.Layer['Surface']['trace'] = np.append(self.Layer['Surface']['trace'], self.Layer['Surface']['trace'][-1] + 1)
                    self.Layer['Surface']['value'] = np.append(self.Layer['Surface']['value'], self.Layer['Surface']['value'][-1])
                    print("Surface picks missing -- appending values")
                    print("Check suface pick...")
            elif missing < 0:
                self.Layer['Surface']['trace'] = self.Layer['Surface']['trace'][:missing]
                self.Layer['Surface']['value'] = self.Layer['Surface']['value'][:missing]
                print("Surface pick has more values than radar data traces -- clipping last values...")
                print("Check suface pick...")
        except:
            print("... check_surface_pick FAILED!")


    def get_layer_idx(self):

        import numpy as np

        layer_list = list(self.Layer.keys())
        time        = self.Time

        for lr in layer_list:
            values_idx = np.array([])
            v_twt      = self.Layer[lr]['value']

            for i in range(len(v_twt)):
                if v_twt[i] == np.nan:
                    v_idx = -9999
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




    ########################################
    # Method: correct_geometric_spreading()
    ########################################

    def correct_geometric_spreading(raw_object, v_ice=1.68914e8):

        from lib.radar_toolbox import correct_geometric_spreading
        import copy

        if raw_object.dB == False:
            pass
        else:
            print('==> Data should not be in dB on input.')

        try:
            new_obj = copy.deepcopy(raw_object)
        except:
            new_obj = copy.copy(raw_object)

        data     = new_obj.Data
        twt      = new_obj.Time
        surf_idx = new_obj.Layer['Surface']['value_idx']

        data_new = correct_geometric_spreading(data, twt, surf_idx, v_ice=v_ice)

        new_obj.Data = data_new
        new_obj.dB   = False

        return new_obj

    ########## END of correct_geometric_spreading() ###########



    ##########################################
    # Method: correct_englacial_attenuation()
    ##########################################

    def correct4constant_englacial_attenuation(raw_object, loss_factor=0, v_ice=1.68914e8):

        from lib.radar_toolbox import correct4constant_englacial_attenuation
        import copy

        if raw_object.dB == True:
            pass
        else:
            print('==> Data should be in dB on input.')

        loss_factor = loss_factor

        try:
            new_obj = copy.deepcopy(raw_object)
        except:
            new_obj = copy.copy(raw_object)

        data     = new_obj.Data
        twt      = new_obj.Time
        surf_idx = new_obj.Layer['Surface']['value_idx']

        data_new = correct4constant_englacial_attenuation(data, twt, surf_idx, v_ice=v_ice, loss_factor=loss_factor)

        new_obj.Data = data_new

        return new_obj


    ########## END of correct_englacial_attenuation() ###########



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

    def twt2surface(twt_object, padding=0):

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

        data_new, time_new =  radar_twt2surface(data=data, twt=twt, surf_idx=surf_idx, padding=padding)


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
        try:
            self.Ice_thickness    = self.Ice_thickness[::-1]
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
        try:
            layer_list = list(self.Layer.keys())
            for lr in layer_list:
                self.Layer[lr]["value"] = self.Layer[lr]["value"][start:end]
                self.Layer[lr]['trace'] = np.arange(int(end - start)).flatten() + 1
        except:
            pass
                
                
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

        try:
            self.Ice_thickness    = self.Ice_thickness[start:end]
        except:
            pass

        print('==> Clipped along-track: traces {}--{}'.format(start, end))


    ########## END of clip_along() ###########





    ##################################
    # Method: clip data in range
    ##################################

    def clip_range(self, start_val, end_val, mode='range_bin'):

        import numpy as np

        '''

        mode = 'range_bin'
        mode = 'elevation'
        mode = 'twt'
        '''

        if mode == 'range_bin':
            start = start_val
            end   = end_val

        elif mode == 'elevation':
            lower_lim = np.abs(self.Z - start_val).argmin()  # np.where(self.Z == start_val)[0][0]
            upper_lim = np.abs(self.Z - end_val).argmin()    # np.where(self.Z == end_val)[0][0]

            if lower_lim > upper_lim:
                start = upper_lim
                end   = lower_lim
            else:
                start = lower_lim
                end   = upper_lim

        elif mode == 'twt':
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
            # print(t_start, start, end, len(self.Time))
            
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
        Data          = []
        Longitude     = []
        Latitude      = []
        Elevation     = []
        GPS_time      = []
        Heading       = []
        Roll          = []
        Pitch         = []
        Spacing       = []
        Frames        = []
        Ice_thickness = []

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
            try:
                Ice_thickness.append(obj.Ice_thickness)
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
            new_obj.Ice_thickness  = np.concatenate(Ice_thickness)
        except:
            pass

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



    ########################
    # Method: relative power
    ########################

    def to_relative_power(self):

        import numpy as np

        self.Data = self.Data - np.nanmin(self.Data)


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

        from lib.radar_toolbox import automatic_gain_control2


        print('==> applying automatic gain control for layer sharpening')
        new_data  = automatic_gain_control2(self.Data, window=window)
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
        Segment = "dummy"                # out_object.Segment
        Season  = "dummy"                # out_object.Season

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

        # out_object.Data  = out_object.Data
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
