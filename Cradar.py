


class Cradar:
    
    #_allObjects = []
    
    # initiate an instance and read all important parameters 
    # from a CReSIS mcords matfile
    def __init__(self):
        #self._allObjects.append(self)
        pass
    
    def load(self, filename, dB=False):
        
        import h5py
        import scipy.io
        import numpy as np
        import pandas as pd
        
        #self._allObjects.append(self)


        ##############################
        # Try with scipy.io.loadmat()
        try:
            self.File      = scipy.io.loadmat(filename)
            self.Frame     = filename.split('.mat')[0]
            self.Reader    = 'scipy'
            
            # Iterate over almost all items in HDF5 File
            for k, v in self.File.items():
                if 'Time' not in k:
                    try:
                        setattr(self, k, np.array(v)[0])
                    except:
                        setattr(self, k, v)

                if 'Time' in k:
                    self.Time = np.array(self.File['Time'])[0]

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
                    try:
                        setattr(self, k, np.array(v).T[0])
                    except:
                        setattr(self, k, v)
                    
                if 'Time' in k:
                    self.Time = np.array(self.File['Time'])[0]

                if 'Time' in k and 'Z' not in k:
                    self.Domain    = 'twt'

            self.Data = pd.DataFrame(np.array(self.File['Data'])).T

            if dB == False:
                self.dB = False
            elif dB == True:
                self.dB = True
            else:
                print('==> dB True or False? Is set to False.')

        # Delete the HDF5 file
        del self.File
        
        return self

    
    ########## END of load() ###########




    #############################
    # Method: add_raster_values
    #############################

    def gridtrack(self, geotif='', geotif_name=''):

        import pygmt
        import rioxarray
        import pandas as pd

        geotif      = geotif
        geotif_name = geotif_name
        Longitude   = self.Longitude
        Latitude    = self.Latitude

        df             = pd.DataFrame(Longitude)
        df['Latitude'] = pd.DataFrame(Latitude)
        df.columns     = ['Longitude', 'Latitude']

        rds = rioxarray.open_rasterio(geotif)
        rds = rds.rio.reproject("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
        rds = rds.squeeze('band')
        rds = rds.astype(float)
        df  = pygmt.grdtrack(df, rds, newcolname=geotif_name)
        
        raster_vals = df[geotif_name].values

        setattr(self, geotif_name, raster_vals)
        print('==> Added {} to the data'.format(geotif_name))

        del df
        del rds
        del raster_vals


    ########## END of rename() ###########



    #############################
    # Method: get_surf_idx()
    #############################

    def get_surf_idx(self):

        import numpy as np

        twt      = self.Time
        twt_surf = self.Surface
        traces   = self.Data.shape[1]

        surf_idx = np.array([])

        for i in range(traces):
            s_idx = (np.abs(np.array(twt) - np.array(twt_surf)[i])).argmin()
            surf_idx = np.append(surf_idx, s_idx)
            #print(surf_idx)

        self.Surface_idx = surf_idx

        del twt, twt_surf, traces, surf_idx


    ########## END of get_surf_idx() ###########




    #############################
    # Method: correct_geom_spreading 
    #############################

    def correct4attenuation(raw_object, mode=0, loss_factor=0):

        from radar_toolbox import correct4attenuation
        import copy

        mode        = mode
        loss_factor = loss_factor

        geom_obj = copy.deepcopy(raw_object)

        geom_obj.get_surf_idx()

        data     = geom_obj.Data
        twt      = geom_obj.Time
        surf_idx = geom_obj.Surface_idx

        data_new = correct4attenuation(data, twt, surf_idx, v_ice=1.68914e8, mode=mode, loss_factor=loss_factor)

        geom_obj.Data = data_new

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
        
        reference       = 'reflection' or 'DEM'
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
        elev_obj = copy.deepcopy(twt_object)
        
        
        print('==> twt2elevation...')
        
        if reference == 'GPS':
            print('==> Using Aircraft GPS and radar surf. reflection to derive elevation')
            
        elif reference == 'DEM':
            print('==> Using a Ice surf. DEM to derive elevation')
        
        elif reference == '':
            reference = 'GPS'
            print("==> !! You didn't define a reference, it is now automatically set to 'GPS'")
        
        
        # Get ice surface elevation values from DEM
        if reference == 'DEM': 
            twt_surface = elev_obj.Surface
            DEM_surface = elev_obj.DEM_surface

        elif reference == 'GPS':
            twt_surface = elev_obj.Surface
            DEM_surface = ''
        
        # imput variables from instance
        data               = elev_obj.Data
        twt                = elev_obj.Time
        aircraft_elevation = elev_obj.Elevation

        # input variables defined in prior steps
        twt_surface   = twt_surface

        # input variables from input of method above
        reference = reference 
        setting   = setting 
        overlap   = overlap 

        df = twt2elevation(data=data,
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
        elev_obj.Z        = df.index.values
        elev_obj.Data     = pd.DataFrame(np.array(df))
        elev_obj.Domain   = 'Z'
        
        # define range resolution depending on the setting
        if setting == 'narrowband':
            elev_obj.Range_Resolution = '0.1 m'
            
        if setting == 'wideband':
            elev_obj.Range_Resolution = '1 m'
            
        if setting == 'snow':
            elev_obj.Range_Resolution = '0.001 m'
            
    
        return elev_obj

        del df
    
    
    ########## END of twt2elevation() ###########
    
    
        
        
        
        
    #############################
    # Method: write matfile
    #############################
    
    def flip_lr(self):
        
        '''
        Flips the whole radar matrix and all its along-track attributes
        '''

        self.Data      = self.Data.T[::-1].T
        self.Longitude = self.Longitude[::-1]
        self.Latitude  = self.Latitude[::-1]
        self.GPS_time  = self.GPS_time[::-1]
        self.Surface   = self.Surface[::-1]
        self.Elevation = self.Elevation[::-1]
        
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
        self.GPS_time  = self.GPS_time[start:end]
        self.Surface   = self.Surface[start:end]
        self.Elevation = self.Elevation[start:end]
        
        # optional
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
        

    ########## END of clip_along() ###########

    
    
    
    
    ##################################
    # Method: clip data in range
    ##################################  
        
    def clip_range(self, start, end):

        '''
        '''
        
        self.Data      = self.Data[start:end]

        domain     = self.Domain
        
        if domain == 'twt':
            self.Time = self.Time[start:end]
            # if not starting at 0, ggf. noch was tun?
    
        if domain == 'Z':
            self.Z = self.Z[start:end]
        
        
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
        new_obj = copy.deepcopy(added_objects[0])
        
        Data      = []
        Longitude = []
        Latitude  = []
        Elevation = []
        GPS_time  = []
        Surface   = []
        Bottom    = []
        Heading   = []
        Roll      = []
        Pitch     = []
        Spacing   = []
        Frames    = []
        
        
        for obj in added_objects:
            
            if obj.Domain == 'Z':
                obj.Data.index = obj.Z
            elif obj.Domain == 'twt':
                obj.Data.index = obj.Time
            
            Data.append(obj.Data)
            Longitude.append(obj.Longitude)
            Latitude.append(obj.Latitude)
            Elevation.append(obj.Elevation)
            GPS_time.append(obj.GPS_time)
            Surface.append(obj.Surface)
            Bottom.append(obj.Bottom)
            Heading.append(obj.Heading)
            Roll.append(obj.Roll)
            Pitch.append(obj.Pitch)
            Frames.append(obj.Frame)
            
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
                
        new_obj.Longitude = np.concatenate(Longitude)
        new_obj.Latitude  = np.concatenate(Latitude)
        new_obj.Elevation = np.concatenate(Elevation)
        new_obj.GPS_time  = np.concatenate(GPS_time)
        new_obj.Surface   = np.concatenate(Surface)
        new_obj.Bottom    = np.concatenate(Bottom)
        new_obj.Heading   = np.concatenate(Heading)
        new_obj.Roll      = np.concatenate(Roll)
        new_obj.Pitch     = np.concatenate(Pitch)
        new_obj.Frames    = Frames
        new_obj.Frame     = added_objects[0].Frame + '_concat'
        
        return new_obj
    
    
    ########## END of concat_frames() ###########
    
    
    
    
    
    ##################
    # Method: rename
    ##################
    
    def rename_frame(self, new_framename):
        self.Frame = new_framename
        
        
    ########## END of rename() ###########



    



    #############################
    # Method: write shape
    #############################
    
    def write_shape(self, out_filename='', out_folder='', out_format='shapefile'):
        
        import numpy as np
        import geopandas
        from geo_toolbox import coords2shape
        import copy
        import os

        out_object = copy.deepcopy(self)

        X = out_object.Longitude
        Y = out_object.Latitude

        out_filename = out_filename
        out_format   = out_format

        if out_filename == '':
            shape_filename = out_object.Frame# + '_' + out_object.Domain

        if out_format == '':
            out_format = 'shapefile'


        out = coords2shape(X, Y, EPSG_in=4326, EPSG_out=4326, geometry='Point', attributes='')

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

        else:
            if out_format == 'shapefile':
                out.to_file(out_folder + '/' + shape_filename + '.shp')
                print('==> Written: {}/{}.shp'.format(out_folder, shape_filename))

            if out_format == 'geojson':
                out.to_file(out_folder + '/' + shape_filename + '.geojson', driver='GeoJSON')
                print('==> Written: {}/{}.geojson'.format(out_folder, shape_filename))
        
        del X, Y, out
        
    ########## END of write_mat() ###########





    #############################
    # Method: write matfile
    #############################
    
    def write_mat(self, out_filename=''):
        
        import scipy.io
        import pandas as pd
        import numpy as np
        import copy

        out_object = copy.deepcopy(self)
                        
        out_object.Data    = out_object.Data.values
        full_dict          = out_object.__dict__
        mat_filename       = out_object.Frame + '_' + out_object.Domain + '.mat'

        if out_filename=='':
            pass
        else: 
            mat_filename = out_filename
        
        scipy.io.savemat(mat_filename, full_dict)
        print('==> Written: {}'.format(mat_filename))

        del out_object
        del full_dict
        del mat_filename
        
        
    ########## END of write_mat() ###########
    
    
    
    
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

        print('')
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
        gps_time = self.GPS_time
        
        # re-project lon, lat to X, Y depending on EPSG
        if region == 'Greenland':
            EPSG = 3413
        elif region == 'Antarctica':
            EPSG = 3031
        else:
            print('define region....')

        transformer = Transformer.from_crs(4326, EPSG)
        Lon, Lat    = self.Longitude, self.Latitude
        # ==> BEWARE I think here is a bug, Lon and Lat switched in the function!!
        X, Y        = transformer.transform(Lat, Lon) 



        # Figure out TWT sampling interval
        diff_domain = np.array([])

        for i in range(len(domain) - 1):
            diff_domain = np.append(diff_domain, domain[i] - domain[i+1])

        # get sample interval | doesn't matter if twt or Z
        sample_interval = np.abs(diff_domain.mean()) * 1000000
        
        # apply radar2segy method
        stream = radar2segy(data=data, 
                            receiver_elevation=receiver_elevation, 
                            num_of_samples=num_of_samples, 
                            sample_interval=sample_interval,
                            X=X, 
                            Y=Y,
                            step=1,
                            time_mode='gmtime',
                            gps_time=gps_time,
                            differenciate=differenciate,
                            to_dB=to_dB
                            )

        if out_filename=='':
            pass
        else: 
            segy_filename = out_filename

        self.Stream = stream

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
    # Converts CReSIS format .mat files
    # to SEGY format
    #####################################

    def plot_overview(self, flight_lines, save_png=True, dpi=100, out_folder='', cmap='binary', divergent=False, show=True):

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

        yticks  = np.array(range(0, len(self.Time)))[0::152]
        ylabels = (self.Time[0::152]*1000000).astype(int)


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
        plt.ylabel('TWT (µs)', fontsize=16)
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

                plt.savefig(out_folder + '/' + figname, dpi=dpi, bbox_inches='tight')
                print('==> Written: {}/{}'.format(out_folder, figname))

        if show == True:
            plt.show()
        else:
            plt.clf()
            plt.close('all')

        del xticks, xlabels, yticks, ylabels
        del frame, first, survey_lines
        del df, flight_lines, Lon, Lat






