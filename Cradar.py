


class Cradar:
    
    #_allObjects = []
    
    # initiate an instance and read all important parameters 
    # from a CReSIS mcords matfile
    def __init__(self):
        #self._allObjects.append(self)
        pass
    
    def load(self, filename):
        
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

        # Delete the HDF5 file
        del self.File
        
        return self

    
    ########## END of load() ###########






    #############################
    # Method: calc_elevation 
    #############################
    
    def calc_elevation(twt_object, 
                       reference='', 
                       geotif='', 
                       region='', 
                       setting='', 
                       speed_of_ice=1.689e8, 
                       overlap=False, 
                       number_of_gaps=100):
        
        '''
        ==> Takes a Cradar object and transforms from twt domain 
            to the elvation domain.
        
        ==> The air-ice interface can either be set by the twt from the
            aircraft to the surface reflection or by a ice surface DEM
        
        reference       = 'reflection' or 'DEM'
        region          = Antarctica or Greenland, this is relevant for the reprojection
                          of the coordinates (lon, lat to X,Y). E.g. for the DEM grdtrack
                          ==> EPSG: 3413 for Greenland
                          ==> EPSG: 3031 for Antarctica
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
        #import pyproj
        #import geopy.distance 

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

            # load the function from geo_toolbox
            from geo_toolbox import extract_geotif_values

            # Check for the region
            if region == 'Greenland':
                EPSG = 3413
            if region == 'Antarctica':       
                EPSG = 3031

            # extract DEM values at shot locations
            elev_obj.DEM_surface = extract_geotif_values(geotif, elev_obj.Longitude, elev_obj.Latitude, EPSG)
            surface      = elev_obj.DEM_surface

        elif reference == 'GPS':
            surface = elev_obj.Surface

            # imput variables from instance
            data      = elev_obj.Data
            twt       = elev_obj.Time
            elevation = elev_obj.Elevation

            # input variables defined in prior steps
            surface   = surface

            # input variables from input of method above
            reference = reference 
            setting   = setting 
            overlap   = overlap 

        df = twt2elevation(data=data,
                           twt=twt,
                           elevation=elevation,
                           surface=surface,
                           speed_of_ice=1.689e8,
                           reference=reference,
                           setting=setting,
                           overlap=overlap,
                           overlap_traces=0
                           )
                
        

        # re-define instance atributes
        elev_obj.Data     = df
        elev_obj.Z        = df.index.values
        elev_obj.Domain   = 'Z'
        
        # define range resolution depending on the setting
        if setting == 'narrowband':
            elev_obj.Range_Resolution = '0.1 m'
            
        if setting == 'wideband':
            elev_obj.Range_Resolution = '1 m'
            
        if setting == 'snow':
            elev_obj.Range_Resolution = '0.001 m'
            
    
        return elev_obj
    
    
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
    
        self.Data      = self.Data.T[start:end].T
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
            self.Time[start:end]
    
        if domain == 'Z':
            self.Z[start:end]
        
        
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
    ################################## 
    
    def rename_frame(self, new_framename):
        self.Frame = new_framename
        
        
    ########## END of rename() ###########






    #############################
    # Method: write matfile
    #############################
    
    def write_shape(self, out_filename='', out_format='shapefile'):
        
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
            shape_filename = out_object.Frame + '_' + out_object.Domain

        if out_format == '':
            out_format = 'shapefile'


        out = coords2shape(X, Y, EPSG_in=4326, EPSG_out=4326, geometry='Point', attributes='')

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
        
        
    ########## END of write_mat() ###########
    
    
    
    
    
    
    #####################################
    # Converts CReSIS format .mat files
    # to SEGY format
    #####################################

    def write_segy(self, region='', out_filename='', differenciate=False, step=1):

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
        #import time
        #import sys

        print('')
        print('==> Processing Frame: {} located in {}'.format(self.Frame, region))
        print('==> This file is in >> {} << domain'.format(self.Domain))
        
        # Check on TWT or Elevation Data
        if self.Domain == 'twt':
            domain              = self.Time
            receiver_elevation  = 0
            num_of_samples      = len(self.Time)
            segy_filename       = self.Frame + '_twt.segy'
        elif self.Domain == 'Z':
            domain              = self.Z
            receiver_elevation  = self.Z.max()
            num_of_samples      = len(self.Time)
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
        X, Y        = transformer.transform(Lon, Lat)

        # Figure out TWT sampling interval
        diff_domain = np.array([])

        for i in range(len(domain) - 1):
            diff_domain = np.append(diff_domain, domain[i] - domain[i+1])

        # get sample interval | doesn't matter if twt or Z
        sample_interval = diff_domain.mean()
        
        # apply radar2segy method
        stream = radar2segy(data=data, 
                            receiver_elevation=receiver_elevation, 
                            num_of_samples=num_of_samples, 
                            sample_interval=sample_interval,
                            X=X, 
                            Y=X,
                            step=1,
                            time_mode='gmtime',
                            gps_time=gps_time,
                            differenciate=differenciate
                            )

        if out_filename=='':
            pass
        else: 
            segy_filename = out_filename


        stream.write(segy_filename, format='SEGY', data_encoding=5, byteorder='>',textual_file_encoding='ASCII')
        print('==> Written: {}'.format(segy_filename))


        