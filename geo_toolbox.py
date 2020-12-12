
'''
 1. Function: extract_geotif_values
 2. Function: radartrack2shape

'''

##############################################################################

#################################
# Function extract_geotif_values
#################################

def extract_geotif_values(geotif, data_frame, EPSG=''):
    
    '''
     
    Required Libraries:
            
        - numpy
        - pandas
        - osgeo / gdal
        - pyproj
        
    
    Import as:
        
        import extract_geotif_values
        
        Usage: df['New Column'] = eextract_geotif_values(geotif, data_frame, EPSG)
        
    
    Valid geotif file and pandas DataFrame with either X and Y as EPSG3413 coordinates
    or lon, lat (lowercase)
    
    '''


    import numpy as np
    import pandas as pd
    from osgeo import gdal
    import pyproj
    
    
    df = data_frame
    
    
    if 'X' and 'Y' in df.columns:
        pass
    
    else:
        if EPSG == 3413:
            EPSG3413=pyproj.Proj("EPSG:3413")
            df['X'], df['Y'] = EPSG3413(np.array(df['Longitude']), np.array(df['Latitude']))
        
        elif EPSG == 3031:
            EPSG3031=pyproj.Proj("EPSG:3031")
            df['X'], df['Y'] = EPSG3031(np.array(df['Longitude']), np.array(df['Latitude']))

        else:
            print('No valid X and Y files or EPSG specified:')
            pass
    
    
    #driver = gdal.GetDriverByName('GTiff')
    filename = geotif
    dataset = gdal.Open(filename)
    #band = dataset.GetRasterBand(1)
    
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    
    transform = dataset.GetGeoTransform()
    
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]
    
    #data = band.ReadAsArray(0, 0, cols, rows)
    data = dataset.ReadAsArray(0, 0, cols, rows)
    
    position = df[['X', 'Y']]
    points_list = [tuple(x) for x in position.values] #list of X,Y coordinates
    
    lst = []
    
    for point in points_list:
        try:
            col = int((point[0] - xOrigin) / pixelWidth)
            row = int((yOrigin - point[1] ) / pixelHeight)

            try:
                lst.append(data[row][col])
            except:
                print('problem at point: {} ==> appending last value: {}'.format(point, lst[-1]))
                try:
                    lst.append(lst[-1])
                except:
                    lst.append(np.nan)
                    
        except: #OverflowError:
            lst.append(np.nan)
        
    out_column = pd.DataFrame(lst)
    
    return out_column













#################################
# Function extract_geotif_values
#################################

def gridtrack(geotif, X, Y, EPSG_xy=0, EPSG_raster=0):
    
    '''
    Required Libraries:
            
        - numpy
        - pandas
        - osgeo / gdal
        - pyproj
        
    
    Import as:
        
        import extract_geotif_values
        
        Usage: df['New Column'] = eextract_geotif_values(geotif, data_frame, EPSG)
        
    
    Valid geotif file and pandas DataFrame with either X and Y as EPSG3413 coordinates
    or lon, lat (lowercase)
    
    '''


    import numpy as np
    import pandas as pd
    from osgeo import gdal
    import pyproj
    
    geotif      = geotif
    X_in        = X
    Y_in        = Y
    EPSG_xy     = EPSG_xy
    EPSG_raster = EPSG_raster

    # convert coordinates from EPSG_in to EPSG_out
    transformer = Transformer.from_crs(EPSG_xy, EPSG_raster)
    X, Y        = transformer.transform(X_in, Y_in)
    
    
    #driver = gdal.GetDriverByName('GTiff')
    filename = geotif
    dataset = gdal.Open(filename)
    #band = dataset.GetRasterBand(1)
    
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    
    transform = dataset.GetGeoTransform()
    
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]
    
    #data = band.ReadAsArray(0, 0, cols, rows)
    data = dataset.ReadAsArray(0, 0, cols, rows)
    
    #position = df[['X', 'Y']]
    points_list = list(tuple(zip(X,Y))) #list of X,Y coordinates
    
    lst = []
    
    for point in points_list:
        try:
            col = int((point[0] - xOrigin) / pixelWidth)
            row = int((yOrigin - point[1] ) / pixelHeight)

            try:
                lst.append(data[row][col])
            except:
                print('problem at point: {} ==> appending last value: {}'.format(point, lst[-1]))
                try:
                    lst.append(lst[-1])
                except:
                    lst.append(np.nan)
                    
        except: #OverflowError:
            lst.append(np.nan)
        
    out_column = np.array(lst)
    
    return out_column



#############################################
# Function: coords2shape
# creates shape files (point, or linestring)
# from a cresis radar .mat file
#############################################

def coords2shape(X, Y, EPSG_in=4326, EPSG_out=4326, geometry='Point', attributes=[]):

    '''
        Usage ==>> coords2shape(x, y, EPSG, attributes=[])
        geometry types:
                        - Point
                        - LineString
                        - Both
                        
        EPSG (CRS):     Set the output EPSG.
                        - 4326 (standard, uses Lat/Lon)
                        - 3413 (for Greenland  - NSIDC Sea Ice Polar Stereographic North)
                        - 3031 (for Antarctica - Antarctic Polar Stereographic)
                        
        Attributes:     List of Attributes to hand over to shape files
                        - Example:                    
                        - ['X', 'Y', 'GPS_Time', 'Aircraft_Elevation'] 
    '''


    import numpy as np
    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import Point, LineString
    from pyproj import Transformer

    X          = X
    Y          = Y
    EPSG_in    = EPSG_in
    EPSG_out   = EPSG_out
    attributes = attributes


    if EPSG_in and EPSG_out == 4326:
        pass
    else:
        X_in = X
        Y_in = Y

        # convert coordinates from EPSG_in to EPSG_out
        transformer = Transformer.from_crs(EPSG_in, EPSG_out)
        X, Y        = transformer.transform(X_in, Y_in)

    # create the data frame with the coordinates as well as the attributes
    df         = pd.DataFrame(X)
    df['Y']    = pd.DataFrame(Y)
    df.columns = ['X', 'Y']

    if attributes:
        for key, value in attributes.items():
            df[key] = value 
    else:
        pass

    # create geopandas data frame with geometry
    gdf = gpd.GeoDataFrame(df, crs=EPSG_out, geometry=gpd.points_from_xy(df['X'], df['Y']))

    return gdf



























###########################################################################
###########################################################################
###########################################################################
#
#                   OLD
#
###########################################################################
###########################################################################
###########################################################################

    
def radartrack2shape(file, geometry='Point', EPSG=4326, radar='uwb', Attributes=[]):
    
    '''
        Usage ==>> radartrack2shape(file, geometry='Point' or 'LineString', 
                                    EPSG=4326/3413/3031, 
                                    Attributes=[] or List with Attributes in matfile):
        
        geometry types:
                        - Point
                        - LineString
                        - Both
                        
        EPSG (CRS):
                        - 4326 (standard, uses Lat/Lon)
                        - 3413 (for Greenland  - NSIDC Sea Ice Polar Stereographic North)
                        - 3031 (for Antarctica - Antarctic Polar Stereographic)
                        
        Attributes:     List of Attributes to hand over to shape files
        
                            - Example:                    
                            - ['X', 'Y', 'GPS_Time', 'Aircraft_Elevation'] 
    '''
        
    
    # Importing Libraries
    import scipy.io
    import h5py
    import os
    import numpy as np
    import pandas as pd
    from collections import OrderedDict
    from shapely.geometry import Point, LineString, mapping
    from fiona import collection
    from fiona.crs import from_epsg
    

    if radar == 'uwb':
        print('===> Dealing with UWB radar data...')

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
            
            # Lat and Lon Required
            Latitude  = np.array(mat['Latitude'])[0]
            Longitude = np.array(mat['Longitude'])[0]


        #loading with h5py
        if reader == 'h5py':
            print('Using h5py to load matfile')

            # Lat and Lon Required
            Latitude  = np.array(mat['Latitude']).T[0]
            Longitude = np.array(mat['Longitude']).T[0]


    elif radar == 'emr':
        print('===> Dealing with EMR radar data...')

        koord       = pd.read_csv(file, delim_whitespace=True)
        Longitude   = np.array(koord['GPSLon'])
        Latitude    = np.array(koord['GPSLat'])


    else:
        print("Pleas provide a radar type ('uwb' or 'emr')")
    
    shape_name  = file.split('.')[0]
    CRS         = from_epsg(EPSG)

    
    # Define X and Y according to EPSG
    if EPSG == 4326:
        
        X_  = 'Longitude'
        Y_  = 'Latitude'
        X   = Longitude
        Y   = Latitude
        
    # Greenland - NSIDC Sea Ice Polar Stereographic North  
    if EPSG == 3413:
        
        import pyproj
        
        Projection  = "EPSG:3413"
        EPSG_       = pyproj.Proj(Projection)
        X, Y        = EPSG_(Longitude, Latitude)
        X_, Y_      = 'X', 'Y'
        
    # Antarctica - Antarctic Polar Stereographic 
    if EPSG == 3031:
        
        import pyproj
        
        Projection  = "EPSG:3031"
        EPSG_       = pyproj.Proj(Projection)
        X, Y        = EPSG_(Longitude, Latitude)
        X_, Y_      = 'X', 'Y'
        
        
    ###################
    # Create Folder
    ###################

    if not os.path.exists('shapes'):
        os.makedirs('shapes')


    ###################
    # Shapefile POINT
    ###################
    
    if geometry == 'Point':
        
        shape_name = shape_name + '_points.shp'
        
        schema = {'geometry'    : geometry, 
                  'properties'  : OrderedDict([
                                             ('file_name', 'str'),
                                             (X_, 'float'),
                                             (Y_, 'float'),
                                             ('Trace', 'int')
                                             ])
                                             }
        
        with collection('shapes/' + shape_name, mode='w', driver='ESRI Shapefile',crs=CRS, schema=schema) as output:
            for i in range(len(Longitude) - 1):
                
                point = Point(Longitude[i], Latitude[i])
                output.write({
                    'properties': OrderedDict([
                                             ('file_name', file),
                                             (X_, X[i]),
                                             (Y_, Y[i]),
                                             ('Trace', i)
                                             ]),
                    'geometry': mapping(point)})
    
            print('')
            print('===>> Converted {} to Point Shapefile'.format(shape_name))
    
    
    #################
    # Shapefile LINE
    #################
    
    if geometry == 'LineString':
        
        shape_name = shape_name + '_line.shp'
        
        schema = {'geometry'    : geometry, 
                  'properties'  : OrderedDict([
                                             ('file_name', 'str')
                                             ])
                                             }
        
        line = LineString(list(zip(X, Y)))
        
        with collection('shapes/' + shape_name, mode='w', driver='ESRI Shapefile',crs=CRS, schema=schema) as output:
            output.write({
                'geometry': mapping(line),
                'properties': {'file_name': file},
            })
            
        print('')
        print('===>> Converted {} to LineString Shapefile'.format(shape_name))


