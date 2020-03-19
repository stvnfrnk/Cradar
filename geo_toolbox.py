
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
    
    
    df= data_frame
    
    
    if 'X' and 'Y' in df.columns:
        pass
    
    else:
        if EPSG == '3413':
            EPSG3413=pyproj.Proj("+init=EPSG:3413")
            df['X'], df['Y'] = EPSG3413(np.array(df['Longitude']), np.array(df['Latitude']))
        
        if EPSG == '3031':
            EPSG3031=pyproj.Proj("+init=EPSG:3031")
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
        col = int((point[0] - xOrigin) / pixelWidth)
        row = int((yOrigin - point[1] ) / pixelHeight)
    
        #print(row,col, data[row][col])
        lst.append(data[row][col])
        
    out_column = pd.DataFrame(lst)
    
    return out_column



#%%
#############################
# Function: radartrack2shape
#############################
    
def radartrack2shape(file, geometry='Point', EPSG=4326, Attributes=[]):
    
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
    from collections import OrderedDict
    from shapely.geometry import Point, LineString, mapping
    from fiona import collection
    from fiona.crs import from_epsg
    


    # depending on the .mat file version
    # either scipy.io (older versions) or h5py (newer versions)
    # will be used to load the file        
    try:
        mat     = scipy.io.loadmat(file)
        reader  = 'scipy'
    except NotImplementedError:
        mat     = h5py.File(file)
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
        
        Projection  = "+init=EPSG:3413"
        EPSG_       = pyproj.Proj(Projection)
        X, Y        = EPSG_(Longitude, Latitude)
        X_, Y_      = 'X', 'Y'
        
    # Antarctica - Antarctic Polar Stereographic 
    if EPSG == 3031:
        
        import pyproj
        
        Projection  = "+init=EPSG:3031"
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


