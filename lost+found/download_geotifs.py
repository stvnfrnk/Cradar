from osgeo import ogr, gdal
import os
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
import datetime
import numpy as np
import fnmatch
import datetime
import netCDF4
import re
from import_export_geotiff import resampleGeoTiff
from pathlib import Path

def FindStringBetweenSign(file,sign):
	value = file.split(sign)
	return value

#####Download and delete files on Google Drive
#gauth = GoogleAuth()
#gauth.LocalWebserverAuth() # client_secrets.json need to be in the same directory as the script
#drive = GoogleDrive(gauth)
#fileList = drive.ListFile({'q': "'root' in parents and trashed=false"}).GetList()
#for file in fileList:
#	print('Title: %s, ID: %s' % (file['title'], file['id']))
#	if file['title'].startswith("S2"):
#		file.GetContentFile(file['title'])
#		file.Delete()

date_list = []
for file in os.listdir("./"):
	if fnmatch.fnmatch(file, 'S2**.tif'):
		os.system("grdcompress.sh "+file)
		date = datetime.datetime.strptime(FindStringBetweenSign(file,"_")[1][:10], '%Y-%m-%d')
		print date
		date_list = np.append(date_list,date)
		date_list = np.append(date_list,file)
		date_list = np.append(date_list,Path(file).stat().st_size)
date_list = date_list.reshape(int(len(date_list)/3),3)
date_list = date_list[np.argsort(date_list[:,0])]

master = date_list[0,1]
ds = gdal.Open(master)
a = ds.ReadAsArray()
ny,nx = np.shape(a)

b = ds.GetGeoTransform() #bbox, interval
x = np.arange(nx)*b[1]+b[0]
y = np.arange(ny)*b[5]+b[3]

basedate = datetime.datetime(1858,11,17,0,0,0)

# create NetCDF file
nco = netCDF4.Dataset('S2_time_series.nc','w',clobber=True)

# chunking is optional, but can improve access a lot: 
# (see: http://www.unidata.ucar.edu/blogs/developer/entry/chunking_data_choosing_shapes)
chunk_x=16
chunk_y=16
chunk_time=12

# create dimensions, variables and attributes:
nco.createDimension('x',nx)
nco.createDimension('y',ny)
nco.createDimension('time',None)
timeo = nco.createVariable('time','f4',('time'))
timeo.units = 'days since 1858-11-17 00:00:00'
timeo.standard_name = 'time'

xo = nco.createVariable('x','f4',('x'))
xo.units = 'm'
xo.standard_name = 'projection_x_coordinate'

yo = nco.createVariable('y','f4',('y'))
yo.units = 'm'
yo.standard_name = 'projection_y_coordinate'

# create container variable for CRS: x/y WGS84 datum
crso = nco.createVariable('crs','i4')
crso.grid_mapping_name='polar_stereographic'
crso.straight_vertical_longitude_from_pole = 0.
crso.latitude_of_projection_origin = 90.
crso.scale_factor_at_projection_origin = 1.0
crso.false_easting = 0.0
crso.false_northing = 0.0
crso.semi_major_axis = 6378137.0
crso.inverse_flattening = 298.257223563
crso.standard_parallel = 70.
crso.longitude_of_prime_meridian = -45.

# create short integer variable for temperature data, with chunking
tmno = nco.createVariable('z', 'f4',  ('time', 'y', 'x'), 
   zlib=True,chunksizes=[chunk_time,chunk_y,chunk_x],fill_value=0.)
tmno.units = 'm'
tmno.scale_factor = 1
tmno.add_offset = 0.00
tmno.long_name = 'surface elevation'
tmno.standard_name = 'r'
tmno.grid_mapping = 'crs'
tmno.set_auto_maskandscale(False)

nco.Conventions='CF-1.6'

#write x,y
xo[:]=x
yo[:]=y

pat = re.compile('us_tmin_[0-9]{4}\.[0-9]{2}')
itime=0

dtime=(date_list[0,0]-basedate).total_seconds()/86400.
for l in np.arange(len(date_list)):
	dataref,dataslave,width,height,match_geotrans,match_proj = resampleGeoTiff(master,date_list[l,1])
	dataref = np.where(dataref==-9999,np.nan,dataref)
	dataslave = np.where(dataslave==0,np.nan,dataslave)
	if not np.isnan(np.nanmean(dataref-dataslave)):
		dtime=(date_list[l,0]-basedate).total_seconds()/86400.
		timeo[itime]=dtime
		tmno[itime,:,:]=dataslave
		itime=itime+1
nco.close()
##os.system("rm *.tif")

