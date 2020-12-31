import ee
from osgeo import ogr
import os
import geojson as geojson
import datetime
import numpy as np

date1 = datetime.datetime.strptime("20200601", '%Y%m%d')
date2 = datetime.datetime.strptime("20201101", '%Y%m%d')

def perdelta(start, end, delta):
    curr = start
    while curr < end:
        yield curr
        curr += delta

date_list = []
for result in perdelta(date1, date2, datetime.timedelta(days=1)):
	date_list = np.append(date_list,result)

ee.Authenticate()
ee.Initialize()
#S2 = ee.ImageCollection("COPERNICUS/S2_SR") ## Level-2 BOA
S2 = ee.ImageCollection("COPERNICUS/S2") ## Level-L1C TOA
gj = geojson.load(open("VG.geojson","r"))
coords = (gj["features"][0]["geometry"]["coordinates"][0])
polygon = ee.Geometry.Polygon(coords)

for i in np.arange(len(date_list)-1):
	print date_list[i].strftime("%Y-%m-%d"),date_list[i+1].strftime("%Y-%m-%d")
	S2_AOI = S2.filterDate(date_list[i].strftime("%Y-%m-%d"),date_list[i+1].strftime("%Y-%m-%d")).filterBounds(polygon)
	min = S2_AOI.min()
	#result = min.select('B4', 'B3', 'B2')
	result = min.select('B2')                                                                                                                                                               
	task_config = {'fileNamePrefix': 'S2'+'_'+date_list[i].strftime("%Y-%m-%d"),'crs': 'EPSG:32607','scale': 10,'fileFormat': 'GeoTIFF','skipEmptyTiles': True,'region': polygon, 'maxPixels': 300000000}
	task = ee.batch.Export.image.toDrive(result, 'S2'+'_'+date_list[i].strftime("%Y-%m-%d"), **task_config)
	task.start()


#date_list = [["2020-06-05","2020-06-06"],["2020-06-29","2020-06-30"],["2020-07-07","2020-07-08"],["2020-07-14","2020-07-15"],["2020-07-15","2020-07-16"],["2020-07-24","2020-07-25"],["2020-08-06","2020-08-07"],["2020-08-21","2020-08-22"],["2020-08-22","2020-08-23"],["2020-09-09","2020-09-10"]]

#for i in np.arange(len(date_list)):
#	print date_list[i][0],date_list[i][1]
#	S2_AOI = S2.filterDate(date_list[i][0],date_list[i][1]).filterBounds(polygon)
#	min = S2_AOI.min()
#	result = min.select('B4', 'B3', 'B2')                                                                                                                                                              
#	task_config = {'fileNamePrefix': 'S2'+'_'+date_list[i][0],'crs': 'EPSG:3413','scale': 30,'fileFormat': 'GeoTIFF','skipEmptyTiles': True,'region': polygon, 'maxPixels': 300000000}
#	task = ee.batch.Export.image.toDrive(result, 'S2'+'_'+date_list[i][0], **task_config)
#	task.start()
