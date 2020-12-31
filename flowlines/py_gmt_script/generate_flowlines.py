# #!/usr/bin/env python2.7

import numpy as np
import os
from osgeo import ogr,osr
from shapely.wkt import loads
import sys

import faulthandler
faulthandler.enable()


def FindStringBetweenSign(file,sign):
	value = file.split(sign)
	return value

if (len(sys.argv) > 1):
	seed_file = sys.argv[1]
	vx = sys.argv[2]
	vy = sys.argv[3]
	outshape = sys.argv[4]

	seedname = FindStringBetweenSign(seed_file,"/")[-1][:-4]
	vxname = FindStringBetweenSign(vx,"/")[-1][:-4]
	vyname = FindStringBetweenSign(vy,"/")[-1][:-4]

	ds=ogr.Open(seed_file)
	lyr=ds.GetLayer()
	proj=lyr.GetSpatialRef()
	coord_list = []
	for feat in lyr:
		geom = feat.GetGeometryRef()
		x,y=geom.GetX(), geom.GetY()  #coord in map units
		coord_list = np.append(coord_list,x)
		coord_list = np.append(coord_list,y)
	#coord_list = coord_list.reshape(len(coord_list)/2,2)
	coord_list = coord_list.reshape(-1, 2)
	seeds = np.sort(coord_list.view('float,float'), order=['f0'], axis=0).view(np.float)
	np.savetxt(seedname+".txt",seeds,fmt='%s')

	#######generate flowlines
	if not os.path.exists(vxname+".grd"):
		os.system("gdalwarp "+vx+" -of netCDF "+vxname+".nc -tr 250 250 -r bilinear -overwrite")
		os.system("gdalwarp "+vy+" -of netCDF "+vyname+".nc -tr 250 250 -r bilinear -overwrite")
		os.system("gmt grdconvert "+vxname+".nc "+vxname+".grd")
		os.system("gmt grdconvert "+vyname+".nc "+vyname+".grd")
	os.system("./grd2stream "+vxname+".grd "+vyname+".grd -f "+seedname+".txt > grd2stream_out")
	
	flowlines = np.genfromtxt("grd2stream_out")
	break_list = np.where(np.isnan(flowlines[:,0]))[0]
	flowlines = np.delete(flowlines,break_list, axis=0)
	break_list_new = break_list-(np.arange(len(break_list)))
	break_list_new = np.delete(break_list_new, 0)
	flowlines_split = np.split(flowlines,break_list_new)
	profile_list = []
	new_i = 0
	for i in np.arange(len(flowlines_split)):
		if len(flowlines_split[i])>1:
			new_i = new_i+1
			line = ogr.Geometry(ogr.wkbLineString)
			tmplist = []
			for j in np.arange(len(flowlines_split[i][:,0])):
				lx = flowlines_split[i][j][0]
				ly = flowlines_split[i][j][1]
				tmplist.append(lx)
				line.AddPoint(lx,ly)
			if len(tmplist)>1:
				lineGeometry = ogr.CreateGeometryFromWkt(line.ExportToWkt())
				lineShapely = loads(line.ExportToWkt())
				profile_list.append(lineShapely)

	####create output shapefile
	dest_srs = proj
	outShapefile = outshape
	outDriver = ogr.GetDriverByName('Esri Shapefile')
	if os.path.exists(outShapefile):
		outDriver.DeleteDataSource(outShapefile)
	outDataSource = outDriver.CreateDataSource(outShapefile)
	layer = outDataSource.CreateLayer('', dest_srs, ogr.wkbLineString)
	layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
	defn = layer.GetLayerDefn()
	for i in np.arange(len(profile_list)):
		feat = ogr.Feature(defn)
		feat.SetField('id', str(i+1))
		geom = ogr.CreateGeometryFromWkb(profile_list[i].wkb)
		feat.SetGeometry(geom)
		layer.CreateFeature(feat)
		feat = geom = None
	outDataSource = layer = feat = geom = None
else:
	print('ok')
	#print("*** N.Neckel 2018 ***")
	#print("*** Generate flowlines via grd2stream from point shapefile and output line shapefile ***\n")
	#print("usage: generate_flowlines.py <input.shp> <vx.tif> <vy.tif> <output.shp>")
    #print("    input.shp         (input)  point shape file with seed points")
    #print("    vx.tif            (input)  vx GeoTIFF file")
	#print("    vy.tif            (input)  vy GeoTIFF file")
	#print("    output.shp        (output) line shape file")
