#!/bin/bash

module load GMT

seed_points="../seed_points/EGRIP_CAMP_3413.shp"
vx="../velo/greenland_vel_mosaic250_vx_v1.tif"
vy="../velo/greenland_vel_mosaic250_vy_v1.tif "
out_file="../out_files/out.shp"

echo ""
echo "==> Starting generate_flowlines.py (incl. grd2stream)"
echo "==> seed_points = ${seed_points}"
echo "==> vx          = ${vx}"
echo "==> vy          = ${vy}"
echo "==> out_file    = ${out_file}"
echo ""

cd py_gmt_script

python generate_flowlines.py ${seed_points} ${vx} ${vy} ${out_file}

echo "Done..."
cd ..

