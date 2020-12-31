#!/bin/bash
season=2018_Greenland_Polar6
#season=2019_Antarctica_Polar6
path=/work/ollie/sfranke/Scratch/rds/${season}
segment=$1
num_apertures=3
CSARP_folder=CSARP_sar_sigma_x_1_sub_aperture

# calculate frames of segment
num_frames=`find ${path}/${CSARP_folder}/${segment}/ -type d | grep _01_01 | wc -l`

# 1. create target folders
echo ""
echo "1. creating directories"

for frame in $(seq -f "%02g" 1 ${num_frames}); do
	for i in {1..3}; do
	mkdir -p ${path}/${CSARP_folder}_${i}/${segment}/fk_data_0${frame}_01_01/; done; done
echo "--> done."
echo ""

# 2. copy sar_coord.mat
echo "2. copying sar_coord.mat"
for i in {1..3}; do
	cp ${path}/${CSARP_folder}/${segment}/sar_coord.mat ${path}/${CSARP_folder}_${i}/${segment}/; done
echo "--> done."
echo ""

# 3. copy fk_data into corresponding sub aperture folder
for frame in $(seq -f "%02g" 1 ${num_frames}); do
	for i in {1..3}; do 
		echo "3. moving fk_data: ${CSARP_folder}/${segment}/fk_data_0${frame}_0${i}_01 --> ${CSARP_folder}_${i}/${segment}/fk_data_0${frame}_01_01"
		mv ${path}/${CSARP_folder}/${segment}/fk_data_0${frame}_0${i}_01/* ${path}/${CSARP_folder}_${i}/${segment}/fk_data_0${frame}_01_01/ 
	done 
done

echo ""
echo "--> done with " $segment


