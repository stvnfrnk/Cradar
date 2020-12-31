#
#
#
#
csarp_path=/work/ollie/sfranke/Scratch/rds/2018_Greenland_Polar6/
hs_path=/hs/gsys/p_radar/arkr2018/CSARP_2018_Greenland_Polar6_rds/

## declare an array variable
declare -a folders=("CSARP_DEM" \
                    "CSARP_combined_3_sub_apertures" \
                    "CSARP_combined_5_sub_apertures" \
                    "CSARP_elevation_standard" \
                    "CSARP_layer" \
                    "CSARP_layerData" \
                    "CSARP_music" \
                    "CSARP_music_final" \
                    "CSARP_music_lr" \
                    "CSARP_music_lr_2" \
                    "CSARP_mvdr" \
                    "CSARP_noise" \
                    "CSARP_post_qlook" \
                    "CSARP_post_standard" \
		            "CSARP_post_sigma_x_1" \
		            "CSARP_qlook" \
                    "CSARP_single_channels" \
                    "CSARP_standard" \
                    "CSARP_sub_aperture_standard_1" \
                    "CSARP_sub_aperture_standard_2" \
                    "CSARP_sub_aperture_standard_3" \
                    "CSARP_sub_aperture_standard_4" \
                    "CSARP_sub_aperture_standard_5" \
                    "CSARP_standard_sigma_x_1" \
                    "CSARP_standard_sigma_x_1_sub_aperture_1" \
                    "CSARP_standard_sigma_x_1_sub_aperture_2" \
                    "CSARP_standard_sigma_x_1_sub_aperture_3" \
                    "CSARP_surfData")

## Sync Results
for folder in "${folders[@]}"
do
    echo ""
    echo "Syncing ===> $folder"
    rsync -avzhe ssh ${csarp_path}${folder} sfranke@hssrv1.awi.de:${hs_path}/
   # or do whatever with individual element of the array
done

## Sync metadata and support files

frames=/work/ollie/sfranke/Scratch/csarp_support/frames/rds/2018_Greenland_Polar6
records=/work/ollie/sfranke/Scratch/csarp_support/records/rds/2018_Greenland_Polar6
gps=/work/ollie/sfranke/Scratch/csarp_support/gps/2018_Greenland_Polar6
metadata=/work/ollie/sfranke/Scratch/metadata/2018_Greenland_Polar6

# csar_support files

echo ""
echo "Syncing ===> Frames"
ssh sfranke@hssrv1.awi.de mkdir -p ${hs_path}/csarp_support/frames/rds/
rsync -avzhe ssh ${frames} sfranke@hssrv1.awi.de:${hs_path}/csarp_support/frames/rds/

echo ""
echo "Syncing ===> Records"
ssh sfranke@hssrv1.awi.de mkdir -p ${hs_path}/csarp_support/records/rds/
rsync -avzhe ssh ${records} sfranke@hssrv1.awi.de:${hs_path}/csarp_support/records/rds/

echo ""
echo "Syncing ===> GPS"
ssh sfranke@hssrv1.awi.de mkdir -p ${hs_path}/csarp_support/gps/
rsync -avzhe ssh ${gps} sfranke@hssrv1.awi.de:${hs_path}/csarp_support/gps/

# metadata
echo ""
echo "Syncing ===> Metadata"
rsync -avzhe ssh ${metadata} sfranke@hssrv1.awi.de:${hs_path}/metadata/
#
#












