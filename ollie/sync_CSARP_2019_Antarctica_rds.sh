#
#
#
#
csarp_path=/work/ollie/sfranke/Scratch/rds/2019_Antarctica_Polar6/
hs_path=/hs/gsys/p_radar/antr2019/CSARP_2019_Antarctica_Polar6_rds/

## declare an array variable
declare -a folders=("CSARP_combined_3_sub_apertures" \
                    "CSARP_combined_5_sub_apertures" \
                    "CSARP_DEM" \
                    "CSARP_layer" \
                    "CSARP_layerData" \
                    "CSARP_music" \
                    "CSARP_noise" \
                    "CSARP_post_combined_5_sub_apertures" \
                    "CSARP_post_qlook" \
                    "CSARP_post_qlook_NMEA" \
                    "CSARP_post_standard" \
                    "CSARP_post_standard_sigma_x_1" \
                    "CSARP_layerData" \
                    "CSARP_qlook" \
                    "CSARP_qlook_NMEA" \
                    "CSARP_standard" \
                    "CSARP_standard_sigma_x_1" \
                    "CSARP_surfData" \
                    "CSARP_noise" )

## Sync Results
for folder in "${folders[@]}"
do
    echo ""
    echo "Syncing ===> $folder"
    rsync -avzhe ssh ${csarp_path}${folder} sfranke@hssrv1.awi.de:${hs_path}/
   # or do whatever with individual element of the array
done

## Sync metadata and support files

frames=/work/ollie/sfranke/Scratch/csarp_support/frames/rds/2019_Antarctica_Polar6
records=/work/ollie/sfranke/Scratch/csarp_support/records/rds/2019_Antarctica_Polar6
gps=/work/ollie/sfranke/Scratch/csarp_support/gps/2019_Antarctica_Polar6
metadata=/work/ollie/sfranke/Scratch/metadata/2019_Antarctica_Polar6

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












