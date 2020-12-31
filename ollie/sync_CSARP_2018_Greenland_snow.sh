#
#
#
#
csarp_path=/work/ollie/rzindler/Scratch/snow/2018_Greenland_Polar6/
#hs_path=/hs/gsys/geophy/projects/arkr2018/CSARP_2018_Greenland_Polar6_rds/
hs_path=/hs/gsys/p_radar/arkr2018/CSARP_2018_Greenland_Polar6_snow_rzindler/

## declare an array variable
declare -a folders=("CSARP_analysis" \
		    "CSARP_analysis_noise_old" \
		    "CSARP_layerData" \
		    "CSARP_layerData_pick" \
                    "CSARP_mvdr" \
                    "CSARP_post" \
		    "CSARP_post_Abbildungen/" \
		    "CSARP_post_deconv" \
		    "CSARP_post_deconv_changedWF" \
                    "CSARP_post_deconv_updated" \
                    "CSARP_post_noise_changedWF" \
		    "CSARP_post_noise_dielectric_test" \
		    "CSARP_post_noise_imgcomb_test" \
		    "CSARP_post_qlook" \
		    "CSARP_post_sar" \
		    "CSARP_post_SAR" \
		    "CSARP_post_sar_noise_changedWF" \
		    "CSARP_post_sar_updatedWF" \
		    "CSARP_post_standard" \
                    "CSARP_post_sar_updateWF" \
                    "CSARP_post_sar_noise_changedWG" \
                    "CSARP_qlook" \
		    "CSARP_qlook_deconv" \
                    "CSARP_qlook_deconv_AJ" \
		    "CSARP_qlook_deconv_changedWF" \
		    "CSARP_qlook_deconv_updated" \
		    "CSARP_qlook_noise" \
		    "CSARP_qlook_noise_changedWF" \
                    "CSARP_standard" \
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

frames=/work/ollie/rzindler/Scratch/csarp_support/frames/snow/2018_Greenland_Polar6
records=/work/ollie/rzindler/Scratch/csarp_support/records/snow/2018_Greenland_Polar6
gps=/work/ollie/rzindler/Scratch/csarp_support/gps/2018_Greenland_Polar6
metadata=/work/ollie/rzindler/Scratch/metadata/2018_Greenland_Polar6

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












