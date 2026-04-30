#!/bin/bash

export SelfDir=$(dirname "$(readlink -f "$0")")

Usage()
{
  echo "Usage: `basename $0` [OPTIONS]"
  echo " "
  echo " "
  echo "OPTIONS:"
  echo " "
  echo " "
  echo " --path <path>                      (Mandatory) BRC Pipeline subjects path"
  echo " --subject <subject id>             (Mandatory) Subject ID"
  echo " --seed-mask <file>                 (Mandatory) Nifti file specifying the seed mask in MNI space" 
  echo " --landmarks-file <file>            (Mandatory) Text file containing the coordinates in mm of the nasion and preauricular points"
  echo " --searchlight-mask <file>          (Optional)  Nifti file specifying the searchlight mask in MNI space (Default: whole brain mask)"
  echo " --background-thresh <value>        (Optional)  Intensity threshold between head and background on subject T1w image (Default: compute automatically)"
  echo " --output-path <path>               (Optional)  Path of the folder where the compressed directory with StimGuide files is written (Default: 'StimGuide' folder inside subject directory)"
  echo " --fmri-folder-name <folder name>   (Optional)  Name of the fMRI folder from the BRC functional analysis pipeline. Default: rfMRI"
  echo " --map-smoothing <fwhm in mm>       (Optional)  Full-width-half-maximum in millimetres of the smoothing kernel passed through the resulting granger/correlation map. Default: 2"
  echo " --do-pos-granger                   (Optional)  Select target based on excitatory signed-path-based Granger causality (Default)"
  echo " --do-neg-granger                   (Optional)  Select target based on inhibitory signed-path-based Granger causality"
  echo " --do-correlation                   (Optional)  Select target based on most positively correlated point"
  echo " --do-anticorrelation               (Optional)  Select target based on most negatively correlated point"
  echo " -h | --help                        help"
  echo " "
  echo " "
  echo "NOTE: If multiple --do-XXXX flags are given, only the last one will be considered."
  echo " "
  echo " "
}

# Just give usage if no arguments specified
if [ $# -eq 0 ] ; then Usage; exit 0; fi

regex_integer='^[1-9][0-9]*$'
regex_float='^[0-9]+\.?[0-9]*$'

# Mandatory arguments
Path=
Sub_ID=
Seed_Mask=
Landmarks_File=
Bkgr_Thresh=

# Default values
Searchlight_Mask=
Output_Path=
fMRI_Folder_Name="rfMRI"
Map_Smoothing="2"
Do_Pos_Granger="yes"
Do_Neg_Granger="no"
Do_Correlation="no"
Do_Anticorrelation="no"

# parse arguments
while [ "$1" != "" ]
do
    case $1 in
        --path )	            shift
                                Path=$1
                                ;;

        --subject )				shift
				                Sub_ID=$1
                                ;;
								
		--seed-mask )			shift
                                Seed_Mask=$1
                                ;;
								
		--searchlight-mask )	shift
                                Searchlight_Mask=$1
                                ;;
								
		--landmarks-file )		shift
                                Landmarks_File=$1
                                ;;
								
		--background-thresh )   shift
				                Bkgr_Thresh=$1
                                ;;
								
		--output-path )	        shift
                                Output_Path=$1
                                ;;
								
		--fmri-folder-name )	shift
                                fMRI_Folder_Name=$1
                                ;;
								
		--map-smoothing )		shift
                                Map_Smoothing=$1
                                ;;
								
		--do-pos-granger )  	Do_Pos_Granger="yes"
								Do_Neg_Granger="no"
								Do_Correlation="no"
								Do_Anticorrelation="no"
                                ;;
								
		--do-neg-granger )  	Do_Pos_Granger="no"
								Do_Neg_Granger="yes"
								Do_Correlation="no"
								Do_Anticorrelation="no"
                                ;;
								
		--do-correlation )  	Do_Pos_Granger="no"
								Do_Neg_Granger="no"
								Do_Correlation="yes"
								Do_Anticorrelation="no"
                                ;;
								
		--do-anticorrelation ) 	Do_Pos_Granger="no"
								Do_Neg_Granger="no"
								Do_Correlation="no"
								Do_Anticorrelation="yes"
                                ;;
								
        -h | --help )           Usage
                                exit
                                ;;

        * )                     echo "Unknown flag: $1"
								echo
								Usage
                                exit 1
    esac
    shift
done

### Sanity checking of arguments
if [ X$Path = X ] || [ X$Sub_ID = X ] || [ X$Seed_Mask = X ] || [ X$Landmarks_File = X ]  ; then
  echo "All of the compulsory arguments --path, --subject, --seed-mask and --landmarks-file MUST be used"
  Usage
  exit 1;
fi

if [[ $Bkgr_Thresh != "" ]] ; then
	if ! [[ $Bkgr_Thresh =~ $regex_integer ]] ; then
		echo "Argument --background-thresh MUST be a valid integer"
		exit 1;
	fi
fi

if ! [[ $Map_Smoothing =~ $regex_float ]] ; then
   echo "Argument --map-smoothing MUST be a valid float number with no sign"
   exit 1;
fi

Scripts_Dir=$SelfDir/scripts
Root_Proc_Dir=$Path/$Sub_ID
Sigma=`echo "scale=5; $Map_Smoothing/2.35482" | bc`

if [ X$Searchlight_Mask = X ] ; then
	Searchlight_Mask=$SelfDir/MNI/MNI152_T1_2mm_brain_mask.nii.gz
fi

if [ X$Output_Path = X ] ; then
	Output_Path=$Root_Proc_Dir/StimGuide/
fi

### Clean previous processed data (if it exists)
rm -rf $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg
rm -rf $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms
mkdir -p $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg
mkdir -p $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms

echo -n "Preparing functional mask ... "
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_brain_mask.nii.gz -kernel sphere 3 -ero $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/func_mask.nii.gz -odt char
echo "done"

echo -n "Masking CSF and WM partial volume estimation images ... "
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF.nii.gz -mas $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/func_mask.nii.gz $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF_masked.nii.gz
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM.nii.gz -mas $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/func_mask.nii.gz $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM_masked.nii.gz
echo "done"

echo -n "Thresholding CSF and WM partial volume estimation images ... "
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF_masked.nii.gz -sub 0.98 -bin $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF_thresh.nii.gz
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM_masked.nii.gz -sub 0.98 -bin $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM_thresh.nii.gz
echo "done"

echo -n "Eroding CSF and WM partial volume estimation images ... "
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF_thresh.nii.gz -kernel sphere 2 -ero $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF_eroded.nii.gz
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM_thresh.nii.gz -kernel sphere 2 -ero $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM_eroded.nii.gz
echo "done"

echo -n "Transforming CSF and WM partial volume estimation images ... "
applywarp -i $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_CSF_eroded.nii.gz -r $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/SBref.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/CSF_highres2func.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/str2rfMRI.nii.gz --interp=nn
applywarp -i $Root_Proc_Dir/analysis/anatMRI/T1/processed/seg/tissue/sing_chan/T1_pve_WM_eroded.nii.gz -r $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/SBref.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/WM_highres2func.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/str2rfMRI.nii.gz --interp=nn
echo "done"

echo -n "Extracting WM and CSF timeseries ... "
fslmeants -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/rfMRI.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/WM_TS.txt -m $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/WM_highres2func.nii.gz
fslmeants -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/rfMRI.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/CSF_TS.txt -m $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/CSF_highres2func.nii.gz
echo "done"

echo -n "Creating regressor file ... "
paste $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/WM_TS.txt $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/CSF_TS.txt > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/WM_CSF_TS.txt
echo "done"

echo -n "Regressing out WM and CSF from fMRI data ... "
fsl_regfilt -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/rfMRI.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/rfMRI_CSFWMregressed.nii.gz -d $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/WM_CSF_TS.txt -f "1,2"
echo "done"

echo -n "Transforming seed mask ... "
applywarp -i $Seed_Mask -r $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/SBref.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/seed_mask.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/std2rfMRI.nii.gz --interp=nn
echo "done"

echo -n "Transforming searchlight mask ... "
applywarp -i $Searchlight_Mask -r $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/SBref.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/searchlight_mask.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/std2rfMRI.nii.gz --interp=nn
echo "done"

if [ $Do_Pos_Granger = "yes" ] || [ $Do_Neg_Granger = "yes" ] ; then

	${MATLABpath}/matlab -nodesktop -r "addpath(genpath('$Scripts_Dir/MATLAB')); compute_coefficient_GCA_maps('$Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/rfMRI_CSFWMregressed.nii.gz', '$Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/seed_mask.nii.gz', '$Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/searchlight_mask.nii.gz', '$Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/GCA_map_', 1); exit"
	
	echo -n "Warping Granger map into T1 space ... "
	applywarp -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/GCA_map_X2Y_1.nii.gz -r $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_unbiased.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/GCA_map_X2Y_1_T1.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/rfMRI2str.nii.gz --interp=trilinear
	echo "done"
	
	echo -n "Smoothing T1 space Granger map ... "
	fslmaths $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/GCA_map_X2Y_1_T1.nii.gz -s $Sigma $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/GCA_map_X2Y_1_T1_smooth.nii.gz
	echo "done"
	
else

	${MATLABpath}/matlab -nodesktop -r "addpath(genpath('$Scripts_Dir/MATLAB')); compute_conn_map('$Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/data/rfMRI_CSFWMregressed.nii.gz', '$Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/seed_mask.nii.gz', '$Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/seg/searchlight_mask.nii.gz', '$Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/'); exit"
	
	echo -n "Warping correlation map into T1 space ... "
	applywarp -i $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/conn_map.nii.gz -r $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_unbiased.nii.gz -o $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/conn_map_T1.nii.gz -w $Root_Proc_Dir/analysis/$fMRI_Folder_Name/preproc/reg/rfMRI2str.nii.gz --interp=trilinear
	echo "done"
	
	echo -n "Smoothing T1 space correlation map ... "
	fslmaths $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/conn_map_T1.nii.gz -s $Sigma $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/conn_map_T1_smooth.nii.gz
	echo "done"
	
fi

if [ $Do_Pos_Granger = "yes" ] ; then

	echo -n "Finding maximum point of Granger map ... "
	fslstats $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/GCA_map_X2Y_1_T1_smooth.nii.gz -l `fslstats $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/GCA_map_X2Y_1_T1_smooth.nii.gz -R | awk -F ' ' '{print $2 - 0.00001}'` -n -a -c | xargs printf "%0.2f %0.2f %0.2f\n" > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/${Sub_ID}_target_point_T1.txt
	echo "done"
	
elif [ $Do_Neg_Granger = "yes" ] ; then

	echo -n "Finding minimum point of Granger map ... "
	fslstats $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/GCA_map_X2Y_1_T1_smooth.nii.gz -u `fslstats $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/GCA_map_X2Y_1_T1_smooth.nii.gz -R | awk -F ' ' '{print $1 + 0.00001}'` -n -a -c | xargs printf "%0.2f %0.2f %0.2f\n" > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/${Sub_ID}_target_point_T1.txt
	echo "done"
	
elif [ $Do_Correlation = "yes" ] ; then

	echo -n "Finding maximum point of correlation map ... "
	fslstats $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/conn_map_T1_smooth.nii.gz -l `fslstats $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/conn_map_T1_smooth.nii.gz -R | awk -F ' ' '{print $2 - 0.00001}'` -n -a -c | xargs printf "%0.2f %0.2f %0.2f\n" > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/${Sub_ID}_target_point_T1.txt
	echo "done"
	
else

	echo -n "Finding minimum point of correlation map ... "
	fslstats $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/conn_map_T1_smooth.nii.gz -u `fslstats $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/conn_map_T1_smooth.nii.gz -R | awk -F ' ' '{print $1 + 0.00001}'` -n -a -c | xargs printf "%0.2f %0.2f %0.2f\n" > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/${Sub_ID}_target_point_T1.txt
	echo "done"
	
fi

echo -n "Creating full landmark file ... "
cat -s $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/${Sub_ID}_target_point_T1.txt $Landmarks_File > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/landmarks.txt
echo "done"

echo -n "Warping MNI head map and strip into T1 space ... "
flirt -in $SelfDir/MNI/headmask_dil.nii.gz -ref $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_unbiased.nii.gz -out $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_headmask_dil.nii.gz -applyxfm -init $Root_Proc_Dir/analysis/anatMRI/T1/preproc/reg/std_2_T1.mat -interp nearestneighbour
flirt -in $SelfDir/MNI/headstrip.nii.gz -ref $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_unbiased.nii.gz -out $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_headstrip.nii.gz -applyxfm -init $Root_Proc_Dir/analysis/anatMRI/T1/preproc/reg/std_2_T1.mat -interp nearestneighbour
echo "done"

if [ X$Bkgr_Thresh = X ]  ; then
	echo -n "Computing background threshold ... "
	Bkgr_Thresh=`fslstats $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_unbiased.nii.gz -k $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_headstrip.nii.gz -P 98 | xargs printf "%.2f"`
	echo "done"
fi

echo -n "Removing background of T1 image (threshold = $Bkgr_Thresh) ... "
fslmaths $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_unbiased.nii.gz -mas $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_headmask_dil.nii.gz -thr $Bkgr_Thresh $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_nobackground.nii.gz
echo "done"

echo -n "Creating head surface mesh ... "
python $Scripts_Dir/Python/generate_surface_mesh.py $Sub_ID $Root_Proc_Dir/analysis/anatMRI/T1/processed/data/T1_nobackground.nii.gz --subjects-dir $Path --no-bias-corr > /dev/null
echo "done"

echo -n "Creating StimGuide files ... "
python $Scripts_Dir/Python/generate_treatment_files.py $Sub_ID $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/landmarks.txt $Output_Path/tmp --subjects-dir $Path > $Root_Proc_Dir/analysis/$fMRI_Folder_Name/processed/tip4tms/target_info.txt
echo "done"

echo -n "Creating compressed files ... "
zip -q -j $Output_Path/${Sub_ID}_target.zip $Output_Path/tmp/*
rm -rf $Output_Path/tmp
echo "done"