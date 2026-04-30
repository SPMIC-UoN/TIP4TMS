# BRIGhTMIND Trial Processing Pipeline

## Table of Contents

* [Introduction](#introduction)
* [Prerequisities](#prerequisities)
* [Script Usage](#script-usage)
* [Landmark Definition](#landmark-definition)
* [Background Threshold Definition](#background-threshold-definition)
* [Register Participant and Target in StimGuide™](#register-participant-and-target-in-stimguide)

<a id="introduction"></a>
## Introduction

The Target Identification Pipeline for Transcranial Magnetic Stimulation (TIP4TMS) is a simplification of the [pipeline](https://github.com/SPMIC-UoN/brightmind_pipeline) used in the BRIGhTMIND study (e.g. no global signal regression is implemented and only one imaging-based target based on fMRI connectivity is generated).

In contrast to the aforementioned pipeline, TIP4TMS can generate a taget based on either functional or effective connectivity.

<a id="prerequisities"></a>
## Prerequisities

TIP4TMS has the following software requirements:

1. A 64-bit Linux Operating System

2. [BRC Pipeline](https://github.com/SPMIC-UoN/BRC_Pipeline).
	* The BRC Pipeline has also its own software requirements.

3. MATLAB R2017b or later.

4. Python version 2.7 (the scripts should in principle also work with modern versions of Python (e.g. Python 3.x), but this has not been tested) with the following packages:
	* transforms3d
	* nibabel
	* numpy
	* nipype

-----

<a id="script-usage"></a>
## Script Usage

Before running this script, the structural and resting-state fMRI datasets of the subject to process must be run through the BRC Pipeline's structural `BRC_structural_pipeline/struc_preproc.sh` and functional `BRC_functional_pipeline/fMRI_preproc.sh` procedures, respectively (in that order).

The following calls to BRC Pipeline were used in the BRIGhTMIND study:

* Structural pipeline (replace italicised parameters with actual values according to your data)
  
<pre>
struc_preproc.sh --input <em>&lt;path_to_T1w_nifti&gt;</em> \
		 --path <em>&lt;brc_output_path&gt;</em> \
		 --subject <em>&lt;subject id&gt;</em> \
		 --strongbias \
		 --nodefacing
</pre>
	
* Functional pipeline (replace italicised parameters with actual values according to your data)

<pre>
fMRI_preproc.sh --fmripath <em>&lt;path_to_rsfMRI_nifti&gt;</em> \
		--path <em>&lt;brc_output_path&gt;</em> \
		--subject <em>&lt;subject id&gt;</em> \
		--mctype MCFLIRT6 \
		--dcmethod TOPUP \
		--slice2vol \
		--slspec <em>&lt;dcm2niix_json_file_for_rsfMRI&gt;</em> \
		--fwhm 5 \
		--echospacing <em>&lt;effective_echo_spacing_of_acquisition&gt;</em> \
		--stcmethod 1 \
		--SEPhasePos <em>&lt;path_to_rsfMRI_anterior_to_posterior_blip_nifti&gt;</em> \
		--SEPhaseNeg <em>&lt;path_to_rsfMRI_posterior_to_anterior_blip_nifti&gt;</em> \
		--unwarpdir y- \
		--biascorrection SEBASED \
		--echodiff NONE \
		--fmapmag NONE \
		--tempfilter 100
</pre>

-----

The usage of tip4tms.sh is as follows:

        tip4tms.sh <arguments>
		 
* Compulsory arguments (You MUST set all of):

		--path <path>                      BRC Pipeline subjects path
        --subject <subject id>             Subject ID as set in BRC Pipeline
        --seed-mask <file>                 Nifti file specifying the seed mask in MNI space 
        --landmarks-file <file>            Text file containing the coordinates in mm of the nasion and preauricular points

* Optional arguments (You may optionally specify one or more of):

        --searchlight-mask <file>          Nifti file specifying the searchlight mask in MNI space (Default: whole brain mask)
        --background-thresh <value>        Intensity threshold between head and background on subject T1w image (Default: compute automatically)
        --output-path <path>               Path of the folder where the compressed directory with StimGuide files is written (Default: 'StimGuide' folder inside subject directory)
        --fmri-folder-name <folder name>   Name of the fMRI folder from the BRC functional analysis pipeline. Default: rfMRI
        --map-smoothing <fwhm in mm>       Full-width-half-maximum in millimetres of the smoothing kernel passed through the resulting granger/correlation map. Default: 2
        --do-pos-granger                   Select target based on excitatory signed-path-based Granger causality (Default)
        --do-neg-granger                   Select target based on inhibitory signed-path-based Granger causality
        --do-correlation                   Select target based on most positively correlated point
        --do-anticorrelation               Select target based on most negatively correlated point
        -h | --help                        help

        **NOTE:** If multiple --do-XXXX flags are given, only the last one will be considered.
		
-----

<a id="landmark-definition"></a>
## Landmark Definition

The BRIGhTMIND processing pipeline requires as input a text file with the location in mm of the nasion, left preauricular point and right preauricular point. In order to create this file, open the file `<brc_output_path>/<subject id>/analysis/anatMRI/T1/processed/data/T1.nii.gz` in your favourite NIFTI image viewer.

### Location of the nasion

Locate the nasion as shown in the image below and take note of the X, Y and Z coordinates **in mm**

![Nasion Landmark](https://github.com/SPMIC-UoN/TIP4TMS/blob/main/img/nasion_landmark.png)

### Location of the left/right preauricular points

Locate the left/right preauricular points as shown in the image below and take note of the X, Y and Z coordinates **in mm**

![Preauricular Landmark](https://github.com/SPMIC-UoN/TIP4TMS/blob/main/img/preauricular_landmark.png)

### Landmark File

The final landmark file should look like the image below. The first line corresponds to the nasion coordinates, the second line to the left preauricular point coordinates and the third line to the right preauricular point coordinates.

![Landmark Coordinates](https://github.com/SPMIC-UoN/TIP4TMS/blob/main/img/landmark_coords.png)

-----

<a id="background-threshold-definition"></a>
## Background Threshold Definition

Another mandatory input to the BRIGhTMIND pipeline is an integer number defining the threshold between background and foreground in the T1-weighted structural image. In general terms, the procedure to determine this number is as follows:

1. Open the file `<brc_output_path>/<subject id>/analysis/anatMRI/T1/processed/data/T1.nii.gz` in your favourite NIFTI image viewer.

2. Window the image by lowering the maximum intensity, such that the skull/scalp has a clear boundary with the background, but at the same time there is not too much noise near the skull/scalp like in the example below:

	![Background Threshold Setting](https://github.com/SPMIC-UoN/TIP4TMS/blob/main/img/background_thresh_setting.png)
	
3. The maximum intensity selected is the number to input into the script.

After the pipeline is run, the resulting GIFTI mesh `<brc_output_path>/<subject id>/surf/lh.seghead.surf.gii` should look smooth and complete like in the example below:

![Smooth Scalp Mesh](https://github.com/SPMIC-UoN/TIP4TMS/blob/main/img/scalp_smooth.png)

If the threshold number is set too low, the resulting scalp mesh may become too noisy like in this example:

![Rough Scalp Mesh](https://github.com/SPMIC-UoN/TIP4TMS/blob/main/img/scalp_rough.png)

If the threshold number is set too high, the resulting scalp mesh may be not complete and have holes in it.

<a id="register-participant-and-target-in-stimguide"></a>
## Register Participant and Target in StimGuide™

1. Open StimGuide™
2. Register a new patient by clicking on "Register New Patient" and typing the same Subject ID used for the pipeline as Patient Identifier
3. Close StimGuide™
4. Create a "Targets" folder inside `Documents\EzGuide\Studies\<subject id>\`, where the "Documents" folder corresponds to the Microsoft Windows "Documents" folder of the user installed on the StimGuide™ computer.
5. Find the pipeline output folder (by default `<brc_output_path>/<subject id>/StimGuide`) and copy the ZIP file to a USB stick or similar.
6. Plug the USB stick on the StimGuide™ computer and unzip the contents of the ZIP file (one `.asc` and one `.xml` file) into `Documents\EzGuide\Studies\<subject id>\Targets`.
