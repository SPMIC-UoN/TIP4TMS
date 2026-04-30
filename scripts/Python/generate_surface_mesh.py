#
# Author: Stefan Pszczolkowski P. <mszspp@nottingham.ac.uk>
#
import nipype.interfaces.freesurfer as freesurf
import argparse
import os
import subprocess


# Call Freesurfer's mkheadsurf (which, for some reason, is not present in nipype)
def create_head_surface(subject_id, subjects_dir, nsmooth):
    return subprocess.call(['mkheadsurf', '-s', subject_id, '-sd', subjects_dir, '-nsmooth', '%d' % nsmooth])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create StimGuide coordinate file for the specified subject', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('subject_id', help='Subject ID')
    parser.add_argument('t1_image_file', help='Input T1 image file')
    parser.add_argument('--subjects-dir', dest='subjects_dir', default=None, help='Specify alternative FreeSurfer subjects directory')
    parser.add_argument('--no-bias-corr', dest='no_bias_corr', default=False, action='store_true', help='Do not perform bias correction')

    args = parser.parse_args()
    
    t1_image_file = args.t1_image_file
    subject_id = args.subject_id
    subjects_dir = args.subjects_dir
    no_bias_corr = args.no_bias_corr

    # Default SUBJECTS_DIR is the environ variable SUBJECTS_DIR or, if this is unavailable, a directory called 'subjects' in the user's home directory
    if subjects_dir is None:
        subjects_dir = os.environ['SUBJECTS_DIR']
        if subjects_dir is None:
            subjects_dir = os.path.join(os.environ['HOME'], 'subjects')

    subject_root_dir = os.path.join(subjects_dir, subject_id)
    subject_mri_dir = os.path.join(subject_root_dir, 'mri')
    subject_surf_dir = os.path.join(subject_root_dir, 'surf')

    if not os.path.exists(subjects_dir):
        os.makedirs(subjects_dir)
        
    if not os.path.exists(subject_root_dir):
        os.makedirs(subject_root_dir)
        
    if not os.path.exists(subject_mri_dir):
        os.makedirs(subject_mri_dir)

    if not os.path.exists(subject_surf_dir):
        os.makedirs(subject_surf_dir)
        
    mgz_t1_file = os.path.join(subject_mri_dir, 'T1.mgz')
    fs_scalp_file = os.path.join(subject_surf_dir, 'lh.seghead')
    gifti_scalp_file = os.path.join(subject_surf_dir, 'lh.seghead.surf.gii')

    # Convert input image to Freesurfer's MGZ
    mri_conv = freesurf.MRIConvert()
    mri_conv.inputs.in_file = t1_image_file
    mri_conv.inputs.out_file = mgz_t1_file
    mri_conv.inputs.out_type = 'mgz'
    mri_conv.run()

    if not no_bias_corr:
        # Perform bias field correction
        bias_corr = freesurf.MNIBiasCorrection()
        bias_corr.inputs.in_file = mgz_t1_file
        bias_corr.inputs.out_file = mgz_t1_file
        bias_corr.run()

    # Compute surface
    create_head_surface(subject_id, subjects_dir, 100)

    # Convert resulting Freesurfer surface into GIFTI format
    mris_conv = freesurf.MRIsConvert()
    mris_conv.inputs.in_file = fs_scalp_file
    mris_conv.inputs.out_file = gifti_scalp_file
    mris_conv.run()