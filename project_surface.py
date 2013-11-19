"""
This script project the coregistered fMRI images to the surface:
the surface is the grey-white matter interface of the subject

The purpose is to perform proper group analysis on the surface on fsaverage,
and use existing  atlases on the surface.

Author: Bertrand Thirion, 2013

Note
----
First run: export SUBJECTS_DIR=''
"""
import os
import glob
import commands
from nibabel import load, save, Nifti1Image

FWHM = 5.

# get subject list and subject informations
subjects = ['aa130114', 'jl120341', 'mp130263', 'aa130169', 'jl130200',
            'mr120371', 'al130244', 'kg120369', 'nl120167', 'bm120103',
            'ld130145',  'cb120288', 'ce130459', 'rm130241', 'cf120444', 
            'll130242', 'el120268', 
            'lr120300', 'tb120212', 'fm120345', 'vb120303', 'hr120357', 
            'mh120250', 'vb120409', 'jc130030', 'mk130199'][:1]

work_dir = '/neurospin/tmp/mathematicians'

for subject in subjects:
    subject_dir = os.path.join(work_dir, subject)
    t1_dir = os.path.join(subject_dir, 't1')
    fmri_dir = os.path.join(subject_dir, 'fmri')
    # anat_image = glob.glob(os.path.join(t1_dir, 'anat*.nii'))[0]
    fmri_images = glob.glob(os.path.join(fmri_dir, 'crvisu*.nii'))
    fmri_images += glob.glob(os.path.join(fmri_dir, 'craudio*.nii.gz'))
    fmri_images += glob.glob(os.path.join(fmri_dir, 'crloc*.nii.gz'))

    fs_dir = os.path.join(t1_dir, subject)

    # --------------------------------------------------------------------
    # run the projection using freesurfer
            
    # take the fMRI series
    for fmri_session in fmri_images:
        # output names
        # the .gii files will be put in the same directory as the input fMRI
        left_fmri_tex = fmri_session[:-4] + '_lh.gii' 
        right_fmri_tex = fmri_session[:-4] + '_rh.gii'
        print left_fmri_tex, right_fmri_tex

        # run freesrufer command
        commands.getoutput(
            '$FREESURFER_HOME/bin/mri_vol2surf --src %s --o %s '\
                '--out_type gii --regheader %s --hemi lh --projfrac 0.5'
            % (fmri_session, left_fmri_tex, fs_dir))

        plop = commands.getoutput(
            '$FREESURFER_HOME/bin/mri_vol2surf --src %s --o %s '\
                '--out_type gii --regheader %s --hemi rh --projfrac 0.5'
            % (fmri_session, right_fmri_tex, fs_dir))
  

