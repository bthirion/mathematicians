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
# from nibabel import load, save, Nifti1Image
import gzip


FWHM = 5.

# get subject list and subject informations
subjects = ['cf120444','jl120341','lr120300','aa130114','aa130169','mk130199',
            'jl130200','mp130263','rm130241','al130244','bm120103','ce130459',
            'of140017','jf140025','cr140040','fm120345','hr120357','kg120369',
            'mr120371','jc130030','ld130145','cf140022','jn140034','mv140024',
            'tj140029','ap140030','af140169','pp140165','eb140248','gq140243'][1:]         

#work_dir = '/neurospin/tmp/mathematicians'

work_dir = '/neurospin/unicog/protocols/IRMf/mathematicians_Amalric_Dehaene2012/Surface_analysis/mathematicians'
spm_dir = os.path.join('/neurospin/unicog/protocols/IRMf', 
                       'mathematicians_Amalric_Dehaene2012/fMRI_data/')

os.environ['SUBJECTS_DIR'] = ""

for subject in subjects:
    print "Subject :", subject
    subject_dir = os.path.join(work_dir, subject)
    t1_dir = os.path.join(subject_dir, 't1')
    fmri_dir = os.path.join(subject_dir, 'fmri')
#    preproc_dir = os.path.join(fmri_dir, 'tmp')
    preproc_dir = os.path.join(spm_dir,subject,'fMRI')
    
#    fmri_images = glob.glob(os.path.join(fmri_dir, 'crvisu*.nii.gz'))
#    fmri_images = glob.glob(os.path.join(preproc_dir, 'craudio*.nii.gz'))
#    fmri_images += glob.glob(os.path.join(preproc_dir, 'crloc*.nii.gz'))    # old names

#    fmri_images = glob.glob(os.path.join(preproc_dir, 'Coregvisu*.nii.gz'))
#    fmri_images += glob.glob(os.path.join(preproc_dir, 'Coregaudio*.nii.gz'))
#    fmri_images += glob.glob(os.path.join(preproc_dir, 'Coregloc*.nii.gz'))   # new names

    fmri_images = glob.glob(os.path.join(preproc_dir,'visualcategs/wavisu*.nii'))
    fmri_images += glob.glob(os.path.join(preproc_dir,'audiosentence/waaudio*.nii'))
    fmri_images += glob.glob(os.path.join(preproc_dir,'localizer/walocalizer*.nii'))


    fs_dir = os.path.join(t1_dir, subject)

    # --------------------------------------------------------------------
    # run the projection using freesurfer
            
    # take the fMRI series
    for fmri_session in fmri_images:
        # output names
        # the .gii files will be put in the same directory as the input fMRI
        left_fmri_tex = fmri_session[:-7] + '_lh.gii' 
        right_fmri_tex = fmri_session[:-7] + '_rh.gii'
        print left_fmri_tex, right_fmri_tex

        # unzip the fMRI data
        fmri_file = fmri_session[:-3]
#        f_in = gzip.open(fmri_session, 'rb')
        f_in = open(fmri_session, 'rb')
        f_out = open(fmri_file, 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()

        # run freesrufer command
        commands.getoutput(
            '$FREESURFER_HOME/bin/mri_vol2surf --src %s --o %s '\
                '--out_type gii --regheader %s --hemi lh --projfrac 0.5'
            % (fmri_session, left_fmri_tex, fs_dir))

        plop = commands.getoutput(
            '$FREESURFER_HOME/bin/mri_vol2surf --src %s --o %s '\
                '--out_type gii --regheader %s --hemi rh --projfrac 0.5'
            % (fmri_session, right_fmri_tex, fs_dir))
        
        # delete the nii file
        os.remove(fmri_file)
     
