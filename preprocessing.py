"""
Pre-process the fMRI data to allow projection onto freesurfer meshes

Author: Bertrand Thirion, 2013
"""
import os
import glob
import shutil

from nipype.caching import Memory
from nipype.interfaces import spm
import nipype.interfaces.matlab as matlab
from nibabel import load, save, Nifti1Image, four_to_three, concat_images

import resample_anat as rsa

###############################################################################
# Set the way matlab should be called at Neurospin:
# we have to hit directly the original mathworks shell script
matlab.MatlabCommand.set_default_matlab_cmd(
    "/neurospin/local/matlab/bin/matlab")
matlab.MatlabCommand.set_default_paths('/i2bm/local/spm8/')
T1_TEMPLATE = '/i2bm/local/spm8/templates/T1.nii'

###############################################################################


##############################################################
# Set the paths, dataset-specific stuff
##############################################################

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
    anat_image = glob.glob(os.path.join(t1_dir, 'anat*.nii'))
    fmri_images = glob.glob(os.path.join(fmri_dir, 'visualcategs/visu*.nii'))
    fmri_images += glob.glob(os.path.join(fmri_dir, 'audiosentence/audio*.nii'))
    fmri_images += glob.glob(os.path.join(fmri_dir, 'localizer/loc*.nii'))

    ##############################################################
    # Preprocessing
    ##############################################################

    ##############################################################
    # Anatomical segmentation (White/Grey matter)
    mem = Memory(base_dir=subject_dir)
    seg = mem.cache(spm.Segment)
    out_seg = seg(data=anat_image,
                  gm_output_type=[True, True, True],
                  wm_output_type=[True, True, True],
                  csf_output_type=[True, True, True])
    sn_file = out_seg.outputs.transformation_mat
    inv_sn_file = out_seg.outputs.inverse_transformation_mat
    gm_image = out_seg.outputs.normalized_gm_image
    native_gm_image = out_seg.outputs.native_gm_image

    shutil.copyfile(native_gm_image, os.path.join(t1_dir,
        '%s_gm_image.nii' % subject))

    ##############################################################
    # realign the data just to get the mean image
    realign = mem.cache(spm.Realign)
    realign_result = realign(
        in_files=fmri_images,
        register_to_mean=True, jobtype='estwrite', out_prefix='r')
    mean_img = realign_result.outputs.mean_image
    shutil.copyfile(mean_img,
                os.path.join(fmri_dir, os.path.basename(mean_img)))
     
    realigned_files = []
    for rf in realign_result.outputs.realigned_files:
        realigned_files.append(
            os.path.join(fmri_dir, os.path.basename(rf)))
        shutil.copyfile(rf, os.path.join(fmri_dir, os.path.basename(rf)))
    
    for rpf in realign_result.outputs.realignment_parameters:
        shutil.copyfile(rf, os.path.join(fmri_dir, os.path.basename(rpf)))

    ##############################################################
    # Coregister the anatomy to the mean functional image 
    #    (only rewrite the header)
    
    coreg = mem.cache(spm.Coregister)
    """
    # get the mean functionnal image to coregister the data on the anat
    mean_image = glob.glob(os.path.join(fmri_dir, 'mean*.nii'))[0]
    
    # does the coregistration of the time corrected+realigned fmri series
    coreg_result = coreg(source=mean_image, target=anat_image,
                         jobtype='estimate')
    """
    
    # coregister different sessions from same subject
    resampled_anat = rsa.resample_anat(t1_dir, anat_image)
    
    # get the mean functionnal image to coregister the data on the anat
    mean_image = glob.glob(os.path.join(fmri_dir, 'mean*.nii'))[0]

    ##### work around SPM bug : convert 4D images to series of 3D
    _realigned_files = []
    realigned_files = glob.glob(os.path.join(fmri_dir, '*/r*.nii'))
    realigned_files.sort()
    for rf in realigned_files:
        image_series = four_to_three(load(rf))
        for idx, nim in enumerate(image_series):
            save(nim, os.path.join(
                    '/tmp/', os.path.basename(rf)[:-4] + '_%04d.nii' % idx))
        for idx in range(len(image_series)):
            _realigned_files.append(os.path.join(
                    '/tmp/', os.path.basename(rf)[:-4] + '_%04d.nii' % idx))
    ####
    
    # does the coregistration of the time corrected+realigned fmri series
    coreg_result = coreg(target=resampled_anat, source=mean_image,
                         jobtype='estwrite', apply_to_files=_realigned_files,
                         out_prefix='c')
        
    ####
    # copy coregistered mean image used to coregister the fmri series
    coreg_source = coreg_result.outputs.coregistered_source
    coregistered_files = []
    for key in []:
        cr_files = [f for f in coreg_result.outputs.coregistered_files
                    if key in f]
        cr_files.sort()
        cr_file = concat_images(cr_files)
        cr_file_path = os.path.join(
            fmri_dir, os.path.basename(cr_files[0][:-9]) + '.nii') 
        coregistered_files.append(cr_file_path)
        save(cr_file, cr_file_path)
        del cr_file
    ####
