"""
Implement the level-1 GLM on a subject by subject basis

For the moment, this is done in the volume

todo
----
* add more contrasts
* fixed effects
* replace motion files by those of pypreprocess

Author: Bertrand Thirion, 2013
"""
import os
import glob
import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
from nibabel import load, save, Nifti1Image
from joblib import Memory

from nipy.modalities.fmri.experimental_paradigm import BlockParadigm
from nipy.modalities.fmri.design_matrix import make_dmtx
from nipy.modalities.fmri.glm import FMRILinearModel

from utils import (
    audiosentence_paradigm, audiosentence_dmtx, audiosentence_contrasts,
    fixed_effects_img, make_mask, localizer_dmtx, localizer_contrasts,
    visualcategs_dmtx, visualcategs_contrasts)

subjects = ['aa130114', 'jl120341', 'mp130263', 'aa130169', 'jl130200',
            'mr120371', 'al130244', 'kg120369', 'nl120167', 'bm120103',
            'ld130145',  'cb120288', 'ce130459', 'rm130241', 'cf120444', 
            'll130242', 'el120268', 
            'lr120300', 'tb120212', 'fm120345', 'vb120303', 'hr120357', 
            'mh120250', 'vb120409', 'jc130030', 'mk130199'][:1]

work_dir = '/neurospin/tmp/mathematicians'
spm_dir = os.path.join('/neurospin/unicog/protocols/IRMf', 
                       'mathematicians_Amalric_Dehaene2012/fMRI_data/')

# some fixed parameters
tr = 1.5 # TR
    
for subject in subjects:
    # necessary paths
    analysis_dir = os.path.join(spm_dir, subject, 'analyses')
    subject_dir = os.path.join(work_dir, subject)
    fmri_dir = os.path.join(subject_dir, 'fmri')
    result_dir = os.path.join(fmri_dir, 'results')
    if os.path.exists(result_dir) == False:
        os.mkdir(result_dir)
    memory = Memory(cachedir=os.path.join(fmri_dir, 'cache_dir'), verbose=0)

    # audiosentence protocol
    # step 1: get the necessary files
    onset_dir = os.path.join(analysis_dir, 'audiosentence')
    onset_files = glob.glob(os.path.join(onset_dir, 'onsetfile*.mat'))
    motion_files = glob.glob(
        os.path.join(spm_dir, subject, 'fMRI/audiosentence/rp*.txt'))
    fmri_files = glob.glob(os.path.join(fmri_dir, 'craudio*.nii.gz'))
    onset_files.sort()
    motion_files.sort()
    fmri_files.sort()
    
    # scan times
    n_scans = 200

    # mask image
    make_mask = memory.cache(make_mask)
    mask = make_mask(fmri_files)
    save(mask, os.path.join(fmri_dir, 'mask.nii'))
    """
    effects = {'visual':[], 'audio':[], 'reflection':[], 'motor':[]}
    variances = {'visual':[], 'audio':[], 'reflection':[], 'motor':[]}
    design_matrices = []
    for (onset_file, motion_file, fmri_file) in zip(
        onset_files, motion_files, fmri_files):
        # Create the design matrix
        dmtx = audiosentence_dmtx(onset_file, motion_file, n_scans, tr)
        ax = dmtx.show()
        ax.set_position([.05, .25, .9, .65])
        ax.set_title('Design matrix')
        design_matrices.append(dmtx.matrix)
        session_contrasts = audiosentence_contrasts(dmtx.names)
         
        # fit the GLM
        fmri_glm = FMRILinearModel(fmri_file, dmtx.matrix, mask=mask)
        fmri_glm.fit(do_scaling=True, model='ar1')
        # Estimate the contrasts
        print('Computing contrasts...')
        for index, contrast_id in enumerate(session_contrasts):
            print('  Contrast % i out of %i: %s' %
                  (index + 1, len(session_contrasts), contrast_id))
            # save the z_image
            con, var = fmri_glm.contrast(
                session_contrasts[contrast_id], con_id=contrast_id, 
                output_z=False, output_effects=True, output_variance=True)
            effects[contrast_id].append(con)
            variances[contrast_id].append(var)
        del fmri_glm
    
    for index, contrast_id in enumerate(session_contrasts):
        _, _, z_map = fixed_effects_img(
            effects[contrast_id], variances[contrast_id], mask)
        z_map_path = os.path.join(result_dir, '%s_z_map.nii' % contrast_id)
        save(z_map, z_map_path)
    """
    
    #########################################################################
    # localizer protocol
    """
    # get the necessary files
    motion_file, = glob.glob(
        os.path.join(spm_dir, subject, 'fMRI/localizer/rp*.txt'))
    fmri_file = glob.glob(os.path.join(fmri_dir, 'crlocalizer*.nii.gz'))[0]
    n_scans = 205

    # Create the design matrix
    dmtx = localizer_dmtx(motion_file, n_scans, tr)
    ax = dmtx.show()
    ax.set_position([.05, .25, .9, .65])
    ax.set_title('Design matrix')
    session_contrasts = localizer_contrasts(dmtx)

    fmri_glm = FMRILinearModel(fmri_file, dmtx.matrix, mask=mask)
    fmri_glm.fit(do_scaling=True, model='ar1')
    # Estimate the contrasts
    print('Computing contrasts...')
    for index, contrast_id in enumerate(session_contrasts):
        print('  Contrast % i out of %i: %s' %
              (index + 1, len(session_contrasts), contrast_id))
        # save the z_image
        z_image, = fmri_glm.contrast(
            session_contrasts[contrast_id], con_id=contrast_id, 
            output_z=True)
        z_map_path = os.path.join(result_dir, '%s_z_map.nii' % contrast_id)
        save(z_image, z_map_path)
    del fmri_glm
    """
    #########################################################################
    # visualcategs protocol
    onset_dir = os.path.join(analysis_dir, 'visualcategs')
    onset_files = glob.glob(os.path.join(onset_dir, 'onsetfile*.mat'))
    motion_files = glob.glob(
        os.path.join(spm_dir, subject, 'fMRI/visualcategs/rp*.txt'))
    fmri_files = glob.glob(os.path.join(fmri_dir, 'crvisu*.nii.gz'))
    onset_files.sort()
    motion_files.sort()
    fmri_files.sort()
    
    # scan times
    n_scans = 185
    cnames = ['symbols-rest', 'pictures-rest', 'symbols-pictures']
    effects = {key:[] for key in cnames}
    variances = {key:[] for key in cnames}
    design_matrices = []
    for (onset_file, motion_file, fmri_file) in zip(
        onset_files, motion_files, fmri_files):
        # Create the design matrix
        dmtx = visualcategs_dmtx(onset_file, motion_file, n_scans, tr)
        ax = dmtx.show()
        ax.set_position([.05, .25, .9, .65])
        ax.set_title('Design matrix')
        design_matrices.append(dmtx.matrix)
        session_contrasts = visualcategs_contrasts(dmtx.names)
         
        # fit the GLM
        fmri_glm = FMRILinearModel(fmri_file, dmtx.matrix, mask=mask)
        fmri_glm.fit(do_scaling=True, model='ar1')
        # Estimate the contrasts
        print('Computing contrasts...')
        for index, contrast_id in enumerate(session_contrasts):
            print('  Contrast % i out of %i: %s' %
                  (index + 1, len(session_contrasts), contrast_id))
            # save the z_image
            con, var = fmri_glm.contrast(
                session_contrasts[contrast_id], con_id=contrast_id, 
                output_z=False, output_effects=True, output_variance=True)
            effects[contrast_id].append(con)
            variances[contrast_id].append(var)
        del fmri_glm
    
    for index, contrast_id in enumerate(session_contrasts):
        _, _, z_map = fixed_effects_img(
            effects[contrast_id], variances[contrast_id], mask)
        z_map_path = os.path.join(result_dir, '%s_z_map.nii' % contrast_id)
        save(z_map, z_map_path)
        
    
plt.show()
