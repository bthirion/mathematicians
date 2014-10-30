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
    visualcategs_dmtx, visualcategs_contrasts, make_ratings)

subjects = ['cf120444','jl120341','lr120300','aa130114','aa130169','mk130199',
            'jl130200','mp130263','rm130241','al130244','bm120103','ce130459',
            'of140017','jf140025','cr140040','fm120345','hr120357','kg120369',
            'mr120371','jc130030','ld130145','cf140022','jn140034','mv140024',
            'tj140029','ap140030','af140169','pp140165','eb140248','gq140243']
subjects = subjects[:1]

work_dir = '/neurospin/tmp/mathematicians'
spm_dir = os.path.join('/neurospin/unicog/protocols/IRMf', 
                       'mathematicians_Amalric_Dehaene2012/fMRI_data/')
behavioral_dir = '/neurospin/unicog/protocols/IRMf/mathematicians_Amalric_Dehaene2012/behavioral_data/'

# some fixed parameters
tr = 1.5 # TR
    
for subject in subjects:
    # necessary paths
    analysis_dir = os.path.join(spm_dir, subject, 'analyses')
    subject_dir = os.path.join(work_dir, subject)
    fmri_dir = os.path.join(spm_dir, subject, 'fMRI')
    result_dir = os.path.join(work_dir, subject, 'results')
    if os.path.exists(result_dir) == False:
        os.mkdir(result_dir)
    memory = Memory(cachedir=os.path.join(fmri_dir, 'cache_dir'), verbose=0)

    # audiosentence protocol
    # step 1: get the necessary files
    final_data = os.path.join(behavioral_dir, subject,
                               'finaldata_%s.mat' %subject)
    ratings = make_ratings(final_data)
    onset_dir = os.path.join(analysis_dir, 'audiosentence')
    onset_files = glob.glob(os.path.join(onset_dir, 'onsetfile*.mat'))
    motion_files = glob.glob(
        os.path.join(spm_dir, subject, 'fMRI/audiosentence/rp*.txt'))
    fmri_files = glob.glob(os.path.join(fmri_dir,
                                        'audiosentence/waaudio*.nii'))
    onset_files.sort()
    motion_files.sort()
    fmri_files.sort()
    
    # scan times
    n_scans = 200

    # mask image
    make_mask = memory.cache(make_mask)
    mask = make_mask(fmri_files)
    save(mask, os.path.join(fmri_dir, 'mask.nii'))
    
    design_matrices = []
    for i, (onset_file, motion_file, fmri_file) in enumerate(zip(
            onset_files, motion_files, fmri_files)):
        # Create the design matrix
        dmtx = audiosentence_dmtx(final_data, motion_file, n_scans, tr, i)
        design_matrices.append(dmtx.matrix)
        session_contrasts = audiosentence_contrasts(dmtx.names, final_data, i)
        if i == 0:
            effects = dict([(con_id, [])
                            for con_id in session_contrasts.keys()])
            variances = dict([(con_id, [])
                            for con_id in session_contrasts.keys()])
     
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
    
    #########################################################################
    # localizer protocol
    # get the necessary files
    motion_file, = glob.glob(
        os.path.join(fmri_dir, 'localizer/rp*.txt'))
    fmri_file = glob.glob(os.path.join(fmri_dir,
                                       'localizer/walocalizer*.nii'))[0]
    n_scans = 205

    # Create the design matrix
    dmtx = localizer_dmtx(motion_file, n_scans, tr)
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

    #########################################################################
    # visualcategs protocol
    onset_dir = os.path.join(analysis_dir, 'visualcategs')
    onset_files = glob.glob(os.path.join(onset_dir, 'onsetfile*.mat'))
    motion_files = glob.glob(
        os.path.join(fmri_dir, 'visualcategs/rp*.txt'))
    fmri_files = glob.glob(os.path.join(fmri_dir, 'visualcategs/wavisu*.nii'))
    onset_files.sort()
    motion_files.sort()
    fmri_files.sort()
    
    # scan times
    n_scans = 185
    design_matrices = []
    for i, (onset_file, motion_file, fmri_file) in enumerate(zip(
        onset_files, motion_files, fmri_files)):
        # Create the design matrix
        dmtx = visualcategs_dmtx(onset_file, motion_file, n_scans, tr)
        design_matrices.append(dmtx.matrix)
        session_contrasts = visualcategs_contrasts(dmtx.names)
        if i == 0:
            effects = dict([(con_id, [])
                            for con_id in session_contrasts.keys()])
            variances = dict([(con_id, [])
                              for con_id in session_contrasts.keys()])
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
