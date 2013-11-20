"""
Implement the level-1 GLM on a subject by subject basis on the cortical surface

Todo: both hemispheres

Author: Bertrand Thirion, 2013
"""
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from joblib import Memory

from nipy.modalities.fmri.glm import GeneralLinearModel
from utils import (
    audiosentence_paradigm, audiosentence_dmtx, audiosentence_contrasts,
    fixed_effects, make_mask)
from nibabel.gifti import read, write, GiftiDataArray, GiftiImage

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
    fmri_files = glob.glob(os.path.join(fmri_dir, 'craudio*_lh.gii'))
    onset_files.sort()
    motion_files.sort()
    fmri_files.sort()
    
    # scan times
    n_scans = 200
    
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
         
        # load the data
        Y = np.array([darrays.data for darrays in read(fmri_file).darrays])
        # fit the GLM
        fmri_glm = GeneralLinearModel(dmtx.matrix)
        fmri_glm.fit(Y, model='ar1')
        # Estimate the contrasts
        print('Computing contrasts...')
        for index, contrast_id in enumerate(session_contrasts):
            print('  Contrast % i out of %i: %s' %
                  (index + 1, len(session_contrasts), contrast_id))
            # save the z_image
            contrast_ = fmri_glm.contrast(session_contrasts[contrast_id])
            effects[contrast_id].append(contrast_.effect.ravel())
            variances[contrast_id].append(contrast_.variance.ravel())
        del fmri_glm
    
    for index, contrast_id in enumerate(session_contrasts):
        _, _, z_map = fixed_effects(
            effects[contrast_id], variances[contrast_id])
        z_texture = GiftiImage(
            darrays=[GiftiDataArray().from_array(z_map, intent='t test')])
        z_map_path = os.path.join(result_dir, '%s_z_map_lh.gii' % contrast_id)
        write(z_texture, z_map_path)
        
        
    
plt.show()
