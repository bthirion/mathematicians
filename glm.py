"""
Implement the level-1 GLM on a subject by subject basis

Author: Bertrand Thirion, 2013
"""
import os
import glob
import numpy as np
from scipy.io import loadmat


subjects = ['aa130114', 'jl120341', 'mp130263', 'aa130169', 'jl130200',
            'mr120371', 'al130244', 'kg120369', 'nl120167', 'bm120103',
            'ld130145',  'cb120288', 'ce130459', 'rm130241', 'cf120444', 
            'll130242', 'el120268', 
            'lr120300', 'tb120212', 'fm120345', 'vb120303', 'hr120357', 
            'mh120250', 'vb120409', 'jc130030', 'mk130199'][:1]

work_dir = '/neurospin/tmp/mathematicians'
spm_dir = '/neurospin/unicog/protocols/IRMf/mathematicians_Amalric_Dehaene2012/fMRI_data/'


for subject in subjects:
    subject_dir = os.path.join(work_dir, subject)
    fmri_dir = os.path.join(subject_dir, 'fmri')
    analysis_dir = os.path.join(spm_dir, subject, 'analyses')
    
    # audiosentence protocol
    # step1: create the paradigm
    onset_dir = os.path.join(analysis_dir, 'audiosentence')
    onset_files = glob.glob(os.path.join(onset_dir, 'onsetfile*.mat'))
    motion_param = glob.glob(
        os.path.join(spm_dir, 'fMRI/audiosentence/rp*.txt'))
    onset_files.sort()
    motion_param.sort()
    for onset_file in onset_files:
        paradigm_data = loadmat(onset_file)
        durations = np.concatenate(
            [x.ravel() for x in paradigm_data['durations'][0]])
        onsets = np.concatenate(
            [x.ravel() for x in paradigm_data['onsets'][0]])
        names = np.concatenate(
            [x.ravel() for x in paradigm_data['names'][0]]).astype(str)
        names = np.concatenate((names[:-3], np.repeat(names[-3:], 15)))

    # step 2: create the design matrix
    
        
    # step 3: fit the GLM
