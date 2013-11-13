"""
Implement the level-1 GLM on a subject by subject basis

For the moment, this is done in the volume

todo
----
* simplify, separate out code 
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

    
def audiosentence_paradigm(onset_file):
    """utility for audiosentence paradigm creation from the matlab onset file"""
    paradigm_data = loadmat(onset_file)
    durations = np.concatenate(
        [x.ravel() for x in paradigm_data['durations'][0]])
    onsets = np.concatenate(
        [x.ravel() for x in paradigm_data['onsets'][0]])
    names = np.concatenate(
        [x.ravel() for x in paradigm_data['names'][0]]).astype(str)
    names = np.concatenate((names[:-3], np.repeat(names[-3:], 16)))
    return BlockParadigm(names, onsets, durations)


def audiosentence_dmtx(onset_file, motion_file, n_scans):
    """Utility for ceating a design matrix from onset and motion files"""
    # Some default parameters
    hrf_model = 'canonical'  # hemodynamic reponse function
    drift_model = 'cosine'   # drift model 
    hfcut = 128              # low frequency cut
    motion_names = ['tx', 'ty', 'tz', 'rx', 'ry', 'rz'] 
    # motion param identifiers
    frametimes = np.linspace(0, (n_scans - 1) * tr, n_scans)

    paradigm = audiosentence_paradigm(onset_file)
    
    # add motion regressors and low frequencies
    # and create the design matrix
    motion_params = np.loadtxt(motion_file)
    dmtx = make_dmtx(frametimes, paradigm, hrf_model=hrf_model,
                     drift_model=drift_model, hfcut=hfcut,
                     add_regs=motion_params, add_reg_names=motion_names)
    return dmtx

def audiosentence_contrasts(names_):
    """Create the contrasts,
    given the names of the columns of the design matrix"""
    contrasts = {}
    audio = np.array([name[-4:] == 'reg1' for name in names_], np.float)
    contrasts['audio'] = audio / audio.sum()
    visual = np.array([name[-6:] == 'signal' for name in names_], np.float)
    contrasts['visual'] = visual/ visual.sum()
    motor = np.array([name == 'key press' for name in names_], np.float)
    contrasts['motor'] = motor
    reflection = np.array([name[-4:] == 'reg2' for name in names_],
                          np.float)
    contrasts['reflection'] = reflection / reflection.sum()
    return contrasts


def make_mask(fmri_files):
    """ Generate a mask from a set of fMRI files"""
    from nipy.labs.mask import compute_mask
    mean = None
    for fmri_file in fmri_files:
        if mean == None:
            mean = load(fmri_file).get_data().mean(-1)
            affine = load(fmri_file).get_affine()
        else:
            mean += load(fmri_file).get_data().mean(-1)

    mask_img = Nifti1Image(compute_mask(mean, opening=3).astype(np.uint8),
                           affine)
    return mask_img



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
    
    for (onset_file, motion_file, fmri_file) in zip(
        onset_files, motion_files, fmri_files):
        # Create the design matrix
        dmtx = audiosentence_dmtx(onset_file, motion_file, n_scans)
        ax = dmtx.show()
        ax.set_position([.05, .25, .9, .65])
        ax.set_title('Design matrix')

        contrasts = audiosentence_contrasts(dmtx.names)
        
        # step 5: fit the GLM
        fmri_glm = FMRILinearModel(fmri_file, dmtx.matrix, mask=mask)
        fmri_glm.fit(do_scaling=True, model='ar1')
    
        # Estimate the contrasts
        print('Computing contrasts...')
        for index, contrast_id in enumerate(contrasts):
            print('  Contrast % i out of %i: %s' %
                  (index + 1, len(contrasts), contrast_id))
            # save the z_image
            z_map_path = os.path.join(result_dir, '%s_z_map.nii' % contrast_id)
            z_map, = fmri_glm.contrast(
                contrasts[contrast_id], con_id=contrast_id, output_z=True)
            save(z_map, z_map_path)
        del fmri_glm

plt.show()
