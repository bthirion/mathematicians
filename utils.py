""" Utility functions for GLM analysis

Author: Bertrand Thirion, 2013
"""

import numpy as np
from scipy.io import loadmat
from nibabel import load, Nifti1Image

from nipy.modalities.fmri.experimental_paradigm import BlockParadigm
from nipy.modalities.fmri.design_matrix import make_dmtx



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


def audiosentence_dmtx(onset_file, motion_file, n_scans, tr):
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


def fixed_effects_img(con_imgs, var_imgs, mask_img):
    """Compute the fixed effets given images of first-level effects and variance
    """
    con, var = [], []
    mask = mask_img.get_data().astype(np.bool)
    for (con_img, var_img) in zip(con_imgs, var_imgs):
        con.append(con_img.get_data()[mask])
        var.append(var_img.get_data()[mask])
    
    arrays = fixed_effects(con, var)
    """
    tiny = 1.e-16
    con, var = np.asarray(con), np.asarray(var)
    var = np.maximum(var, tiny)
    prec = 1./ var
    ffx_con = np.sum(con * prec, 0) * 1./ np.sum(prec, 0)
    ffx_var = 1./ np.sum(prec, 0)
    ffx_stat = ffx_con / np.sqrt(ffx_var)
    arrays = [ffx_con, ffx_var, ffx_stat]
    """
    outputs = []
    for array in arrays:
        vol = mask.astype(np.float)
        vol[mask] = array
        outputs.append(Nifti1Image(vol, mask_img.get_affine()))
    return outputs

def fixed_effects(contrasts, variances):
    """Compute the fixed effets given arrays of first-level effects and variance
    """
    tiny = 1.e-16
    con, var = np.asarray(contrasts), np.asarray(variances)
    var = np.maximum(var, tiny)
    prec = 1./ var
    ffx_con = np.sum(con * prec, 0) * 1./ np.sum(prec, 0)
    ffx_var = 1./ np.sum(prec, 0)
    ffx_stat = ffx_con / np.sqrt(ffx_var)
    arrays = [ffx_con, ffx_var, ffx_stat]
    return arrays

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
