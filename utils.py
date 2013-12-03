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


def audiosentence_contrasts(names_, ratings, n_session):
    """Create the contrasts,
    given the names of the columns of the design matrix

    names_: list of strings
    ratings: list of arrays,
             high-level variables describing the trials 
    """
    vrai, faux, meaningless = ratings
    n_reg = len(vrai) / 6  # assumes that  there are 6 eqaully length sessions

    # keep only the values relavant for the session that you are considering
    session_indexes = np.arange(n_session * n_reg, (n_session + 1)* n_reg)
    vrai, faux, meaningless = (
        vrai[session_indexes], faux[session_indexes],
        meaningless[session_indexes])
    # Caveat: assumes that  names_ order is:
    # 'alert signal', 'key press', 'response signal',  'trial_01_ ...
    reordering = np.concatenate((np.arange(30, 33), np.arange(30), 
                                np.arange(33, 39)))

    vrai, faux, meaningless = (vrai[reordering], faux[reordering], 
                               meaningless[reordering])
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
    
    contrasts['true-false'] = np.zeros(len(names_))
    contrasts['true-false'][:n_reg] = vrai * faux.sum() - faux * vrai.sum()
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


def formatted_matrix(boolean, understood_level, intuition_level,
                     immediacy_level, visual_imagery_level, mode):
    """ 
    Attempt to translate the formatted_matrix_v2 matlab function
    
    Caveat: may be buggy. Check it

    mode define what part of stimulus is analyzed: (mode = 1 pour ecoute et
    mode = 2 pour reflexion)
    """
    n_sentences = 90
    n_conds = 15
    n_sessions = 6
    n_reg = 39
    n_total = n_reg * n_sessions

    cond = np.zeros(n_total)
    positive_intuitive_cond = np.zeros(n_total)
    negative_intuitive_cond = np.zeros(n_total)
    positive_visu_cond = np.zeros(n_total)
    negative_visu_cond = np.zeros(n_total)
    
    # load data corresponding to the boolean condition "boolean"
    understood_level_cond = understood_level[boolean]
    intuition_level_cond = intuition_level[boolean]
    immediacy_level_cond = immediacy_level[boolean]
    visual_imagery_level_cond = visual_imagery_level[boolean]
    
    # indexes de 1 to 15 les indices des phrases correspondant au boolean 
    indices_for_boolean_condition = np.where(boolean)[0]
    indices_for_boolean_condition_understood = indices_for_boolean_condition[
        understood_level_cond != 0]
    
    # remove not understood sentences
    intuition_understood = intuition_level_cond[understood_level_cond != 0]
    visual_imagery_understood = visual_imagery_level_cond[
        understood_level_cond != 0]
    immediacy_understood = immediacy_level_cond[understood_level_cond != 0]
    
    # also remove all immediate responses from the intuition analysis:
    intuition_understood_non_immediate = intuition_understood[
        immediacy_understood != 1]
    visual_imagery_understood_immediate = visual_imagery_understood[
        immediacy_understood != 1]

    #compute the mean and std of intuition and visual_imagery only on the
    # sentences of interest:
    positive_mean_int = np.mean(intuition_understood_non_immediate )
    negative_mean_int = np.mean(7 - intuition_understood_non_immediate)
    positive_mean_visu = np.mean(visual_imagery_understood_immediate)
    negative_mean_visu = np.mean(7 - visual_imagery_understood_immediate)

    for i_sentence in range(n_sentences):
        if mode == 1:
            i_regressor = 2 * i_sentence + 9 * np.floor((i_sentence - 1) / n_conds)
        elif mode == 2:
            i_regressor = 2 * i_sentence + 1 + (
                9 * np.floor((i_sentence - 1) / n_conds))
        if i_sentence in indices_for_boolean_condition_understood:
            cond[i_regressor] = 1
            if i_sentence in indices_for_boolean_condition_understood[immediacy_understood != 1]:
                positive_intuitive_cond[i_regressor] =\
                    intuition_level[i_sentence] - positive_mean_int
                negative_intuitive_cond[i_regressor] =\
                    7 - intuition_level[i_sentence] - negative_mean_int
                positive_visu_cond[i_regressor] =\
                    visual_imagery_level[i_sentence] - positive_mean_visu
                negative_visu_cond[i_regressor] = 7 -\
                    visual_imagery_level[i_sentence] - negative_mean_visu
    return (cond, positive_intuitive_cond, negative_intuitive_cond, 
            positive_visu_cond, negative_visu_cond)




def define_contrast_audiosentence(response_questionnaire, correspondence):
    """
    """
    # read 4 variables in reponse_questionnaire
    understood_level = response_questionnaire[:, 1]
    intuition_level = 7 - response_questionnaire[:, 4]
    immediacy_level = response_questionnaire[:, 5]
    visual_imagery_level = response_questionnaire[:, 6]
    
    # define types
    sentence_type = np.mod(correspondence[:,1], 15) 
    # number in [0;14] for the 15 different types
    value_type = np.mod(sentence_type, 3)
    # number in [0;2] for the 3 different values of truth
    
    mode = 2
    (vrai, positive_intuitive_vrai, negative_intuitive_vrai,
     positive_visu_vrai,negative_visu_vrai) = formatted_matrix(
        value_type == 1, understood_level, intuition_level, immediacy_level,
        visual_imagery_level, mode)

    (faux, positive_intuitive_faux, negative_intuitive_faux,
     positive_visu_faux, negative_visu_faux) = formatted_matrix(
        value_type == 2, understood_level, intuition_level, immediacy_level,
        visual_imagery_level, mode)

    (meaningless, positive_intuitive_meaningless,
     negative_intuitive_meaningless, positive_visu_meaningless,
     negative_visu_meaningless) = formatted_matrix(
        value_type == 0, understood_level, intuition_level, immediacy_level,
        visual_imagery_level, mode)

    # FXME: should be: 
    vrai = np.array([
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
        1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    faux = np.array([
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
        1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    meaningless = np.array([
        0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    
    return vrai, faux, meaningless