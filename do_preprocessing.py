""" Perform the preprcessing in a pure Python implementation

Note 
----
* requires joblib 0.6.5

Author: Bertrand thirion, Michael Eickenberg, 2013
"""
import os
import glob
from pypreprocess.purepython_preproc_utils import do_subject_preproc
import warnings

subjects = ['aa130114', 'jl120341', 'mp130263', 'aa130169', 'jl130200',
            'mr120371', 'al130244', 'kg120369', 'nl120167', 'bm120103',
            'ld130145',  'cb120288', 'ce130459', 'rm130241', 'cf120444', 
            'll130242', 'el120268', 
            'lr120300', 'tb120212', 'fm120345', 'vb120303', 'hr120357', 
            'mh120250', 'vb120409', 'jc130030', 'mk130199'][16:]
dont = ['nl120167', 'll130242', 'tb120212']
subjects = [subject for subject in subjects if subject not in dont]

work_dir = '/neurospin/tmp/mathematicians'

if __name__ == "__main__":

    for subject in subjects:
        subject_dir = os.path.join(work_dir, subject)
        t1_dir = os.path.join(subject_dir, 't1')
        fmri_dir = os.path.join(subject_dir, 'fmri')
        anat_image = glob.glob(os.path.join(t1_dir, 'anat*.nii'))[0]
        fmri_images = glob.glob(os.path.join(fmri_dir,
                                             'visualcategs/visu*.nii'))
        fmri_images += glob.glob(os.path.join(fmri_dir,
                                              'audiosentence/audio*.nii'))
        fmri_images += glob.glob(os.path.join(fmri_dir,
                                              'localizer/loc*.nii'))
        
        preproc_dict = {
            'n_sessions': len(fmri_images),
            'func': fmri_images,
            'anat': anat_image,
            'subject_id': subject,
            'output_dir': fmri_dir
        }

        do_subject_preproc(preproc_dict, concat=False, do_coreg=True,
                           do_stc=False, do_cv_tc=True, do_realign=True,
                           do_report=False)
