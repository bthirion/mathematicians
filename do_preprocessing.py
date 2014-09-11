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

#subjects = ['aa130114', 'jl120341', 'mp130263', 'aa130169', 'jl130200',
#            'mr120371', 'al130244', 'kg120369', 'nl120167', 'bm120103',
#            'ld130145',  'cb120288', 'ce130459', 'rm130241', 'cf120444', 
#            'll130242', 'el120268', 
#            'lr120300', 'tb120212', 'fm120345', 'vb120303', 'hr120357', 
#            'mh120250', 'vb120409', 'jc130030', 'mk130199'][16:]
#dont = ['nl120167', 'll130242', 'tb120212']

#subjects = ['fm120345', 'hr120357', 'jc130030',  'mk130199']   

#subjects = ['cf140022', 'jf140025', 'of140017', 'tj140029', 'jn140034', 'cr140040', 'mv140024', 'ap140030'][6:]

subjects = ['cf120444','jl120341','lr120300','aa130114','aa130169','mk130199',
            'jl130200','mp130263','rm130241','al130244','bm120103','ce130459',
            'of140017','jf140025','cr140040','fm120345','hr120357','kg120369',
            'mr120371','jc130030','ld130145','cf140022','jn140034','mv140024',
            'tj140029','ap140030','af140169','pp140165','eb140248','gq140243'][13:14]

#subjects = [subject for subject in subjects if subject not in dont]

work_dir = '/neurospin/unicog/protocols/IRMf/mathematicians_Amalric_Dehaene2012/Surface_analysis/mathematicians'


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

#        do_subject_preproc(preproc_dict, concat=False, do_coreg=True,
#                           do_stc=False, do_cv_tc=True, do_realign=True,
#                           do_report=True)
        do_subject_preproc(preproc_dict, 
                             coregister=True,
                             cv_tc=0, 
                             realign=True,   
                             report=True,) 