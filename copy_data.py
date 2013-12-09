""" 
This scripts copies data to some working directory

Author: Bertrand Thirion, 2013
"""
import os
import glob
import shutil

source_dir = '/neurospin/unicog/protocols/IRMf/mathematicians_Amalric_Dehaene2012/fMRI_data/'
write_dir = '/neurospin/tmp/mathematicians'
subjects =  ['aa130114', 'jl120341', 'mp130263', 'aa130169', 'jl130200',
             'mr120371', 'al130244', 'kg120369', 'nl120167', 'bm120103',
             'ld130145',  'cb120288', 'ce130459', 'rm130241', 'cf120444', 
             'll130242', 'el120268', 
             'lr120300', 'tb120212', 'fm120345', 'vb120303', 'hr120357', 
             'mh120250', 'vb120409', 'jc130030', 'mk130199']

if os.path.exists(write_dir) == False:
    os.mkdir(write_dir)

for subject in subjects:
    # create directories in output_dir
    subject_dir = os.path.join(write_dir, subject)
    t1_dir = os.path.join(subject_dir, 't1')
    fmri_dir = os.path.join(subject_dir, 'fmri')
    for dir_ in [subject_dir, t1_dir, fmri_dir]:
        if os.path.exists(dir_) == False:
            os.mkdir(dir_)

    # copy the t1 data
    src_ =  os.path.join(source_dir, subject)
    src_file = glob.glob(os.path.join(src_, 'anat/anat*.nii'))[0]
    shutil.copy(src_file, t1_dir)

    # copy the fMRI data
    
    # audiosentence
    src_ =  os.path.join(source_dir, subject, 'fMRI/audiosentence')
    src_files = glob.glob(os.path.join(src_, 'audio*.nii'))
    out_ = os.path.join(fmri_dir, 'audiosentence')
    if os.path.exists(out_) == False:
            os.mkdir(out_)
    for src_file in src_files:
        shutil.copy(src_file, out_)
    
    # localizer
    src_ =  os.path.join(source_dir, subject, 'fMRI/localizer')
    src_files = glob.glob(os.path.join(src_, 'localizer*.nii'))
    out_ = os.path.join(fmri_dir, 'localizer')
    if os.path.exists(out_) == False:
            os.mkdir(out_)
    for src_file in src_files:
        shutil.copy(src_file, out_)
    
    # visualcategs
    src_ =  os.path.join(source_dir, subject, 'fMRI/visualcategs')
    src_files = glob.glob(os.path.join(src_, 'visu*.nii'))
    out_ = os.path.join(fmri_dir, 'visualcategs')
    if os.path.exists(out_) == False:
            os.mkdir(out_)
    for src_file in src_files:
        shutil.copy(src_file, out_)
                          



    
