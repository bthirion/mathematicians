""" This scripts copies data to some working directory
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

    """
    # copy the data
    src_ =  os.path.join(source_dir, subject)
    src_file = glob.glob(os.path.join(src_, 'anat/anat*.nii'))[0]
    shutil.copy(src_file, t1_dir)
    """

    # should be put elsewhere, but launch recon_all now and here
    # export SUBJECTS_DIR=''
    anat_image = glob.glob(os.path.join(t1_dir, 'anat*.nii'))[0]
    from nipype.caching import Memory
    mem = Memory(base_dir=subject_dir)
    from nipype.interfaces.freesurfer import ReconAll
    reconall =  mem.cache(ReconAll)
    recon_result = reconall(subject_id = subject, 
                            directive='all', 
                            subjects_dir = t1_dir,
                            T1_files = anat_image)
    
