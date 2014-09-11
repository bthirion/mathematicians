""" 
Launch freesurfer segmentation on the subjects of the dataset

Author: Bertrand Thirion, 2013

export SUBJECTS_DIR=''
"""
import os
import glob
from nipype.caching import Memory
from nipype.interfaces.freesurfer import ReconAll

source_dir = '/neurospin/unicog/protocols/IRMf/mathematicians_Amalric_Dehaene2012/fMRI_data/'
write_dir = '/neurospin/unicog/protocols/IRMf/mathematicians_Amalric_Dehaene2012/Surface_analysis/mathematicians'
#subjects =  ['aa130114', 'jl120341', 'mp130263', 'aa130169', 'jl130200',
#             'mr120371', 'al130244', 'kg120369', 'nl120167', 'bm120103',
#             'ld130145',  'cb120288', 'ce130459', 'rm130241', 'cf120444', 
#             'll130242', 'el120268', 
#             'lr120300', 'tb120212', 'fm120345', 'vb120303', 'hr120357', 
#             'mh120250', 'vb120409', 'jc130030', 'mk130199']

#subjects = ['cf140022', 'jf140025', 'of140017', 'tj140029', 'jn140034', 'cr140040', 'mv140024', 'ap140030']

subjects = ['af140169', 'eb140248','gq140243', 'pp140165']


for subject in subjects:
    # create directories in output_dir
    subject_dir = os.path.join(write_dir, subject)
    t1_dir = os.path.join(subject_dir, 't1')

    anat_image = glob.glob(os.path.join(t1_dir, 'anat*.nii'))[0]
    mem = Memory(base_dir=subject_dir)
    reconall =  mem.cache(ReconAll)
    recon_result = reconall(subject_id = subject, 
                            directive='all', 
                            subjects_dir = t1_dir,
                            T1_files = anat_image)
    
