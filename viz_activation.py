"""Visualize the mesh using Mayavi

Author: Bertrand Thirion, 2013

export SUBJECTS_DIR=''
"""
import os
import nibabel.freesurfer as nf
import mayavi.mlab as mlab
import numpy as np
from nibabel.gifti import read
from nibabel import load
from nipy.labs import viz3d
from nipy.labs.viz import cm
import commands

subjects = ['aa130114', 'jl120341', 'mp130263', 'aa130169', 'jl130200',
            'mr120371', 'al130244', 'kg120369', 'nl120167', 'bm120103',
            'ld130145',  'cb120288', 'ce130459', 'rm130241', 'cf120444', 
            'll130242', 'el120268', 
            'lr120300', 'tb120212', 'fm120345', 'vb120303', 'hr120357', 
            'mh120250', 'vb120409', 'jc130030', 'mk130199']

work_dir = '/neurospin/tmp/mathematicians'

subject = subjects[0]

subject_dir = os.path.join(work_dir, subject)
t1_dir = os.path.join(subject_dir, 't1')
surf_dir = os.path.join(t1_dir, subject, 'surf')
fun_dir = os.path.join(subject_dir, 'fmri/results')
contrast = 'right-left' # 'visual' # 'true-false' #'audio' # 'motor' # 'reflection' #

THRESHOLD = 3.

mlab.figure(bgcolor=(1, 1, 1))
for hemisphere in ['l', 'r']:
    mesh_file = os.path.join(surf_dir, '%sh.white' % hemisphere)
    curv_file = os.path.join(surf_dir, '%sh.curv' % hemisphere)
    fun_file = os.path.join(fun_dir, '%s_z_map_%sh.gii' %
                            (contrast, hemisphere))

    coords, triangles = nf.read_geometry(mesh_file)
    x, y, z = coords.T
    curv =  nf.read_morph_data(curv_file).astype(np.float)
    tex = np.array([darrays.data for darrays in read(fun_file).darrays]).ravel()

    name = ''
    cmin = -1
    cmax = 1
    mlab.triangular_mesh(x, y, z, triangles, transparent=True, opacity=1.,
                         name=name, scalars=curv, colormap="bone",
                         vmin=cmin, vmax=cmax)
    func_mesh = mlab.pipeline.triangular_mesh_source(
        x, y, z, triangles, scalars=tex)
    thresh = mlab.pipeline.threshold(func_mesh, low=3.)
    mlab.pipeline.surface(thresh, colormap="hot", vmin=THRESHOLD, vmax=8)


# Plot the MRI image on top of the meshes
ref = os.path.join(t1_dir, subject, 'mri/orig.mgz')
norig = commands.getoutput('mri_info --vox2ras %s' %ref)
norig = np.array([x.split() for x in norig.split('\n')]).astype(np.float)
torig = commands.getoutput('mri_info --vox2ras-tkr %s' %ref)
torig = np.array([x.split() for x in torig.split('\n')]).astype(np.float)
aff = np.linalg.inv(np.dot(norig, np.linalg.inv(torig)))

vol_file = os.path.join(fun_dir, '%s_z_map.nii' % contrast)
data_ = load(vol_file).get_data()
fun_affine = load(vol_file).get_affine()
aff_ = np.dot(aff, fun_affine)
viz3d.plot_map_3d(data_, aff_, cmap=cm.cold_hot,
                  vmin=-8, vmax=8, anat=False, threshold=THRESHOLD)

# mlab.show()

