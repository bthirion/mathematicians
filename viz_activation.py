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

subjects = ['cf120444','jl120341','lr120300','aa130114','aa130169','mk130199',
            'jl130200','mp130263','rm130241','al130244','bm120103','ce130459',
            'of140017','jf140025','cr140040','fm120345','hr120357','kg120369',
            'mr120371','jc130030','ld130145','cf140022','jn140034','mv140024',
            'tj140029','ap140030','af140169','pp140165','eb140248','gq140243'][:1]            


work_dir = '/neurospin/unicog/protocols/IRMf/mathematicians_Amalric_Dehaene2012/Surface_analysis/mathematicians'

subject = subjects[0]

os.environ['SUBJECTS_DIR'] = ""

subject_dir = os.path.join(work_dir, subject)
t1_dir = os.path.join(subject_dir, 't1')
surf_dir = os.path.join(t1_dir, subject, 'surf')
fun_dir = os.path.join(subject_dir, 'fmri/results')
#contrast = 'math - nonmath' # 'visual' # 'true-false' #'audio' # 'motor' # 'reflection' #
contrast = 'faces-others'

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
norig = commands.getoutput('$FREESURFER_HOME/bin/mri_info --vox2ras %s' %ref)
print norig
norig = np.array([x.split() for x in norig.split('\n')]).astype(np.float)
torig = commands.getoutput('$FREESURFER_HOME/bin/mri_info --vox2ras-tkr %s' %ref)
print torig
torig = np.array([x.split() for x in torig.split('\n')]).astype(np.float)
aff = np.linalg.inv(np.dot(norig, np.linalg.inv(torig)))


# Plot the MRI image on top of the meshes
#ref = os.path.join(t1_dir, subject, 'mri/orig.mgz')
#norig = commands.getoutput('mri_info --vox2ras %s' %ref)
#norig = np.array([x.split() for x in norig.split('\n')]).astype(np.float)
#torig = commands.getoutput('mri_info --vox2ras-tkr %s' %ref)
#torig = np.array([x.split() for x in torig.split('\n')]).astype(np.float)
#aff = np.linalg.inv(np.dot(norig, np.linalg.inv(torig)))

vol_file = os.path.join(fun_dir, '%s_z_map.nii' % contrast)
data_ = load(vol_file).get_data()
fun_affine = load(vol_file).get_affine()
aff_ = np.dot(aff, fun_affine)
viz3d.plot_map_3d(data_, aff_, cmap=cm.cold_hot,
                  vmin=-8, vmax=8, anat=False, threshold=THRESHOLD)

# mlab.show()

