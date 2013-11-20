"""Visualize the mesh using Mayavi

Author: Bertrand Thirion, 2013
"""
import os
import nibabel.freesurfer as nf
import mayavi.mlab as mlab
import numpy as np
from nibabel.gifti import read

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
contrast =   'motor' # 'reflection' # 'visual' # 'audio' # 

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
    mlab.triangular_mesh(x, y, z, triangles, transparent=False, opacity=1.,
                         name=name, scalars=curv, colormap="bone",
                         vmin=cmin, vmax=cmax)
    func_mesh = mlab.pipeline.triangular_mesh_source(
        x, y, z, triangles, scalars=tex)
    thresh = mlab.pipeline.threshold(func_mesh, low=3.)
    mlab.pipeline.surface(thresh, colormap="hot", vmin=3, vmax=8)


mlab.show()
