import numpy as np
from nipy.labs.datasets import VolumeImg
from nibabel import load, save, Nifti1Image

def resample_anat(t1_dir, anatomical_image):
	anat_affine = load(anatomical_image).get_affine()

	template_affine = anat_affine.copy()
	# change the voxel size
	template_affine[:3, :3] *= 2
	template_shape = np.array(load(anatomical_image).get_shape())[:3] / 2


	data = load(anatomical_image).get_data().copy()
	data[np.isnan(data)] = 0
	input_image = VolumeImg(data, anat_affine, 'arbitrary')
	resampled_image = input_image.as_volume_img(template_affine,
	                                            template_shape)
	# write the output
	wim = Nifti1Image(resampled_image.get_data().astype('int16'),
	                                          template_affine)
	wim.set_data_dtype('int16')
	resampled_anat = anatomical_image + '_2mm.nii'
	save(wim, resampled_anat)

	return resampled_anat
