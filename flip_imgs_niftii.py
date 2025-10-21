import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt


def check_orientation(ct_image, ct_arr):
    """
    Check the NIfTI orientation, and flip to  'RPS' if needed.
    :param ct_image: NIfTI file
    :param ct_arr: array file
    :return: array after flipping
    """
    x, y, z = nib.aff2axcodes(ct_image.affine)
    if x != 'R':
        ct_arr = np.flip(ct_arr, axis=0)
    if y != 'P':
        ct_arr = np.flip(ct_arr, axis=1)
    if z != 'S':
        ct_arr = np.flip(ct_arr, axis=2)
    return ct_arr


ct_image= nib.nifti1.load('CT_Abdo')
ct_arr = ct_image.get_fdata()


ct_arr_flipped=check_orientation(ct_image,ct_arr)

fig, axs= plt.subplots(1, 2, figsize=[10, 10])

axs.flat[0].imshow(ct_arr[90], cmap='gray')
axs.flat[1].imshow(ct_arr_flipped[90], cmap='gray')

plt.tight_layout()
plt.show()
