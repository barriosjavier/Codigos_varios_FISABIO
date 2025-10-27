import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import pathlib as path
import dicom_parser
import os
import pydicom
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('path', type=str, help='Route to nifti or dicom dataset')
arg = parser.parse_args()
path_to_dataset=arg.path

def check_orientation_niftii(path_to_dataset):
    """
    Check the NIfTI orientation, and flip if needed.
    """
    ct_image= nib.nifti1.load(path_to_dataset)
    ct_arr = ct_image.get_fdata()
    x, y, z = nib.aff2axcodes(ct_image.affine)
    
    if x != 'L':
        ct_arr = np.flip(ct_arr, axis=0)
    if y != 'P':
        ct_arr = np.flip(ct_arr, axis=1)
    if z != 'S':
        ct_arr = np.flip(ct_arr, axis=2)
    
    new_nifti = nib.Nifti1Image(ct_arr.astype(float), ct_image.affine)
    nib.save(new_nifti,f'new_nifti.nii')
    



def check_orientation_dicom(path_to_dicom_series):
    """
    If Anatomical Orientation Type (0010,2210) is absent or has a value of BIPED
    X-axis is increasing to the left hand side of the patient
    Y-axis is increasing to the posterior side of the patient
    Z-axis is increasing toward the head of the patient.
    """

    PathDicom = path_to_dicom_series
    lstFilesDCM = []  # create an empty list
    for dirName, subdirList, fileList in os.walk(PathDicom):
        for filename in fileList:
            if ".dcm" in filename.lower():  # check whether the file's DICOM
                lstFilesDCM.append(os.path.join(dirName,filename))

    # Get ref file
    RefDs = pydicom.dcmread(lstFilesDCM[0])
    # Get orientation of the image 
    orientation=RefDs.data_element("ImageOrientationPatient")
    ref=RefDs.data_element("ImagePositionPatient")
    print(ref)
    # Load dimensions based on the number of rows, columns, and slices (along the Z axis)
    ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns), len(lstFilesDCM))
    # Load spacing values (in mm)
    ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(RefDs.SliceThickness))

    # The array is sized based on 'ConstPixelDims'
    ArrayDicom = np.zeros(ConstPixelDims, dtype=np.uint16)

    # loop through all the DICOM files
    for filenameDCM in lstFilesDCM:
        # read the file
        ds = pydicom.dcmread(filenameDCM)
        # store the raw image data
        ArrayDicom[:, :, lstFilesDCM.index(filenameDCM)] = ds.pixel_array

    #girar eje X
    Dicom_flipped = np.flip(ArrayDicom, axis=0)

    for i in range(Dicom_flipped.shape[2]):

        ds.PixelData=Dicom_flipped.tobytes()
        filename = os.path.join("out",f"slice_{i:03d}.dcm")
        ds.save_as(filename)
    
    return ds


extension=path.Path(path_to_dataset).suffix

if ".nii" in extension:
    check_orientation_niftii(path_to_dataset)
else:
    check_orientation_dicom(path_to_dataset)
    

#fig, axs= plt.subplots(1, 2, figsize=[10, 10])
#axs.flat[0].imshow(ct_arr[90], cmap='gray')
#axs.flat[1].imshow(ct_arr_flipped[90], cmap='gray')
#plt.tight_layout()
#plt.show()