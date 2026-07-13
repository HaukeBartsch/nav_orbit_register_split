#!/usr/bin/env python3
import sys
import os
import argparse
import json
import SimpleITK as sitk
import numpy as np

def create_mask_with_intensity_constraint(image, min_intensity=600, max_intensity=None):
    """
    Create a binary mask based on intensity constraints.
    """
    max_value = None
    if max_intensity is None:
        stats = sitk.StatisticsImageFilter()
        stats.Execute(image)
        max_value = int(stats.GetMaximum())
    else:
        max_value = int(max_intensity)
    # could we just use an upperThreshold of None? to speed it up?
    if min_intensity > max_value:
        raise ValueError(f"min_intensity ({min_intensity}) is greater than the maximum intensity in the image ({max_value}).")

    return sitk.BinaryThreshold(image, 
                                lowerThreshold=min_intensity,
                                upperThreshold=max_value,
                                insideValue=1,
                                outsideValue=0
    )

def dilateIt(image, filter_size=3):
    filter = sitk.BinaryDilateImageFilter()
    filter.SetKernelRadius([filter_size,filter_size,filter_size])
    filter.SetForegroundValue(1)
    filter.SetBackgroundValue(0)

    return filter.Execute(image)


def elastix_3d_rigid_to_matrix(params):
    if len(params) != 6:
        raise ValueError("Elastix 3D rigid transform requires exactly 6 parameters.")
    
    tx, ty, tz, rx, ry, rz = params
    
    # Rotation matrices (column-vector convention)
    def rot_x(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    
    def rot_y(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    
    def rot_z(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    
    # ITK/Elastix applies rotations in Z-Y-X order: R = Rz * Ry * Rx
    R = rot_z(rz) @ rot_y(ry) @ rot_x(rx)
    
    # Homogeneous 4x4 matrix (row-major for Python/NumPy)
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = [tx, ty, tz]
    
    return T

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='Analyze two ct images with masks')
    
    # Add named arguments
    parser.add_argument('--inputCT1', '-i1', required=True, help='CT image fixed')
    parser.add_argument('--inputCT2', '-i2', required=True, help='CT image moving')
    parser.add_argument('--inputMask1', '-m1', required=True, help='Mask image fixed')
    parser.add_argument('--inputMask2', '-m2', required=True, help='Mask image moving')

    parser.add_argument('--output', '-o', required=True, help='Output folder')

    # Parse arguments
    args = parser.parse_args()
    
    # Access the named arguments
    ct1_name = args.inputCT1
    ct2_name = args.inputCT2
    mask1_name = args.inputMask1
    mask2_name = args.inputMask2
    output = args.output

    fixed_image_name = None
    moving_image_name = None

    if not os.path.exists(output):
        os.makedirs(output)

    # if the input is a folder we should assume we read DICOM data (all files in all directories)
    ct1 = None
    if os.path.isdir(ct1_name):
        # read as dicom
        reader = sitk.ImageSeriesReader()
        reader.MetaDataDictionaryArrayUpdateOn();
        dicom_names = reader.GetGDCMSeriesFileNames(ct1_name, recursive=True)
        # print(dicom_names)
        # Read the series
        reader.SetFileNames(dicom_names)
        ct1 = reader.Execute()
        fixed_image_name = reader.GetMetaData(0,"0010|0020") 
    else:
        ct1 = sitk.ReadImage(ct1_name)

    ct2 = None
    if os.path.isdir(ct2_name):
        # read as dicom
        reader = sitk.ImageSeriesReader()
        reader.MetaDataDictionaryArrayUpdateOn();
        dicom_names = reader.GetGDCMSeriesFileNames(ct2_name, recursive=True)

        # Read the series
        reader.SetFileNames(dicom_names)
        ct2 = reader.Execute()
        moving_image_name = reader.GetMetaData(0,"0010|0020")
    else:
        ct2 = sitk.ReadImage(ct2_name)

    mask1 = None
    if os.path.isdir(mask1_name):
        # read as dicom
        reader = sitk.ImageSeriesReader()
        dicom_names = reader.GetGDCMSeriesFileNames(mask1_name, recursive=True)

        # Read the series
        reader.SetFileNames(dicom_names)
        mask1 = reader.Execute()
    else:
        mask1 = sitk.ReadImage(mask1_name)

    mask2 = None
    if os.path.isdir(mask2_name):
        # read as dicom
        reader = sitk.ImageSeriesReader()
        dicom_names = reader.GetGDCMSeriesFileNames(mask2_name, recursive=True)

        # Read the series
        reader.SetFileNames(dicom_names)
        mask2 = reader.Execute()
    else:
        mask2 = sitk.ReadImage(mask2_name)



    def _compare_geometry(name_a, img_a, name_b, img_b) -> bool:
        """Check that two images share the same physical geometry; warn if not."""
        mismatches = []
        if img_a.GetSize() != img_b.GetSize():
            mismatches.append(
                f"sizes differ: {name_a}={img_a.GetSize()}, {name_b}={img_b.GetSize()}"
            )
        if img_a.GetOrigin() != img_b.GetOrigin():
            mismatches.append(
                f"origins differ: {name_a}={img_a.GetOrigin()}, {name_b}={img_b.GetOrigin()}"
            )
        if img_a.GetSpacing() != img_b.GetSpacing():
            mismatches.append(
                f"spacings differ: {name_a}={img_a.GetSpacing()}, {name_b}={img_b.GetSpacing()}"
            )
        if img_a.GetDirection() != img_b.GetDirection():
            mismatches.append(
                f"direction cosines differ: {name_a}={img_a.GetDirection()}, {name_b}={img_b.GetDirection()}"
            )
        if mismatches:
            print(f"WARNING: {name_a} and {name_b} do not share the same geometry:")
            for m in mismatches:
                print(f"  - {m}")
            print(
                f"  The mask-based registration will produce a transform in a "
                f"different coordinate system than the CT registration. "
                f"Consider resampling the masks to the CT geometry first."
            )
            return False
        else:
            print(f"OK: {name_a} and {name_b} share the same geometry.")
        return True

    if not(_compare_geometry("mask1", mask1, "ct1", ct1)):
        exit(-1)
    if not(_compare_geometry("mask2", mask2, "ct2", ct2)):
        exit(-1)

    # To make the registration better we should estimate a translation of the
    # mask2 onto mask1 first.  That transformation would overlay the eyeballs
    # sufficiently for a later stage registration using elastix.
    maskElastixImageFilter = sitk.ElastixImageFilter()
    output_dir = os.path.join(output, "mask_registration_output")
    os.makedirs(output_dir, exist_ok=True)
    maskElastixImageFilter.SetOutputDirectory(output_dir)
    maskElastixImageFilter.LogToFileOn()

    parameterMapVector = []
    parameterMap = sitk.GetDefaultParameterMap('rigid')
    parameterMap["Transform"] = ("EulerTransform", )
    parameterMap["UseDirectionCosines"] = ("true", )
    #parameterMap["Metric"] = ("AdvancedMeanSquares", )
    parameterMap["DefaultPixelValue"] = ("0",)
    parameterMap["NumberOfSpatialSamples"] = ("2000",)
    parameterMap["MaximumNumberOfIterations"] = ("500",)
    parameterMap["AutomaticTransformInitialization"] = ("true",)
    parameterMap["AutomaticTransformInitializationMethod"] = ("GeometricalCenter",)
    parameterMapVector.append(parameterMap)

    maskElastixImageFilter.SetParameterMap(parameterMapVector)

    maskElastixImageFilter.SetFixedImage(mask1)
    maskElastixImageFilter.SetMovingImage(mask2)
    maskElastixImageFilter.Execute()

    mask2OnMask1Image = maskElastixImageFilter.GetResultImage()
    sitk.WriteImage(mask2OnMask1Image, output + "/mask2OnMask1.nii.gz")
    sitk.WriteImage(mask1, output + "/mask1.nii.gz")

    print("We have found a mask registration")

    # the resulting transformation is this:
    transformParameterMapVector = maskElastixImageFilter.GetTransformParameterMap()
    if len(transformParameterMapVector) > 1:
        print("Warning, we only use the first transformation here")
    maskTransformParameterMap = maskElastixImageFilter.GetTransformParameterMap(0)
    mask_tx_file = os.path.join(output, "mask_initial_transform.txt")
    sitk.WriteParameterFile(maskTransformParameterMap, mask_tx_file)


    elastixImageFilter = sitk.ElastixImageFilter()
    output_dir = "/tmp/registration_output"
    os.makedirs(output_dir, exist_ok=True)
    elastixImageFilter.SetOutputDirectory(output_dir)
    elastixImageFilter.LogToFileOn()

    # now register the CT volumes using that as an initial registration
    parameterMapVector = []
    # copy all transformations over from the transformParameterMapVector
    #for trans in transformParameterMapVector:
    #    parameterMapVector.append(trans)

    #parameterMapVector.append(maskTransformParameterMap)

    parameterMap = sitk.GetDefaultParameterMap('rigid')
    #parameterMap["InitialTransformParameterFileName"] = (mask_tx_file,)
    parameterMap["UseDirectionCosines"] = ("true", )
    #parameterMap["Transform"] = ("EulerTransform", )
    parameterMap["Metric"] = ("AdvancedMattesMutualInformation", )
    #parameterMap["Optimizer"] = ("AdaptiveStochasticGradientDescent", )
    #parameterMap["Metric"] = ("AdvancedMeanSquares", )
    parameterMap["Optimizer"] = ("AdaptiveStochasticGradientDescent", )
    parameterMap["ASGDParameterEstimationMethod"] = ("DisplacementDistribution", )
    parameterMap["Registration"] = ("MultiResolutionRegistration", )

    parameterMap["ImageSampler"] = ("RandomCoordinate", )
    parameterMap["NewSamplesEveryIteration"] = ("true", )

    parameterMap["FixedImagePyramid"] = ("FixedSmoothingImagePyramid", )
    parameterMap["MovingImagePyramid"] = ("MovingSmoothingImagePyramid", )
    #parameterMap["FixedImagePyramid"] = ("FixedRecursiveImagePyramid", )
    #parameterMap["MovingImagePyramid"] = ("MovingRecursiveImagePyramid", )
    #(ImagePyramidSchedule 8 8 8  4 4 4  2 2 2  1 1 1)

    parameterMap["NumberOfResolutions"] = ("7", )
    #parameterMap["ImagePyramidSchedule"] = ("15", "15", "15", "4", "4", "4", "2", "2", "2", "1", "1", "1", )
    #parameterMap["NumberOfResolutions"] =  ("6", )
    #parameterMap["ImagePyramidSchedule"] = ("11", "11", "4", "10", "10", "4", "8", "8", "2", "4", "4", "1", "2", "2", "1", "1", "1", "1", )
    #parameterMap["MaximumNumberOfIterations"] =  ("10000", "10000", "10000", "500", )
    #parameterMap["MinimumStepLength"] = ("0.1", "0.05", "0.01", "0.005", "0.0001", )

    #parameterMap["NumberOfHistogramBins"] = ("128", "64", "64", )
    parameterMap["DefaultPixelValue"] = ("-1024",)
    parameterMap["NumberOfSpatialSamples"] = ("2000",)
    #parameterMap["NewSamplesEveryIteration"] = ("true", )
    #parameterMap["MaximumNumberOfIterations"] = ("600",)
    parameterMap["NumberOfHistogramBins"] = ("64", )
    parameterMap["ErodeFixedMask"] = ("false",)
    parameterMap["ErodeMovingMask"] = ("false",)
    parameterMap["MaximumNumberOfSamplingAttempts"] = ("8",)
    parameterMap["AutomaticScalesEstimation"] = ("true", )
    #parameterMap["Scales"] = ["1.0", "1.0", "1.0", "10000.0", "10000.0", "10000.0" ]
    parameterMap["AutomaticTransformInitialization"] = ("false",)
    parameterMap["ResampleInterpolator"] = ("FinalBSplineInterpolator", )
    parameterMap["HowToCombineTransforms"] = ("Compose", )
    #parameterMap["AutomaticTransformInitializationMethod"] = ("GeometricalCenter",)
    parameterMapVector.append(parameterMap)

    #parameterMapVector = [
    #    #sitk.GetDefaultParameterMap("translation"),
    #    parameterMap
    #]
    elastixImageFilter.SetParameterMap(parameterMapVector)
    elastixImageFilter.SetInitialTransformParameterFileName(mask_tx_file)

    elastixImageFilter.SetFixedImage(ct1)
    elastixImageFilter.SetMovingImage(ct2)

    # The essential problem here was that the intensity range of CT can be
    # strongly influenced by metal artifacts (streaks). We need to use a mask that
    # removes the metal artifacts especially if they only appear in one of the two volumes.
    # Operators might angle the CT scans so that teeth are not part of the scan. If the
    # operator does this only in one of the scans we end up with miss-alignment if we 
    # do not remove the streaks using our mask. Its sufficient to do this on the fixed image.
    bone_mask = create_mask_with_intensity_constraint(ct1, 200, 600)
    #all_mask = create_mask_with_intensity_constraint(ct1, 200)
    sitk.WriteImage(dilateIt(bone_mask, 4), output + "/bone_mask_200_fixed.nii.gz")
    elastixImageFilter.SetFixedMask(dilateIt(bone_mask, 4))

    # use a region around the eyeballs as the mask
    #fixedMask = dilateIt(mask1, 7)
    # convert this mask to sitkUInt8 using Cast
    #fixedMask = sitk.Cast(fixedMask, sitk.sitkUInt8)
    #elastixImageFilter.SetFixedMask(fixedMask)
    #elastixImageFilter.SetFixedMask(bone_mask)
    
    bone_mask2 = create_mask_with_intensity_constraint(ct2, 200, 600)
    #sitk.WriteImage(bone_mask2, output + "/bone_mask_200_moving.nii.gz")
    #elastixImageFilter.SetMovingMask(dilateIt(bone_mask2, 4))
    #movingMask = dilateIt(mask2, 4)
    # convert this mask to sitkUInt8 using Cast
    #movingMask = sitk.Cast(movingMask, sitk.sitkUInt8)
    #elastixImageFilter.SetMovingMask(movingMask)

    elastixImageFilter.Execute()

    # Get the mutual information metric value from the log file
    final_metric = None
    log_path = os.path.join(output_dir, "elastix.log")
    if os.path.exists(log_path):
        with open(log_path, 'r') as f:
            for line in f:
                # Elastix records: "Final metric value = -0.592680"
                if "Final metric value" in line:
                    final_metric = float(line.split("=")[-1].strip())
                    break

    #print("Final Metric Value:", final_metric)

    resultImage = elastixImageFilter.GetResultImage()
    transformParameterMap = elastixImageFilter.GetTransformParameterMap(0)
    sitk.WriteImage(resultImage, output + "/ct_moved_resampled.nii.gz")
    sitk.WriteImage(ct1, output + "/ct_fixed.nii.gz")

    # Extract the scale factors from the affine transform (params 3-5 = scale X, Y, Z)
    transform_params = transformParameterMap["TransformParameters"]
    parameters = np.array(transformParameterMap["TransformParameters"], dtype=float)
    A = None
    if len(parameters) == 12:
        A = np.array(parameters[:9]).reshape(3, 3)
    elif len(parameters) == 6:
        A =  np.eye(3) # transformParameterMap.GetMatrix()

    affine_scale = { "x": 0, "y": 0, "z": 0 }
    if transformParameterMap["Transform"][0] == "AffineTransform":        
        # The physical scale per axis is the L2 norm of each column/row (depending on ITK convention)
        affine_scale = {
            "x": round(float(A[0][0]), 4),
            "y": round(float(A[1][1]), 4),
            "z": round(float(A[2][2]), 4),
        }
    #affine_scale = {
    #    "x": round(float(transform_params[3]), 4),
    #    "y": round(float(transform_params[4]), 4),
    #    "z": round(float(transform_params[5]), 4),
    #}

    transformixImageFilter = sitk.TransformixImageFilter()
    transformParameterMap["ResampleInterpolator"] = ("FinalNearestNeighborInterpolator",)
    transformParameterMap["DefaultPixelValue"] = ("0",)
    transformixImageFilter.SetTransformParameterMap(transformParameterMap)

    transformixImageFilter.SetMovingImage(mask2)
    transformixImageFilter.Execute()
    sitk.WriteImage(transformixImageFilter.GetResultImage(), output + "/mask_moved_resampled.nii.gz")
    mask2_registered = transformixImageFilter.GetResultImage()

    # Compute a dice coefficient between the two masks for quality control
    mask1_array = sitk.GetArrayFromImage(mask1)
    mask2_registered_array = sitk.GetArrayFromImage(mask2_registered)

    # Compute Dice between fixed mask and registered moving mask
    fixed_binary = (mask1_array > 0).astype(float)
    moved_binary = (mask2_registered_array > 0).astype(float)
    intersection = np.sum(fixed_binary * moved_binary)
    full_masks_dice = 2.0 * intersection / (np.sum(fixed_binary) + np.sum(moved_binary))
    # print(f"Mask Dice coefficient: {full_masks_dice:.4f}")

    # we should split the transformed masks now, both on the mask1 and on mask2
    # get center of mass for the mask

    mask_array = sitk.GetArrayFromImage(mask2_registered)
        
    # Calculate center of mass
    coords = np.where(mask_array > 0)
        
    if len(coords[0]) == 0:
        print("Error: No mask pixel found in the resampled mask, registration error!")
        exit(-1)
            
    # Calculate center of mass for each dimension
    center_x = np.mean(coords[1])  # x coordinate (column index)
    center_y = np.mean(coords[0])  # y coordinate (row index)
    center_z = np.mean(coords[2])  # z coordinate (image index)

    min_x = np.min(coords[1])
    min_y = np.min(coords[0])
    min_z = np.min(coords[2])

    max_x = np.max(coords[1])
    max_y = np.max(coords[0])
    max_z = np.max(coords[2])

    # one of these should be the sagittal direction, cut along it a couple of times
    # cut into 4 pieces
    x_coords = coords[0]
    y_coords = coords[1]
    z_coords = coords[2]

    #
    # Split eyeballs
    #
    condition_mask = z_coords > center_z
    mask_array[x_coords[condition_mask], y_coords[condition_mask], z_coords[condition_mask]] = 2

    # we need to split again, with a given distance away from the center
    # we could use the last slice in the axial direction, split at its center
    coords_half1 = np.where(mask_array == 1)
    coords_half2 = np.where(mask_array == 2)
    # max in y direction
    x_coords = coords_half1[0]
    y_coords = coords_half1[1]
    z_coords = coords_half1[2]
    max_x = np.max(x_coords)
    condition_mask = x_coords == max_x

    # the last slice
    # mask_array[x_coords[condition_mask], y_coords[condition_mask], z_coords[condition_mask]] = 3

    #
    # now split laterally again
    #
    split1 = np.mean(z_coords[condition_mask])
    condition_mask = z_coords < split1
    mask_array[x_coords[condition_mask], y_coords[condition_mask], z_coords[condition_mask]] = 3

    #
    # do the same with the other eye 
    #
    x_coords = coords_half2[0]
    y_coords = coords_half2[1]
    z_coords = coords_half2[2]
    max_x = np.max(x_coords)
    condition_mask = x_coords == max_x

    # the last slice
    # mask_array[x_coords[condition_mask], y_coords[condition_mask], z_coords[condition_mask]] = 5
    #
    # now split laterally again
    #
    split2 = np.mean(z_coords[condition_mask])
    condition_mask = z_coords > split2
    mask_array[x_coords[condition_mask], y_coords[condition_mask], z_coords[condition_mask]] = 4

    # now save that output mask for the moved image
    # We should do the same splits with the fixed mask...
    #
    #
    split_mask = sitk.GetImageFromArray(mask_array)
    split_mask.CopyInformation(mask2_registered)
    volumes_per_region = {
        "final_mutual_information": -final_metric,
        "full_masks_dice": full_masks_dice,
        "affine_scale": affine_scale,
        "volume_change_ratio": {},
        "moved": {},
        "fixed": {}
    }
    # compute stats
    # 1 - Left eye / medial
    # 2 - Right eye / medial
    # 3 - Left eye / lateral
    # 4 - Right eye / lateral
    names_dict = [ 
        "background", 
        "Left eye / medial", 
        "Right eye / medial", 
        "Left eye / lateral", 
        "Right eye / lateral"
    ]
    shape_stats = sitk.LabelShapeStatisticsImageFilter()
    # convert the split_mask into uint8
    m = sitk.Cast(split_mask, sitk.sitkInt8)
    shape_stats.Execute(m)
    labels = shape_stats.GetLabels()
    for label in labels:
        # GetPhysicalSize returns the volume in the image's physical units
        volume_mm3 = shape_stats.GetPhysicalSize(label)
        nam = names_dict[label]
        volumes_per_region["fixed"][nam] = volume_mm3 / 1000
        print(f"{nam}: Volume = {volumes_per_region["fixed"][nam]:.2f} cm³")

    sitk.WriteImage(split_mask, output + "/split_mask_fixed_resampled.nii.gz")

    # Apply the same splits with the mask1 and save that also
    mask_array = sitk.GetArrayFromImage(mask1)
    coords = np.where(mask_array > 0)
    # split left / right
    x_coords = coords[0]
    y_coords = coords[1]
    z_coords = coords[2]

    #
    # Split eyeballs (use the existing center_z)
    #
    condition_mask = z_coords > center_z
    mask_array[x_coords[condition_mask], y_coords[condition_mask], z_coords[condition_mask]] = 2

    condition_mask = z_coords < split1
    mask_array[x_coords[condition_mask], y_coords[condition_mask], z_coords[condition_mask]] = 3

    condition_mask = z_coords > split2
    mask_array[x_coords[condition_mask], y_coords[condition_mask], z_coords[condition_mask]] = 4

    split_mask = sitk.GetImageFromArray(mask_array)
    split_mask.CopyInformation(mask1)
    # compute stats
    # 1 - Left eye / medial
    # 2 - Right eye / medial
    # 3 - Left eye / lateral
    # 4 - Right eye / lateral
    names_dict = [ 
        "background", 
        "Left eye / medial", 
        "Right eye / medial", 
        "Left eye / lateral", 
        "Right eye / lateral"
    ]
    shape_stats = sitk.LabelShapeStatisticsImageFilter()
    # convert the split_mask into uint8
    m = sitk.Cast(split_mask, sitk.sitkInt8)
    shape_stats.Execute(m)
    labels = shape_stats.GetLabels()
    for label in labels:
        # GetPhysicalSize returns the volume in the image's physical units
        volume_mm3 = shape_stats.GetPhysicalSize(label)
        nam = names_dict[label]
        volumes_per_region["moved"][nam] = volume_mm3 / 1000
        print(f"{nam}: Volume = {volumes_per_region["moved"][nam]:.2f} cm³")

    sitk.WriteImage(split_mask, output + "/split_mask_moved_resampled.nii.gz")

    # add volume change ratios
    for region in names_dict[1:]:  # skip "background"
        fv = volumes_per_region["fixed"][region]
        mv = volumes_per_region["moved"][region]
        if fv > 0:
            volumes_per_region["volume_change_ratio"][region] = round(mv / fv, 4)
        else:
            volumes_per_region["volume_change_ratio"][region] = None

    volumes_per_region["fixed_image_name"] = fixed_image_name.strip() if fixed_image_name else "unknown"
    volumes_per_region["moving_image_name"] = moving_image_name.strip() if moving_image_name else "unknown"

    # save the volume info
    with open(output + "/volumes.json", "w") as file:
        json.dump(volumes_per_region, file, indent=2)

    # create a csv version of the same file
    import pandas as pd
    last_folder_name = os.path.basename(os.path.normpath(output))

    # Flatten the nested JSON structure
    df = pd.json_normalize(volumes_per_region, sep='_')

    # Clean up column names by removing spaces and slashes
    df.columns = df.columns.str.replace(' / ', '_').str.replace(' ', '_')

    df.to_csv(output + ("/%s.csv" % (last_folder_name)), index=False)


if __name__ == "__main__":
    main()
