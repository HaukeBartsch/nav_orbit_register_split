#!/usr/bin/env python3
import sys
import os
import argparse
import json
import SimpleITK as sitk
import numpy as np

def create_mask_with_intensity_constraint(image, min_intensity=600):
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    max_value = int(stats.GetMaximum())
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

    if not os.path.exists(output):
        os.makedirs(output)

    # if the input is a folder we should assume we read DICOM data (all files in all directories)
    ct1 = None
    if os.path.isdir(ct1_name):
        # read as dicom
        reader = sitk.ImageSeriesReader()
        dicom_names = reader.GetGDCMSeriesFileNames(ct1_name, recursive=True)
        # print(dicom_names)
        # Read the series
        reader.SetFileNames(dicom_names)
        ct1 = reader.Execute()
    else:
        ct1 = sitk.ReadImage(ct1_name)

    ct2 = None
    if os.path.isdir(ct2_name):
        # read as dicom
        reader = sitk.ImageSeriesReader()
        dicom_names = reader.GetGDCMSeriesFileNames(ct2_name, recursive=True)

        # Read the series
        reader.SetFileNames(dicom_names)
        ct2 = reader.Execute()
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



    elastixImageFilter = sitk.ElastixImageFilter()
    output_dir = "/tmp/registration_output"
    os.makedirs(output_dir, exist_ok=True)
    elastixImageFilter.SetOutputDirectory(output_dir)
    elastixImageFilter.LogToFileOn()
    # parameterMap = sitk.GetDefaultParameterMap('translate')
    # parameterMap = sitk.GetDefaultParameterMap('affine')
    parameterMap = sitk.GetDefaultParameterMap('affine')
    parameterMap["DefaultPixelValue"] = ("-1024",)
    parameterMap["NumberOfSpatialSamples"] = ("1000",)
    parameterMap["MaximumNumberOfIterations"] = ("10000",)
    parameterMap["ErodeFixedMask"] = ("false",)
    parameterMap["MaximumNumberOfSamplingAttempts"] = ("8",)
    parameterMap["AutomaticTransformInitialization"] = ("true",)
    parameterMap["AutomaticTransformInitializationMethod"] = ("CenterOfGravity",)
    #parameterMap["ShowExactMetricValue"] = ["true"] 
    #parameterMap["InitialTransformParametersFileName"] = (initial_transform_file,)
    #parameterMap["NumberOfHistogramBins"] = ("256",)

    parameterMapVector = [
        #sitk.GetDefaultParameterMap("translation"),
        parameterMap
    ]
    elastixImageFilter.SetParameterMap(parameterMapVector)

    elastixImageFilter.SetFixedImage(ct1)
    # we can only use the mask if we dilate it as well, what about smoothing?
    bone_mask = create_mask_with_intensity_constraint(ct1, 200)
    elastixImageFilter.SetMovingImage(ct2)
    elastixImageFilter.SetFixedMask(dilateIt(bone_mask, 7))
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
        print("No coordinates found in the mask")
            
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

    # save the volume info
    with open(output + "/volumes.json", "w") as file:
        json.dump(volumes_per_region, file, indent=2)

    # create a csv version of the same file
    import pandas as pd
    data = volumes_per_region

    # Flatten the nested JSON structure
    df = pd.json_normalize(data, sep='_')

    # Clean up column names by removing spaces and slashes
    df.columns = df.columns.str.replace(' / ', '_').str.replace(' ', '_')

    df.to_csv(output + ("/%s.csv" % (output.split('/')[-1])), index=False)


if __name__ == "__main__":
    main()
