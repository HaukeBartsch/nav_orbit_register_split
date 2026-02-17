#!/usr/bin/env python3
import sys
import argparse
import json
import SimpleITK as sitk
import numpy as np

def create_mask_with_intensity_constraint2(image, min_intensity=600):
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    max_value = int(stats.GetMaximum())
    # could we just use an upperThreshold of None? to speed it up?

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

def create_mask_with_intensity_constraint(image, min_intensity=600):
    # Create the BinaryThresholdImageFilter
    threshold_filter = sitk.BinaryThresholdImageFilter()
        
    # Set the threshold values
    threshold_filter.SetLowerThreshold(min_intensity)
    max_value = sitk.GetArrayFromImage(image).max()
    threshold_filter.SetUpperThreshold(int(max_value))
        
    # Set the values for pixels inside and outside the threshold range
    threshold_filter.SetInsideValue(1)
    threshold_filter.SetOutsideValue(0)
        
    # Execute the filter to get the mask image
    return threshold_filter.Execute(image)


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

    ct1 = sitk.ReadImage(ct1_name)
    ct2 = sitk.ReadImage(ct2_name)
    mask1 = sitk.ReadImage(mask1_name)
    mask2 = sitk.ReadImage(mask2_name)

    elastixImageFilter = sitk.ElastixImageFilter()
    # parameterMap = sitk.GetDefaultParameterMap('translate')
    parameterMap = sitk.GetDefaultParameterMap('affine')
    parameterMap["DefaultPixelValue"] = ("-1024",)
    parameterMap["NumberOfSpatialSamples"] = ("1000",)
    parameterMap["MaximumNumberOfIterations"] = ("10000",)
    parameterMap["ErodeFixedMask"] = ("false",)
    parameterMap["MaximumNumberOfSamplingAttempts"] = ("8",)
    #parameterMap["NumberOfHistogramBins"] = ("256",)

    parameterMapVector = [
        #sitk.GetDefaultParameterMap("translation"),
        parameterMap
    ]
    elastixImageFilter.SetParameterMap(parameterMapVector)

    elastixImageFilter.SetFixedImage(ct1)
    # we can only use the mask if we dilate it as well, what about smoothing?
    bone_mask = create_mask_with_intensity_constraint2(ct1, 200)
    elastixImageFilter.SetMovingImage(ct2)
    elastixImageFilter.SetFixedMask(dilateIt(bone_mask, 7))
    elastixImageFilter.Execute()
    resultImage = elastixImageFilter.GetResultImage()
    transformParameterMap = elastixImageFilter.GetTransformParameterMap(0)
    sitk.WriteImage(resultImage, output + "/ct_moved_resampled.nii.gz")

    transformixImageFilter = sitk.TransformixImageFilter()
    transformParameterMap["ResampleInterpolator"] = ("FinalNearestNeighborInterpolator",)
    transformParameterMap["DefaultPixelValue"] = ("0",)
    transformixImageFilter.SetTransformParameterMap(transformParameterMap)

    transformixImageFilter.SetMovingImage(mask2)
    transformixImageFilter.Execute()
    sitk.WriteImage(transformixImageFilter.GetResultImage(), output + "/mask_moved_resampled.nii.gz")
    mask2_registered = transformixImageFilter.GetResultImage()

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

    sitk.WriteImage(split_mask, output + "/split_mask_moved_resampled.nii.gz")

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

    sitk.WriteImage(split_mask, output + "/split_mask_fixed.nii.gz")

    # save the volume info
    with open(output + "/nav_orbit_volumes.json", "w") as file:
        json.dump(volumes_per_region, file, indent=2)


if __name__ == "__main__":
    main()
