"""
Frozen reference implementation of ``synim.utils.rotshiftzoom_array``.

This is a verbatim copy of ``rotshiftzoom_array`` as of commit 12ea3de
(the implementation based on ``affine_transform``), kept here so that the
equivalence tests keep comparing against the original behaviour even if the
library code changes in the future.

DO NOT EDIT the function body. The only differences from the original are
the function name and that ``xp`` and ``affine_transform`` are read from
``synim`` at call time, so the reference always runs on the active backend.
"""
import synim


def rotshiftzoom_array_reference(input_array, dm_translation=(0.0, 0.0),
                       dm_rotation=0.0, dm_magnification=(1.0, 1.0),
                       wfs_translation=(0.0, 0.0), wfs_rotation=0.0,
                       wfs_magnification=(1.0, 1.0),
                       wfs_anamorphosis_45=1.0,
                       output_size=None):
    """
    This function applies magnification, rotation, shift and resize of a
    2D or 3D numpy/cupy array using affine transformation.
    Rotation is applied in the same direction as the first function.

    Parameters:
    - input_array: numpy/cupy array, input data to be transformed
    - dm_translation: tuple, translation for DM (x, y)
    - dm_rotation: float, rotation angle for DM in degrees
    - dm_magnification: tuple, magnification factors for DM (x, y)
    - wfs_translation: tuple, translation for WFS (x, y)
    - wfs_rotation: float, rotation angle for WFS in degrees
    - wfs_magnification: tuple, magnification factors for WFS (x, y)
    - wfs_anamorphosis_45 (float): Diagonal anamorphosis factor.
      Values > 1 stretch along +45° diagonal and compress along -45° diagonal.
      Implemented as a shear transformation.
    - output_size: tuple, desired output size (height, width)

    Returns:
    - output: numpy/cupy array, transformed data
    """
    xp = synim.xp
    affine_transform = synim.affine_transform


    # Parameter handling: conversion of single values to tuples
    try:
        if not hasattr(dm_translation, '__len__') or len(dm_translation) != 2:
            dm_translation = (float(dm_translation), float(dm_translation))
    except (TypeError, ValueError):
        dm_translation = (0.0, 0.0)

    try:
        if not hasattr(wfs_translation, '__len__') or len(wfs_translation) != 2:
            wfs_translation = (float(wfs_translation), float(wfs_translation))
    except (TypeError, ValueError):
        wfs_translation = (0.0, 0.0)

    try:
        if not hasattr(dm_magnification, '__len__'):
            # If it is a single value, we create a tuple with two identical elements
            dm_magnification = (float(dm_magnification), float(dm_magnification))
        elif len(dm_magnification) != 2:
            # if it is a sequence but not of length 2
            dm_magnification = (float(dm_magnification[0]), float(dm_magnification[0]))
    except (TypeError, ValueError):
        dm_magnification = (1.0, 1.0)

    try:
        if not hasattr(wfs_magnification, '__len__'):
            # If it is a single value, we create a tuple with two identical elements
            wfs_magnification = (float(wfs_magnification), float(wfs_magnification))
        elif len(wfs_magnification) != 2:
            # if it is a sequence but not of length 2
            wfs_magnification = (float(wfs_magnification[0]), float(wfs_magnification[0]))
    except (TypeError, ValueError):
        wfs_magnification = (1.0, 1.0)

    if xp.isnan(input_array).any():
        input_array = xp.nan_to_num(input_array, copy=True, nan=0.0, posinf=None, neginf=None)

    # Check if array is 2D or 3D
    is_3d = len(input_array.shape) == 3

    # resize
    if output_size is None:
        output_size = input_array.shape[:2]  # Only take the first two dimensions

    # Center of the input array
    center = xp.array(input_array.shape[:2]) / 2.0
    # Convert rotations to radians
    # Note: Inverting the sign of rotation to match the first function's direction
    dm_rot_rad = xp.deg2rad(-dm_rotation)  # Negative sign to reverse direction
    wfs_rot_rad = xp.deg2rad(-wfs_rotation)  # Negative sign to reverse direction
    # Initialize the output array
    if is_3d:
        output = xp.zeros((output_size[0], output_size[1], input_array.shape[2]),
                          dtype=input_array.dtype)
    else:
        output = xp.zeros(output_size, dtype=input_array.dtype)

    # Create the transformation matrices
    # For DM transformation
    dm_scale_matrix = xp.array(
        [[1.0/dm_magnification[0], 0], [0, 1.0/dm_magnification[1]]]
    )
    dm_rot_matrix = xp.array(
        [[xp.cos(dm_rot_rad), -xp.sin(dm_rot_rad)], [xp.sin(dm_rot_rad), xp.cos(dm_rot_rad)]]
    )
    dm_matrix = xp.dot(dm_rot_matrix, dm_scale_matrix)

    # For WFS transformation
    wfs_scale_matrix = xp.array(
        [[1.0/wfs_magnification[0], 0], [0, 1.0/wfs_magnification[1]]]
    )
    wfs_rot_matrix = xp.array(
        [[xp.cos(wfs_rot_rad), -xp.sin(wfs_rot_rad)],
         [xp.sin(wfs_rot_rad), xp.cos(wfs_rot_rad)]]
    )

    # OPTICAL ORDER: The Shack-Hartmann sensor rotation is the FINAL step 
    # of the physical optical path. Therefore, in Inverse Mapping, we MUST un-rotate FIRST.
    # xp.dot(A, B) applies B first. So wfs_rot_matrix MUST be on the right!
    if wfs_anamorphosis_45 != 1.0:
        k = wfs_anamorphosis_45
        anam_45_matrix = xp.array([
            [(1.0 + k) / 2.0, (1.0 - k) / 2.0],
            [(1.0 - k) / 2.0, (1.0 + k) / 2.0]
        ])
        wfs_matrix = xp.dot(anam_45_matrix, xp.dot(wfs_scale_matrix, wfs_rot_matrix))
    else:
        wfs_matrix = xp.dot(wfs_scale_matrix, wfs_rot_matrix)

    # Combine transformations
    # INVERSE MAPPING ALGEBRA: To match a 2-step process (DM then WFS),
    # the inverse matrix must be M_DM_inv * M_WFS_inv.
    combined_matrix = xp.dot(dm_matrix, wfs_matrix)

    # For 3D arrays, extend the transformation matrix to 3x3
    if is_3d:
        combined_matrix_3d = xp.eye(3)
        combined_matrix_3d[:2, :2] = combined_matrix
        combined_matrix = combined_matrix_3d

    # Calculate offset
    output_center = xp.array(output_size) / 2.0

    # The mathematical proof of the 2-step equivalence requires:
    # 1. WFS translation is multiplied by the FULL combined matrix.
    # 2. DM translation is multiplied ONLY by the DM rotation matrix.
    rotated_dm_translation = xp.dot(dm_rot_matrix, xp.array(dm_translation)[:2])
    spatial_combined = combined_matrix[:2, :2] if is_3d else combined_matrix
    scaled_wfs_translation = xp.dot(spatial_combined, xp.array(wfs_translation)[:2])

    if is_3d:
        offset_2d = center[:2] - xp.dot(spatial_combined, output_center) \
            - rotated_dm_translation - scaled_wfs_translation
        offset = xp.zeros(3, dtype=offset_2d.dtype)
        offset[:2] = offset_2d
    else:
        offset = center - xp.dot(combined_matrix, output_center) \
            - rotated_dm_translation - scaled_wfs_translation

    # Apply transformation (scipy requires numpy)
    output = affine_transform(
        input_array,
        combined_matrix,
        offset=offset,
        output_shape=output_size if not is_3d else output_size + (input_array.shape[2],),
        order=1
    )

    return output
