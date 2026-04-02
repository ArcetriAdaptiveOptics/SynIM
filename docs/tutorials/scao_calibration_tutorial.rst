.. _scao_calibration_tutorial:

=====================================
SCAO Calibration Tutorial
=====================================

This tutorial demonstrates how to calibrate an interaction matrix for a **Single Conjugate Adaptive Optics (SCAO)** system using SynIM, integrated with a SPECULA simulation.

.. note::
   This tutorial focuses on **interaction matrix computation**. It assumes you have already configured **Modal basis or influence functions** and ``ShSubapCalibrator`` (see `SPECULA calibration documentation <https://specula.readthedocs.io/>`_)

Prerequisites
=============

Before starting this tutorial, you should have:

1. **SPECULA simulation configuration**: A working YAML file defining your SCAO system (telescope, atmosphere, WFS, DM, etc.)
2. **Python environment** with SPECULA and SynIM installed
3. **Basic knowledge** of adaptive optics concepts

System Overview
===============

We'll use a simple SCAO configuration with:

- **Telescope**: 1m diameter (64 pixels @ 0.015625 m/pixel)
- **Deformable Mirror**: Ground-conjugated (height = 0 m), 40 Zernike modes
- **Wavefront Sensor**: 8×8 Shack-Hartmann, 600 nm wavelength
- **Guide Star**: On-axis Natural Guide Star (NGS)

YAML Configuration Files
=========================

The SPECULA SCAO calibration workflow uses three YAML files:

1. **Base simulation** (``params_scao.yml``): Main system configuration
2. **Subaperture calibration** (``params_scao_subap.yml``): WFS geometry calibration
3. **Interaction matrix calibration** (``params_scao_rec.yml``): IM computation settings

Here we will focus on SynIM which requires only the first parameter file for interaction matrix generation.
A prerequisite is the subaperture geometry calibration (Step 1) which produces a data object that is referenced in the main YAML file.
For a description of the subaperture calibration, see the `SPECULA SCAO basic tutorial <https://specula.readthedocs.io/en/latest/tutorials/scao_basic_tutorial.html/>`_.

Base Simulation Configuration
------------------------------

The main configuration file defines the complete SCAO system. Below is an excerpt showing the key components for interaction matrix calibration:

.. code-block:: yaml

   # params_scao.yml (excerpt)
   main:
     class:             'SimulParams'
     root_dir:          './calib/'
     pixel_pupil:       64
     pixel_pitch:       0.015625  # meters
     total_time:        0.010
     time_step:         0.001
   
   on_axis_source:
     class:             'Source'
     polar_coordinates:  [0.0, 0.0]  # On-axis
     magnitude:         5
     wavelengthInNm:    600
   
   pupilstop:
     class: 'Pupilstop'
     simul_params_ref: 'main'
   
   sh:
     class:             'SH'
     subap_on_diameter: 8
     subap_wanted_fov:  4.0  # arcsec
     sensor_pxscale:    0.5  # arcsec/pix
     subap_npx:         8
     wavelengthInNm:    600
     inputs:
       in_ef: 'prop.out_on_axis_source_ef'
   
   detector:
     class:             'CCD'
     simul_params_ref:  'main'
     size:              [64, 64]
     dt:                0.001
     bandw:             300
     quantum_eff:       0.3
     photon_noise:      True
     readout_noise:     True
     readout_level:     1.0
     inputs:
       in_i: 'sh.out_i'
   
   slopec:
     class:             'ShSlopec'
     subapdata_object:  'scao_subaps_n8_th0.5'  # From Step 1
     weightedPixRad:    4.0
     inputs:
       in_pixels: 'detector.out_pixels'
   
   dm:
     class:             'DM'
     simul_params_ref:  'main'
     type_str:          'zernike'
     nmodes:            40
     obsratio:          0.1
     height:            0  # Ground conjugated
     inputs:
       in_command: 'control.out_comm'

.. note::
   This is a simplified configuration showing only the components needed for interaction matrix computation. A complete SCAO simulation would also include:
   
   - ``atmo``: Atmospheric turbulence layers
   - ``prop``: Wavefront propagation
   - ``rec``: Mode reconstructor
   - ``control``: AO loop controller (e.g., integrator)
   - ``psf``: PSF computation for performance metrics
   
   See the `full example file <https://github.com/INAF-OAA/SynIM/blob/main/test/params_scao_sh_test.yml>`_ for a complete configuration.

Workflow Steps
==============

Step 1: Generate Subaperture Geometry
--------------------------------------

First, calibrate the WFS subaperture geometry using SPECULA.
For a description of this step, see the `SPECULA SCAO basic tutorial <https://specula.readthedocs.io/en/latest/tutorials/scao_basic_tutorial.html/>`_.


Step 2: Compute Interaction Matrix with SynIM
----------------------------------------------

Now compute the same IM using SynIM's synthetic approach:

.. code-block:: python

   from synim.params_manager import ParamsManager

   base_yml = 'params_scao.yml'
   
   # Initialize ParamsManager
   params_mgr = ParamsManager(
       base_yml,
       root_dir=calib_dir,
       verbose=True
   )
   
   # Compute interaction matrix
   print("Computing IM with SynIM...")
   synim_im = params_mgr.compute_interaction_matrix(
       wfs_type='ngs',      # WFS type: 'ngs', 'lgs', or 'ref'
       wfs_index=None,      # Auto-detect first NGS WFS
       dm_index=1,          # DM index from YAML
       verbose=True,
       display=False        # Set True to visualize
   )
   
   print(f"Interaction matrix shape: {synim_im.shape}")
   print(f"  Slopes: {synim_im.shape[0]}")
   print(f"  Modes:  {synim_im.shape[1]}")

**Output:**

.. code-block:: text

   Interaction matrix shape: (100, 40)
     Slopes: 100  # 2 × (number of valid subapertures)
     Modes:  40   # Number of DM modes

.. note::
   SynIM automatically selects the optimal computation workflow (SEPARATED or COMBINED) based on your system geometry. For on-axis SCAO with no WFS transformations, the SEPARATED workflow is typically used.
   
   For details on workflow selection logic and when each is optimal, see :ref:`Computation Workflows <computation_workflows>` in the General Documentation.


Step 3: Generate Reconstruction Matrix
---------------------------------------

Create the reconstructor from the interaction matrix:

.. code-block:: python

   from specula.data_objects.intmat import Intmat
   from specula.data_objects.recmat import Recmat
   
   # Create Intmat object
   intmat_obj = Intmat(
       synim_im,
       pupdata_tag='scao_synim_n8_th0.5',
       norm_factor=1.0
   )
   
   # Save interaction matrix
   im_output = os.path.join(calib_dir, 'im', 'scao_im_synim.fits')
   intmat_obj.save(im_output, overwrite=True)
   print(f"Saved IM: {im_output}")
   
   # Generate reconstruction matrix (pseudoinverse)
   recmat_obj = intmat_obj.generate_rec(
       nmodes=40,          # Number of modes to reconstruct
       cut_modes=0,        # Modes to exclude (e.g., piston)
       interactive=False   # Set True for interactive mode selection
   )
   
   # Save reconstruction matrix
   rec_output = os.path.join(calib_dir, 'rec', 'scao_rec_synim.fits')
   recmat_obj.save(rec_output, overwrite=True)
   print(f"Saved reconstructor: {rec_output}")
   print(f"Reconstructor shape: {recmat_obj.recmat.shape}")

**Use in SPECULA Simulation**

Update your YAML to reference the computed reconstructor:

.. code-block:: yaml

   rec:
     class: 'Modalrec'
     recmat_object: 'scao_rec_synim'  # Without .fits extension
     inputs:
       in_slopes: 'slopec.out_slopes'
     outputs: ['out_modes', 'out_pseudo_ol_modes']


Advanced Topics
===============

WFS Transformations
-------------------

For systems with WFS rotation or pupil shifts:

.. code-block:: yaml

   sh:
     class: 'SH'
     # ... other parameters ...
     rotAnglePhInDeg: 15.0      # WFS rotation [degrees]
     xShiftPhInPixel: 5.0       # X shift [pixels]
     yShiftPhInPixel: 2.0       # Y shift [pixels]

SynIM automatically handles these transformations during IM computation.

**Example**: Test with rotation

.. code-block:: python

   # Use rotated configuration
   base_yml_rot = 'params_scao_rot.yml'
   
   params_mgr_rot = ParamsManager(
       base_yml_rot,
       root_dir=calib_dir,
       verbose=True
   )
   
   synim_im_rot = params_mgr_rot.compute_interaction_matrix(
       wfs_type='ngs',
       dm_index=1
   )

GPU Acceleration
----------------

For large systems, enable GPU:

.. code-block:: python

   import synim
   synim.init(device_idx=0, precision=1)  # GPU 0, single precision
   
   # ParamsManager will automatically use GPU
   params_mgr = ParamsManager(base_yml, root_dir=calib_dir)


Troubleshooting
===============

Debug Mode
----------

Enable detailed logging:

.. code-block:: python

   params_mgr = ParamsManager(
       base_yml,
       root_dir=calib_dir,
       verbose=True  # Detailed progress
   )
   
   im = params_mgr.compute_interaction_matrix(
       wfs_type='ngs',
       dm_index=1,
       verbose=True,
       display=True  # Show intermediate plots
   )

Summary
=======

This tutorial covered:

1. **YAML configuration** for SCAO system
2. **Subaperture calibration** with SPECULA
3. **Interaction matrix computation** with SynIM
4. **Reconstructor generation** for closed-loop control
5. **Advanced topics**: transformations, GPU acceleration

For more complex systems (MCAO, LTAO), see:

- :ref:`MORFEO workflow example <morfeo_workflow>` for multi-WFS tomographic calibration
- :doc:`../api/params_manager` for API details
- `SPECULA Documentation <https://specula.readthedocs.io/>`_ for calibration prerequisites
