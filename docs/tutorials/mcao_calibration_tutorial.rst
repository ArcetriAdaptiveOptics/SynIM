.. _morfeo_workflow:
.. _mcao_calibration_tutorial:

===========================================================================
Calibrating an MCAO System with SynIM: The MORFEO Methodology
===========================================================================

Introduction
============

Calibrating a Multi-Conjugate Adaptive Optics (MCAO) system for an Extremely Large Telescope (ELT) is one of the most challenging tasks in modern astronomy. MORFEO (Multi-conjugate Adaptive Optics Relay For ELT Observations) is a major scale-up from previous systems, with 12 wavefront sensors (WFSs), 3 deformable mirrors (DMs), and more than 6,000 controlled actuators.

This tutorial explains the calibration logic and the sequence of mathematical products required to close the loop with SynIM and SPECULA.

.. note::
   This is a **methodology tutorial**, not a usage guide for ``run_calib_workflow_morfeo.py``.
   That script can be used as a reference implementation of the same sequence, but the goal here is to let you implement your own calibration files/scripts for complex MCAO systems.
   
   The slope computation method used in the examples is the default 'derivatives', but you can also use 'gtilt' by changing the ``slope_method`` argument in the relevant function calls.
   'gtilt' provides the G-tilt (Gradient tilt) estimation, and 'derivatives' provides the Z-tilt (Zernike tilt, that is the global angular orientation of a wavefront) estimation.
   'gtilt' can provide more accurate slope estimation because the Shack-Hartmann spot displacements are not always well approximated by local derivatives, especially for high-order modes.

1. Theoretical Framework: The Three WFS Groups
==============================================

To deliver high and uniform wavefront correction across the science field of view, MORFEO uses a split-tomography approach, dividing sensing into three WFS groups:

* **Laser Guide Star (LGS) WFSs:** Six Shack-Hartmann sensors on sodium return, used as the primary high-order turbulence sensors.
* **Low-Order (LO) WFSs (NGS):** Three infrared sensors on Natural Guide Stars, used for global low-order terms (tip, tilt, and average focus control).
* **Reference (REF) WFSs (NGS):** Three visible/near-infrared truth sensors, run at low bandwidth to detrend pseudo-static aberrations and compensate LGS sensing artifacts (for example truncation and NCPA effects).

The LGS Filtering Challenge
---------------------------
On an ELT, the conversion from sodium-layer altitude variation to focus aberration is very large. Sodium dynamics can dominate the atmospheric focus signal.

Because of this, **MORFEO uses LGS measurements only for modes above tip-tilt-focus**. Tip, tilt, and focus are filtered out of the LGS path, while NGS-based channels provide robust low-order estimation.

2. Preparing the SPECULA Configuration
======================================

Before generating matrices with SynIM, your base SPECULA YAML (for example ``params_morfeo_simplified.yml``) must include MCAO-specific calibration sections.

Virtual DMs (Atmospheric Layers)
--------------------------------
Tomographic reconstruction requires estimating turbulence at multiple altitudes. In SPECULA, these are represented as virtual DMs (``layer1``, ``layer2``, ...). SynIM uses them to generate altitude-specific interaction matrices.

.. code-block:: yaml

    # Example from params_morfeo_simplified.yml
    layer1:
      class: 'DM'
      simul_params_ref: 'main'
      height: 0

    layer2:
      class: 'DM'
      simul_params_ref: 'main'
      height: 2000
    
    # ... additional layers up to the maximum sensing altitude (e.g., 22000 m)

Influence Functions and Modal Bases
-----------------------------------
Before interaction-matrix generation, ensure each physical DM and each virtual layer has:

* an influence-function object,
* an ``m2c_object``,
* a consistent ``nmodes`` value.

Without this step, matrix dimensions in the later calibration chain can become inconsistent.

Reconstruction Parameters
-------------------------
Define the parameters needed for minimum-mean-square-error (MMSE) reconstructors, including noise models and mode selection for each loop.

.. code-block:: yaml

    reconstruction:
      lgs:
        sigma2inNm2: 2e4
        noise_elong_model: false
        naThicknessInM: 10000.0
        tGparameter: 0.5
      ngs:
        sigma2inNm2: 1e4
      ref:
        sigma2inNm2: 1.11e5
      
      # Modes to use for NGS reconstruction (tip-tilt/focus split on DMs)
      ngs_n_modes_dm:  [2, 0, 3]
      
      # Modes to use for REF reconstruction (truth sensing)
      ref_n_modes_dm:  [50, 45, 48]

Projection Parameters
---------------------
The projection operator finds DM commands that maximize correction over the target field. Define regularization and optimization directions (with weights).

.. code-block:: yaml

    projection:
      reg_factor: 1e-4
      opt_sources:
        polar_coordinates: 
          - [0.0, 0.0]
          - [30.0, 0.0]
          # ... further coordinates for inner/outer rings
        weights: 
          - 0.5
          - 1.0
          # ... corresponding weights to prioritize specific regions

3. The SynIM Calibration Sequence
=================================

With YAML prepared, the SynIM calibration follows a fixed chain of products.
Each step is typically implemented as a short standalone script that you run once
before starting SPECULA simulations.

Step 0: Influence Functions and Modal Bases
-------------------------------------------
Compute/assign influence functions and modal bases for:

* physical DMs,
* virtual atmospheric layers,
* (optionally) inverse modal basis used by downstream analysis/projection blocks.

This is a prerequisite for stable and dimensionally consistent IM/reconstructor generation.

Create a script ``step0_influence_functions.py``:

.. code-block:: python

   import specula
   specula.init(-1)  # -1 for CPU

   from specula.lib.make_mask import make_mask
   from specula.lib.compute_zonal_ifunc import compute_zonal_ifunc
   from specula.lib.modal_base_generator import make_modal_base_from_ifs_fft
   from specula.data_objects.ifunc import IFunc
   from specula.data_objects.ifunc_inv import IFuncInv
   from specula.data_objects.m2c import M2C
   from specula.data_objects.pupilstop import Pupilstop
   from specula.data_objects.simul_params import SimulParams
   from specula import np, xp, to_xp, cpuArray

   CALIB_DIR = '/raid1/guido/PASSATA/MORFEOtest/'   # must match main.root_dir in your YAML
   pupil_pixels = 160
   telescope_diameter_m = 38.5        # ELT diameter [m]
   pixel_size = telescope_diameter_m / pupil_pixels
   obsratio = 0.283
   n_petals = 6
   fov_arcsec = 160
   overwrite = True

   # DMs then virtual layers: (altitude [m], n_actuators)
   altitudes_m = [600, 7500, 17500,   0, 2000, 4500, 11000, 15000, 18000, 22000]
   nacts        = [ 41,   37,    37,  41,   41,   41,    37,    37,    37,    37]

   ARCSEC2RAD = np.pi / (180 * 3600)
   meta_size_m  = fov_arcsec * ARCSEC2RAD * np.array(altitudes_m) + telescope_diameter_m
   meta_size_px = (np.ceil(meta_size_m / pixel_size / 2) * 2).astype(int)

   # Pupil mask with petals (used as ground-layer mask)
   pupil_mask = make_mask(pupil_pixels, obsratio=obsratio, spider=True,
                          spider_width=1, n_petals=n_petals, xp=np)
   simul_params = SimulParams(pixel_pupil=pupil_pixels, pixel_pitch=pixel_size)
   Pupilstop(simul_params, input_mask=pupil_mask).save(
       f'{CALIB_DIR}/pupilstop/mask_{pupil_pixels}px_{n_petals}petals.fits')

   for alt, n_act, meta_px in zip(altitudes_m, nacts, meta_size_px):
       base_name = f"base_{meta_px}px_{n_act}acts"
       mask = None
       if alt == 0:
           base_name += f"_{n_petals}petals"
           mask = to_xp(xp, cpuArray(pupil_mask))

       influence_functions, meta_pupil_mask, _, _ = compute_zonal_ifunc(
           int(meta_px), n_act, circ_geom=True, angle_offset=0,
           mask=mask, xp=xp, dtype=xp.float32)

       kl_basis, m2c, _ = make_modal_base_from_ifs_fft(
           meta_pupil_mask, telescope_diameter_m, influence_functions,
           r0=0.15, L0=25.0, zern_modes=5, oversampling=2,
           if_max_condition_number=1e4, xp=xp, dtype=xp.float32, verbose=True)

       n_modes = kl_basis.shape[0]
       IFunc(ifunc=influence_functions, mask=meta_pupil_mask).save(
           f"{CALIB_DIR}/ifunc/{base_name}.fits", overwrite=overwrite)
       M2C(m2c=m2c).save(
           f"{CALIB_DIR}/m2c/{base_name}_{n_modes}modes.fits", overwrite=overwrite)

       if alt == 0:
           # Inverse modal basis: ground layer only — used by projection/analysis blocks
           kl_basis_inv = np.linalg.pinv(cpuArray(kl_basis))
           IFuncInv(ifunc_inv=kl_basis_inv, mask=pupil_mask).save(
               f"{CALIB_DIR}/ifunc/{base_name}_{n_modes}modes_inv.fits",
               overwrite=overwrite)

.. code-block:: bash

   python step0_influence_functions.py

After this step, update the ``ifunc_object``, ``m2c_object``, and ``nmodes`` fields for
every DM and layer entry in your YAML, and set ``pupilstop.tag`` and
``projection.ifunc_inverse_tag`` accordingly.

Step A: Subaperture Calibration
-------------------------------
Calibrate subaperture geometry for all WFS families (LGS, NGS, REF). This defines
valid pixels and thresholds for slope extraction, and produces ``subapdata`` objects
to assign to slope processors.

SPECULA handles this step. Create a dedicated calibration YAML
(e.g. ``calib_morfeo_simplified_subaps.yml``) and run:

.. code-block:: bash

   specula params_morfeo_simplified.yml calib_morfeo_simplified_subaps.yml

After this step, update ``subapdata_object`` in every ``slopec_*`` section of your
YAML with the generated tags before proceeding.

Step B: Interaction Matrices (IMs)
----------------------------------
Compute separate interaction matrices for each control branch:

* **LGS IMs:** over all virtual layers (tomographic sensing volume).
* **NGS (LO) IMs:** directly on physical DMs (low-order control path).
* **REF IMs:** typically on Layer 0 for low-order truth/focus handling.

Create a script ``step2_compute_ims.py``:

.. code-block:: python

   import specula
   specula.init(device_idx=0, precision=1)
   import synim
   synim.init(device_idx=0, precision=1)
   from synim.params_manager import ParamsManager

   yaml_file = 'params_morfeo_simplified.yml'
   root_dir  = '/raid1/guido/PASSATA/MORFEOtest/'
   slope_method = 'derivatives'  # default method for slope extraction (the other option is 'gtilt')

   params_mgr = ParamsManager(yaml_file, root_dir=root_dir, verbose=True)

   # LGS IMs — one per WFS × layer combination
   lgs_ims = params_mgr.compute_interaction_matrices(
       wfs_type='lgs',
       output_im_dir=root_dir + 'synim/',
       output_rec_dir=root_dir + 'synrec/',
       slope_method=slope_method,
       overwrite=False, verbose=True, display=False)
   print(f"LGS: {len(lgs_ims)} IMs")

   # NGS (LO) IMs — one per WFS × DM
   ngs_ims = params_mgr.compute_interaction_matrices(
       wfs_type='ngs',
       output_im_dir=root_dir + 'synim/',
       output_rec_dir=root_dir + 'synrec/',
       slope_method=slope_method,
       overwrite=False, verbose=True, display=False)
   print(f"NGS: {len(ngs_ims)} IMs")

   # REF IMs — Layer 0 only (focus / truth sensing)
   ref_ims = params_mgr.compute_interaction_matrices(
       wfs_type='ref',
       output_im_dir=root_dir + 'synim/',
       output_rec_dir=root_dir + 'synrec/',
       slope_method=slope_method,
       overwrite=False, verbose=True, display=False)
   print(f"REF: {len(ref_ims)} IMs")

.. code-block:: bash

   python step2_compute_ims.py

Step C: Filter Matrices for LGS
--------------------------------
For ELT-scale systems, LGS measurements must exclude tip-tilt-focus. Build filter
matrices from Layer 0 LGS IMs and apply them to LGS slope streams.

On smaller telescopes, filtering may be limited to tip-tilt only.

.. note::
   **RTC vs. Simulation Implementation**

   In a physical real-time computer (RTC), this filter is often embedded in the final
   command matrix. In SPECULA simulation, it is typically applied at slope level via
   ``filtmat_data`` in LGS slope objects.

MORFEO has six LGS WFSs at three unique rotation angles (6.2°, −6.2°, 14.2°), so
only three filter matrices are needed — one per unique configuration. Check the
``synim/`` directory for the exact Layer-0 IM filenames produced by Step B.

Create a script ``step3_generate_filtmat.py``:

.. code-block:: python

   from synim.params_utils import generate_filter_matrix_from_intmat_file

   root_dir  = '/raid1/guido/PASSATA/MORFEOtest/'
   synim_dir = root_dir + 'synim/'
   out_dir   = root_dir + 'data/'

   n_modes          = 1000  # modes used for the filter projection
   n_modes_filtered = 3     # tip, tilt, focus

   # Map each unique rotation to its Layer-0 IM and desired output name.
   # Replace the IM filenames with the exact names produced by Step B.
   configs = [
       ('IM_syn_..._rot6.2_layH0.0_....fits',  'filtmat_rot6.2_mn1000_3.fits'),
       ('IM_syn_..._rot-6.2_layH0.0_....fits', 'filtmat_rot-6.2_mn1000_3.fits'),
       ('IM_syn_..._rot14.2_layH0.0_....fits', 'filtmat_rot14.2_mn1000_3.fits'),
   ]

   for im_name, filtmat_name in configs:
       generate_filter_matrix_from_intmat_file(
           intmat_filename=synim_dir + im_name,
           n_modes=n_modes,
           n_modes_filtered=n_modes_filtered,
           output_filename=out_dir + filtmat_name,
           smooth_tt=True,
           overwrite=True,
           verbose=True)

.. code-block:: bash

   python step3_generate_filtmat.py

After this step, update ``filtmat_data`` in every ``slopec_lgs*`` section of your YAML
with the matching filtmat filename before computing reconstructors.

Step D: Tomographic Reconstructors
------------------------------------
Using IMs and ``reconstruction`` parameters, compute reconstructors:

* **LGS Reconstructor:** filtered LGS slopes → virtual layers.
* **NGS Reconstructor:** LO slopes → physical DMs.
* **Focus Reconstructor:** built from REF Layer-0 IMs, stacking all REF channels and
  applying pseudoinverse.

Create a script ``step4_compute_reconstructors.py``:

.. code-block:: python

   import numpy as np
   import specula
   specula.init(device_idx=0, precision=1)
   import synim
   synim.init(device_idx=0, precision=1)
   from synim.params_manager import ParamsManager
   from specula.data_objects.intmat import Intmat
   from specula.data_objects.recmat import Recmat

   yaml_file = 'params_morfeo_simplified.yml'
   root_dir  = '/raid1/guido/PASSATA/MORFEOtest/'
   rec_dir   = root_dir + 'synrec/'
   synim_dir = root_dir + 'synim/'
   r0, L0    = 0.15, 25.0
   slope_method = 'derivatives'  # default method for slope extraction (the other option is 'gtilt')

   params_mgr = ParamsManager(yaml_file, root_dir=root_dir, verbose=True)

   # --- LGS tomographic reconstructor (filtered slopes → virtual layers) ---
   result_lgs = params_mgr.compute_tomographic_reconstructor(
       r0=r0, L0=L0,
       wfs_type='lgs', component_type='layer',
       noise_variance=None,   # computed automatically from detector parameters
       slope_method=slope_method,  # default method for slope extraction
       output_dir=rec_dir,
       save=True, verbose=True)
   print(f"LGS reconstructor: {result_lgs['reconstructor'].shape}")

   # Save the assembled (tip-tilt-focus filtered) LGS IM for POLC reference
   params_mgr.save_assembled_interaction_matrix(
       wfs_type='lgs', component_type='dm',
       output_dir=synim_dir, overwrite=True,
       apply_filter=True, slope_method=slope_method,
       verbose=True)

   # --- NGS reconstructor (LO slopes → DMs) ---
   result_ngs = params_mgr.compute_tomographic_reconstructor(
       r0=r0, L0=L0,
       wfs_type='ngs', component_type='dm',
       noise_variance=None,
       output_dir=rec_dir,
       slope_method=slope_method,
       save=True, verbose=True)
   print(f"NGS reconstructor: {result_ngs['reconstructor'].shape}")

   # --- Focus reconstructor (REF slopes → focus modes, 3 WFSs stacked) ---
   # Replace the filename with the exact REF Layer-0 IM produced by Step B.
   ref_im_path   = synim_dir + 'IM_syn_..._ref1_layH0.0_....fits'
   n_focus_modes = 30
   n_ref_wfs     = 3

   intmat_obj = Intmat.restore(ref_im_path)
   stacked    = np.vstack([intmat_obj.intmat] * n_ref_wfs)
   recmat_obj = Intmat(intmat=stacked).generate_rec(
       nmodes=n_focus_modes, cut_modes=0, interactive=False)
   Recmat(recmat=recmat_obj.recmat,
          norm_factor=recmat_obj.norm_factor).save(
       rec_dir + f'focus_recmat_nmodes{n_focus_modes}_nwfs{n_ref_wfs}',
       overwrite=True)
   print(f"Focus recmat: {recmat_obj.recmat.shape}")

.. code-block:: bash

   python step4_compute_reconstructors.py

After this step, update ``recmat_object`` and ``intmat_object`` in your YAML
(``tomo_polc_lgs``, ``tomo_ngs``, ``rec_focus``) with the filenames just produced.

Step E: Projection Matrices
----------------------------
The final step maps layered turbulence estimates back to physical DMs. Compute
per DM-layer projection terms, then combine them into a single tomographic
projection matrix optimized on the field and weights defined in
``projection.opt_sources``.

Create a script ``step5_compute_projection.py``:

.. code-block:: python

   import specula
   specula.init(device_idx=0, precision=1)
   import synim
   synim.init(device_idx=0, precision=1)
   from synim.params_manager import ParamsManager

   yaml_file = 'params_morfeo_simplified.yml'
   root_dir  = '/raid1/guido/PASSATA/MORFEOtest/'
   pm_dir    = root_dir + 'synpm/'

   params_mgr = ParamsManager(yaml_file, root_dir=root_dir, verbose=True)

   # Individual per DM-layer projection matrices
   pm_paths = params_mgr.compute_projection_matrices(
       output_dir=pm_dir, overwrite=False)
   print(f"Computed {len(pm_paths)} projection matrices")

   # Combine into optimised tomographic projection matrix
   p_opt, _, _, info = params_mgr.compute_tomographic_projection_matrix(
       output_dir=pm_dir, save=True, verbose=True)

   print(f"Tomographic projection matrix: {p_opt.shape}")
   print(f"  Optical sources: {info['n_opt_sources']}")
   print(f"  Regularization:  {info['reg_factor']}")

.. code-block:: bash

   python step5_compute_projection.py

After this step, update ``projmat_object`` in ``tomo_polc_lgs`` with the
filename just produced.

Summary of Outputs
==================
A successful SynIM MCAO calibration provides:

1. Valid subaperture maps for all WFS families.
2. Interaction matrices for LGS tomography, NGS low-order path, and REF truth/focus path.
3. LGS filter matrices (tip-tilt-focus removal for ELT use case).
4. Reconstructors for LGS, NGS, and focus branches.
5. A tomographic projection matrix from virtual layers to physical DMs.

These products are the minimum set required to build and close a complete MCAO control chain in SPECULA.
