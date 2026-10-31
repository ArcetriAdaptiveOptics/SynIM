Overview
--------

SynIM (Synthetic Interaction Matrix) is a Python package designed for adaptive optics (AO) system analysis and design. It provides tools for computing:

- **Interaction matrices** between deformable mirrors (DMs) and wavefront sensors (WFS)
- **Projection matrices** for multi-conjugate adaptive optics (MCAO) tomography
- **Atmospheric covariance matrices** for optimal reconstructor design
- **MMSE reconstructors** with proper noise and turbulence covariance

Main Features
~~~~~~~~~~~~~

**GPU Acceleration**
   Built-in support for GPU acceleration via CuPy, with automatic fallback to CPU when GPU is unavailable or memory-limited.

**Multi-WFS Optimization**
   Intelligent batch processing that shares computations when multiple WFS view the same DM, providing 2-5× speedup.

**Flexible Configuration**
   Support for both YAML and PRO (IDL-style) configuration files with automatic parameter parsing and validation.

**SPECULA Integration**
   Compatible with the SPECULA adaptive optics simulation framework, using standard FITS file formats.

**Smart Caching**
   Automatic caching of computed matrices and intermediate results to avoid redundant computations.

**Workflow Optimization**
   Intelligent selection of computation workflows to minimize interpolation steps and maximize accuracy.

Architecture
------------

Core Components
~~~~~~~~~~~~~~~

SynIM is organized into several main modules:

**synim.py**
   Core functions for interaction matrix computation:
   
   - ``interaction_matrix()``: Main function, computed as a sparse operator applied to
     the DM cube (see :ref:`computation_workflows`)
   - ``interaction_matrices_multi_wfs()``: Optimized multi-WFS computation
   - ``apply_dm_transformations_combined()`` / ``apply_wfs_transformations_combined()``:
     explicit step-by-step computation, kept as reference and for debugging
   - ``compute_derivatives_with_extrapolation()``: Numerical derivatives with edge handling
   - ``compute_telsum_with_extrapolation()``: Telescoping sum estimation from phase differences per subaperture

   Slope extraction supports two methods controlled by ``slope_method`` argument of the main functions:

   - ``derivatives`` (default): numerical derivatives with extrapolation-aware edge handling
   - ``telsum`` (optional): average phase differences (telescoping sum) in each subaperture

**synpm.py**
   Core functions for projection matrix computation:
   
   - ``projection_matrix()``: Main function for DM-Layer projection
   - Similar workflow optimization as interaction matrices

**params_manager.py**
   High-level interface via the ``ParamsManager`` class:
   
   - Centralized parameter management
   - Batch computation of matrices
   - Automatic file naming and organization
   - MMSE reconstructor computation

**params_utils.py**
   Parameter handling utilities:
   
   - Configuration file parsing (YAML/PRO)
   - Parameter validation
   - Filename generation
   - Array transformations

**utils.py**
   General utility functions:

   - Array rebinning and masking
   - Geometric transformations
   - Zernike polynomials
   - FITS I/O helpers

**registration/**
   Unlike every module above, this does not compute a pixel/influence-
   function interaction or projection matrix. It is a purely geometric
   (2D affine) model of how a DM actuator grid lines up with a WFS
   sub-aperture grid in a system with several DMs and WFSs at different
   conjugation altitudes (MCAO/GLAO), following
   `Agapito, Plantet & Heritier (2024) <https://doi.org/10.48550/arXiv.2406.15336>`_:

   - ``geometry.py``: affine transforms (shift, rotation, magnification, anamorphosis)
   - ``model.py``: ``GuideStar``, ``DM``, ``WFS`` and ``System`` - the global mis-registration model
   - ``reconstruction.py``: global parameters from local measurements, by Gauss-Newton inversion
   - ``analysis.py``: sensitivity, degeneracy and noise-propagation analysis
   - ``viz.py``: altitude schematic, mis-registration table, top-down footprint views

   See :doc:`registration` for a full introduction and
   :doc:`tutorials/registration_morfeo_tutorial` for a worked example.

.. _computation_workflows:

Computation Workflows
~~~~~~~~~~~~~~~~~~~~~

**Interaction matrices**
   All the geometric transformations (DM shift, rotation and magnification from
   the guide star direction and height, WFS shift, rotation, magnification and
   anamorphosis) are combined into a single affine transformation, applied to
   the influence functions in one interpolation step. Conceptually, the slopes
   of one DM mode :math:`\varphi` are computed as follows:

   1. Combine all DM and WFS transformations into a single composite operation
      and apply it to the phase in one interpolation step
   2. Multiply by the transformed DM mask and repair the edges (the pixels just
      outside the valid region are extrapolated from the inside)
   3. Compute derivatives or G-tilts (telescoping sum) on the transformed grid
   4. Average over each subaperture and extract the slopes of the valid
      subapertures, in slope units

   The slope of a subaperture is the mean over its pupil pixels (derivatives)
   or over its pairs of adjacent pupil pixels (telescoping sum), so partially
   illuminated subapertures are averaged over their illuminated area. When the
   number of pixels is not a multiple of the number of subapertures (e.g. 480
   pixels and 68 subapertures), a pixel shared by two subapertures contributes
   to each with the fraction of its area inside it (``utils.rebin_matrix``):
   no interpolation of the phase or of the derivatives is needed.

   WFS transformations are never applied to derivative maps: the gradient of
   the transformed phase includes the Jacobian of the transformation (rotation
   and 45° anamorphosis mix the x and y components, magnification and 90°
   anamorphosis scale them), which a resampling of the derivative maps does not
   account for.

**The interaction matrix as a sparse operator**
   The masks do not depend on the DM mode. Once they are fixed, every step
   above is a *linear* function of the phase :math:`\varphi` (raveled to a
   vector of :math:`n_{pix}` pixels), so the whole chain is a product of
   sparse matrices:

   .. math::

      s = C \, P_{valid} \, M_{sa} \, G \, R \, D_{mask} \, T \, \varphi = A \, \varphi

   :math:`T`
      Bilinear interpolation of the combined DM + WFS transformation: each
      output pixel is a weighted sum of at most 4 input pixels. Built by
      ``utils.affine_interpolation_matrix`` from the matrix and offset of
      ``utils.rotshiftzoom_matrix`` (``_transformation_operator``).
   :math:`D_{mask}`
      Multiplication by the transformed DM mask (diagonal,
      ``_slope_operators``).
   :math:`R`
      Edge repair: the pixels inside the strict DM mask are kept, the edge
      pixels around it are extrapolated from up to 8 neighbours, all the
      others are 0 (``_repair_operator``). Derivatives only: the telescoping
      sum uses only pairs of valid pixels.
   :math:`G`
      x and y derivatives (central differences) or differences of adjacent
      pixels (telescoping sum), restricted to the valid pixels
      (``_slope_operators``).
   :math:`M_{sa}`
      Mean over each subaperture, with the pixel fractions of
      ``utils.rebin_matrix`` for non-integer ratios, and the subapertures
      outside the masks set to 0 (``_slope_operators``).
   :math:`P_{valid}`, :math:`C`
      Selection of the valid subapertures, x/y order of the convention and
      conversion to slope units (``_wfs_operator``).

   With the DM cube seen as a matrix :math:`\Phi` of size
   :math:`n_{pix} \times n_{modes}` (a reshape, not a copy), the interaction
   matrix is a single product:

   .. math::

      IM = A \, \Phi

   Why this is efficient:

   - :math:`A` does not depend on the number of modes: it is built once from
     the masks, and the memory used besides the input cube is that of
     :math:`A` and of the result. No transformed copy of the cube, or of its
     derivatives, is ever allocated.
   - :math:`A` is very sparse. The mean of the derivatives over a subaperture
     *telescopes*: along each row of a subaperture, the sum of the differences
     of adjacent pixels is the difference between the last and the first
     pixel, so the interior pixels cancel and only the pixels at the
     subaperture border keep a non-zero weight (the cancellation is partial
     only where the masks cut the subaperture). Each slope therefore depends
     on :math:`O(N)` pixels for :math:`N \times N` pixels per subaperture,
     not :math:`O(N^2)`. For 480×480 pixels and 40×40 subapertures,
     :math:`A` has about 200 000 non-zero elements (about 1.6 MB).
   - :math:`A` is built on CPU in float64 (``scipy.sparse``), including the
     transformed masks, and rounded once to the working precision, so CPU and
     GPU use the same operator. The product with the cube runs on the target
     device. On GPU it uses a small custom CuPy kernel (one thread per element
     of the result, summing each row in a fixed order) instead of cuSPARSE:
     the result is bit-wise reproducible, and the cube is read in its C order,
     without the Fortran-ordered copy that cuSPARSE would need.

   **Relation with the explicit computation.**
   ``apply_dm_transformations_combined`` (transformed phase and derivatives)
   and ``apply_wfs_transformations_combined`` (slopes from the derivatives)
   compute the same chain step by step on the transformed cube. They are kept
   permanently, as the independent reference for the sparse operator
   (``test/test_sparse_operator.py`` checks that the two agree for each
   geometric term, both slope methods and both conventions), and to inspect
   the transformed phase and derivatives when debugging.

   With several WFS (``interaction_matrices_multi_wfs``), the WFS with the same
   guide star, WFS transformations, number of subapertures and slope method
   share the transformation and slope operators, which are computed once for
   the group (see ``derivatives_info['groups']``). The operators of all the WFS
   are stacked and applied to the cube with a single product, so the cube is
   read once; each row is computed independently, so the results are identical
   to ``interaction_matrix``.

**Projection matrices**
   Projection matrices transform scalar fields (DM modes and the projection
   basis), so the DM and basis transformations can be applied separately
   (SEPARATED workflow) or combined in one interpolation step (COMBINED workflow,
   used when both the DM and the basis are transformed).

GPU Architecture
~~~~~~~~~~~~~~~~

GPU support is implemented through a flexible backend system:

**Automatic Backend Selection**
   .. code-block:: python
   
      import synim
      
      # CPU backend (numpy + scipy)
      synim.init(device_idx=-1, precision=1)
      
      # GPU backend (cupy + cupyx.scipy)
      synim.init(device_idx=0, precision=1)

   ``synim.init()`` must be called once, before importing the SynIM submodules
   (``synim.synim``, ``synim.synpm``, ``synim.utils``, ``synim.params_utils``,
   ``synim.params_manager``): they bind the array library and the data types
   when they are imported. Importing them before ``synim.init()`` raises a
   ``RuntimeError``, as does calling ``synim.init()`` again with a different
   configuration after they have been imported.

**Memory Management**
   - Automatic data transfer between CPU and GPU
   - Graceful fallback to CPU on memory errors
   - Efficient caching of GPU arrays

**Precision Control**
   - ``precision=0``: Double precision (float64)
   - ``precision=1``: Single precision (float32, ~2× faster on GPU)

Configuration Files
===================

SynIM uses **YAML** configuration files for defining AO system parameters, ensuring full compatibility with SPECULA. 

Legacy IDL-style parameter files (``.pro``) from `PASSATA <https://arxiv.org/abs/1607.07624>`_ are fully supported via an offline conversion tool. 

For details on how to format your YAML or how to convert your old ``.pro`` files, see the dedicated :ref:`Configuration <configuration>` page.

Loading Configurations
----------------------

Using ParamsManager
~~~~~~~~~~~~~~~~~~~

The recommended way to load configurations:

.. code-block:: python

   from synim.params_manager import ParamsManager
   
   # Load YAML file
   pm = ParamsManager('params_mcao.yml', verbose=True)
   
   # Access configuration
   print(f"Telescope diameter: {pm.params['telescope']['diameter']} m")
   print(f"Number of DMs: {len(pm.dm_list)}")
   print(f"Number of WFS: {len(pm.wfs_list)}")

Direct Parsing
~~~~~~~~~~~~~~

For custom workflows:

.. code-block:: python

   from synim.params_utils import parse_params_file
   
   # Parse any supported format
   params = parse_params_file('config.yml')

Validation
----------

SynIM automatically validates configurations:

.. code-block:: python

   from synim.params_utils import validate_opt_sources
   
   pm = ParamsManager('params.yml')
   
   # Validate optical sources for tomography
   validate_opt_sources(pm.params, verbose=True)
   # Output: Validating optical source configurations...
   #         ✓ All optical sources properly configured

Common validation checks:

- Required parameters present
- Valid ranges for physical quantities
- Consistent array dimensions
- Valid source types and positions
- DM-layer altitude compatibility

File Organization
-----------------

SynIM follows SPECULA directory structure for seamless integration. The ``ParamsManager`` automatically creates and manages these directories based on the ``root_dir`` parameter in your configuration:

.. code-block:: python

   # In ParamsManager initialization:
   self.im_dir = root_dir + '/synim/'       # Interaction matrices
   self.pm_dir = root_dir + '/synpm/'       # Projection matrices  
   self.rec_dir = root_dir + '/synrec/'     # Reconstructors
   self.cov_dir = root_dir + '/covariance/' # Covariance matrices

.. code-block:: text

   project/
   ├── config/
   │   ├── params_scao.yml
   │   ├── params_mcao.yml
   │   └── params_ltao.yml
   │
   └── calib/                   # Calibration data (root_dir)
       ├── synim/               # Interaction matrices (.fits)
       │   ├── intmat_wfs1_dm0.fits
       │   └── intmat_wfs2_dm0.fits
       │
       ├── synpm/               # Projection matrices (.fits)
       │   ├── projmat_dm0_layer0.fits
       │   └── projmat_dm1_layer1.fits
       │
       ├── synrec/              # Reconstructors (.fits)
       │   ├── rec_mmse.fits
       │   └── rec_lsq.fits
       │
       ├── covariance/          # Covariance matrices (.fits)
       │   ├── cov_atm.fits
       │   └── cov_noise.fits
       │
       ├── ifunc/               # Influence functions (SPECULA format)
       │   ├── dm0_ifunc.fits
       │   └── dm1_ifunc.fits
       ├── im/          # Interaction matrices (SPECULA format)
       └── rec/         # Reconstructors and Projection matrices (SPECULA format)


**Note:** When using SynIM with SPECULA, both tools can share the same ``root_dir`` and ``ifunc/`` directory. SPECULA uses additional directories (``im/``, ``rec/``) which can coexist alongside SynIM directories.

Filename Conventions
~~~~~~~~~~~~~~~~~~~~

SynIM automatically generates descriptive filenames based on component tags and parameters:

**Interaction Matrices:**
   ``intmat_{wfs_tag}_{dm_tag}.fits``
   
   Example: ``intmat_WFS_LGS1_DM0.fits``

**Projection Matrices:**
   ``projmat_{dm_tag}_layer{layer_idx}.fits``
   
   Example: ``projmat_DM0_layer0.fits``

**Reconstructors:**
   ``rec_{type}_{wfs_tags}_{dm_tags}.fits``
   
   Example: ``rec_mmse_WFS1-WFS2_DM0-DM1.fits``

**Covariance Matrices:**
   ``cov_atm_{layer_config}.fits`` or ``cov_noise_{wfs_config}.fits``

Custom Directories
~~~~~~~~~~~~~~~~~~

You can override the default directory structure:

.. code-block:: python

   pm = ParamsManager('params.yml')
   
   # Override specific directories
   pm.im_dir = '/custom/path/interaction_matrices/'
   pm.rec_dir = '/custom/path/reconstructors/'
   
   # Compute with custom paths
   pm.compute_all_interaction_matrices()

See Also
--------

- :doc:`installation` - Installation guide with GPU setup
- :doc:`api/index` - Complete API reference