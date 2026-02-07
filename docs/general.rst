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
   
   - ``interaction_matrix()``: Main function with intelligent workflow selection
   - ``interaction_matrices_multi_wfs()``: Optimized multi-WFS computation
   - ``compute_derivatives_with_extrapolation()``: Numerical derivatives with edge handling

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

.. _computation_workflows:

Computation Workflows
~~~~~~~~~~~~~~~~~~~~~

SynIM implements two main workflows for interaction matrix computation, automatically selected based on system geometry to optimize accuracy and performance.

**SEPARATED Workflow**
   Used when transformations exist only in DM or WFS (not both):
   
   **Process:**
   
   1. Apply DM transformations to influence functions (rotation, magnification, source altitude projection)
   2. Compute numerical derivatives with edge extrapolation
   3. Apply WFS transformations to derivatives (rotation, translation, magnification)
   4. Bin derivatives to subaperture resolution
   5. Extract slopes for valid subapertures
   
   **Advantages:**
   
   - Single interpolation step per transformation type
   - Maximum accuracy (no cumulative interpolation errors)
   - Can reuse derivatives for multiple WFS with same geometry

**COMBINED Workflow**
   Used when both DM and WFS have transformations:
   
   **Process:**
   
   1. Combine all DM and WFS transformations into single composite operation
   2. Apply combined transformation to influence functions in one interpolation step
   3. Compute derivatives on final transformed grid
   4. Bin and extract slopes
   
   **Advantages:**
   
   - Avoids double interpolation artifacts that occur when transformations are coupled
   - Single interpolation step when both DM and WFS rotate/translate
   - Consistent handling of complex geometric configurations

**Automatic Workflow Selection**

The workflow is automatically selected based on system geometry analyzing transformation parameters.

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

**Memory Management**
   - Automatic data transfer between CPU and GPU
   - Graceful fallback to CPU on memory errors
   - Efficient caching of GPU arrays

**Precision Control**
   - ``precision=0``: Double precision (float64)
   - ``precision=1``: Single precision (float32, ~2× faster on GPU)

Configuration Files
===================

SynIM supports both YAML and PRO (IDL-style) configuration files for defining AO system parameters.
IDL-style parameter files are supported for compatibility with `PASSATA <https://arxiv.org/abs/1607.07624>`_.

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
   
   # Or PRO file
   params = parse_params_file('config.pro')

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

SynIM follows SPECULA's directory structure for seamless integration. The ``ParamsManager`` automatically creates and manages these directories based on the ``root_dir`` parameter in your configuration:

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


**Note:** When using SynIM with SPECULA, both tools can share the same ``root_dir`` and ``ifunc/`` directory. SPECULA uses additional directories (``im/``, ``rec/``) which can coexist alongside SynIM's directories.

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