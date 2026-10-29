.. _installation:

Installation
============

SynIM requires Python 3.8 or higher and we strongly recommend using conda for package management.

Prerequisites
-------------

**System Requirements:**
   * Python 3.8 or higher
   * Git (for repository cloning)
   * CUDA-compatible GPU (optional, for acceleration)

**Recommended Setup:**
   * Anaconda or Miniconda
   * 8GB+ RAM
   * 8GB+ GPU memory (if using GPU acceleration)

Step 1: Create Conda Environment
---------------------------------

Create a dedicated conda environment for SynIM (here with python 3.11):

.. code-block:: bash

   # Create environment with Python 3.11
   conda create --name synim python=3.11
   
   # Activate the environment
   conda activate synim

Step 2: GPU Support (Optional but Recommended)
----------------------------------------------

If you have a CUDA-compatible GPU and want to benefit from GPU acceleration, install CuPy:

.. code-block:: bash

   # Install CuPy for GPU acceleration
   conda install -c conda-forge cupy

**GPU Benefits:**
   * 2-10× faster matrix computations
   * Enables processing of larger systems (higher resolution pupils)

**Without GPU:**
   SynIM will automatically fall back to CPU computation using NumPy. Performance will be slower but all functionality remains available.

Step 3: Install SynIM
---------------------

**A) Install from PyPI (recommended for most users):**

.. code-block:: bash

   # Install SynIM from PyPI
   pip install synim

**B) Install from source (for development or latest features):**

Clone the SynIM repository from GitHub:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/ArcetriAdaptiveOptics/SynIM.git
   
   # Navigate to the directory
   cd SynIM

Then install in development mode:

.. code-block:: bash

   pip install -e .

This installs SynIM in "editable" mode, allowing you to modify the code and see changes immediately.

**Required Dependencies:**
All required dependencies will be installed automatically, including:

* **numpy**: Numerical computing foundation
* **scipy**: Scientific computing and interpolation
* **matplotlib**: Plotting and visualization

**Additional Requirements for ParamsManager:**

If you plan to use the ``ParamsManager`` class for interaction matrix computation, 
you must install SPECULA separately:

.. code-block:: bash

   pip install specula

Or install from source:

.. code-block:: bash

   git clone https://github.com/ArcetriAdaptiveOptics/SPECULA.git
   cd SPECULA
   pip install -e .

.. note::
   ``ParamsManager`` provides high-level utilities for computing interaction matrices
   in the context of adaptive optics simulations. The core SynIM functionality
   (``synim.synim``, ``synim.synpm``, ``synim.utils``) does not require SPECULA.

.. warning::
   ``ParamsManager`` will be moved to the SPECULA package in a future release (v2.0.0).
   Users are encouraged to migrate to ``specula.utils.ParamsManager`` when available.


**Optional Dependencies:**

* **cupy**: GPU acceleration (installed in Step 2)
* **pyyaml**: YAML configuration file support (recommended)

Verification
------------

Test your installation by computing the interaction matrix of a tilt, whose
slopes are known analytically: a tilt of 1 nm at the edge of a pupil of
diameter :math:`D` gives an x slope of :math:`10^{-9}/(D/2)` rad on every
illuminated subaperture (expressed in units of half the subaperture field of
view) and a zero y slope.

.. code-block:: python

   import numpy as np
   import synim
   synim.init(device_idx=-1, precision=1)   # device_idx=0 to test the GPU
   from synim.synim import interaction_matrix

   # Circular pupil (64 px, D = 8 m) and one mode: x tilt of 1 nm at the pupil edge
   n, D = 64, 8.0
   y, x = (np.mgrid[0:n, 0:n] - (n - 1) / 2) / (n / 2)
   pup = (x**2 + y**2 <= 1).astype(np.float32)
   tilt = (x * pup)[..., None].astype(np.float32)

   im = synim.cpuArray(interaction_matrix(
       pup_diam_m=D, pup_mask=pup, dm_array=tilt, dm_mask=pup,
       dm_height=0.0, dm_rotation=0.0, wfs_nsubaps=8, wfs_fov_arcsec=4.0,
       gs_pol_coo=(0.0, 0.0), gs_height=np.inf, wfs_rotation=0.0,
       wfs_translation=(0.0, 0.0), wfs_mag_global=1.0))

   # Expected x slope: (1 nm / (D/2)) rad, in units of half the subaperture FoV (2 arcsec)
   expected = 1e-9 / (D / 2) / np.deg2rad(1 / 3600) / 2.0
   sx, sy = np.split(im[:, 0], 2)
   lit = sx != 0
   assert lit.sum() > 0 and np.allclose(sx[lit], expected, rtol=1e-3) and np.allclose(sy, 0)
   print("✓ SynIM installation successful!")

Once a SynIM submodule has been imported, ``synim.init`` cannot switch to a
different device or precision: run the snippet in a new Python session to test
another device.


**Environment Variables:**

You can control GPU usage with environment variables:

.. code-block:: bash

   # Disable GPU globally (force CPU)
   export SYNIM_DISABLE_GPU=TRUE
   
   # Select specific GPU device
   export CUDA_VISIBLE_DEVICES=0
   
   # Enable CUDA memory pool for better performance
   export CUPY_CACHE_SAVE_CUDA_SOURCE=1

Environment Management
----------------------

**Useful conda commands:**

.. code-block:: bash

   # List environments
   conda env list
   
   # Activate SynIM environment
   conda activate synim
   
   # Deactivate environment
   conda deactivate
   
   # Update all packages
   conda update --all
   
   # Remove environment (if needed)
   conda env remove --name synim

**Updating SynIM:**

.. code-block:: bash

   # Navigate to SynIM directory
   cd SynIM
   
   # Pull latest changes
   git pull origin main
   
   # Reinstall if needed
   pip install -e .