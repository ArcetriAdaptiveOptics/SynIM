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

Test your installation:

.. code-block:: python

   import synim
   import numpy as np
   
   # Test CPU installation
   synim.init(device_idx=-1, precision=1)
   print("✓ CPU backend initialized")
   
   # Test GPU installation (if available)
   try:
       synim.init(device_idx=0, precision=1)
       print("✓ GPU backend initialized")
   except Exception:
       print("⚠ GPU not available (CPU only)")
   
   # Test basic functionality
   pup_mask = np.ones((128, 128))
   print("✓ SynIM installation successful!")


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