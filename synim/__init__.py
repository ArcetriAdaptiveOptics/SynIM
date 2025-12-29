import numpy as np
import os

# Global variables for array library configuration
xp = None
cp = None
gpuEnabled = False
default_target_device_idx = None
default_target_device = None

# *** Add scipy modules as global ***
affine_transform = None
binary_dilation = None
rotate = None
shift = None
zoom = None

# *** Precision management ***
global_precision = None
float_dtype = None
complex_dtype = None
cpu_float_dtype_list = [np.float64, np.float32]
cpu_complex_dtype_list = [np.complex128, np.complex64]
gpu_float_dtype_list = cpu_float_dtype_list
gpu_complex_dtype_list = cpu_complex_dtype_list

def init(device_idx=-1, precision=1):
    """
    Initialize SynIM with numpy or cupy backend.
    
    Args:
        device_idx (int): GPU device index (-1 for CPU, >=0 for GPU)
        precision (int): 0 for double precision, 1 for single precision
    
    Returns:
        None
    """
    global xp, cp, gpuEnabled, default_target_device_idx, default_target_device
    global global_precision, float_dtype, complex_dtype
    global gpu_float_dtype_list, gpu_complex_dtype_list
    # *** Declare scipy globals ***
    global affine_transform, binary_dilation

    default_target_device_idx = device_idx
    global_precision = precision

    # Check if GPU is disabled by environment variable
    systemDisable = os.environ.get('SYNIM_DISABLE_GPU', 'FALSE')

    if systemDisable == 'FALSE':
        try:
            import cupy as cp_module
            print(f"CuPy import successful. Installed version: {cp_module.__version__}")
            gpuEnabled = True
            cp = cp_module
            gpu_float_dtype_list = [cp.float64, cp.float32]
            gpu_complex_dtype_list = [cp.complex128, cp.complex64]
        except ImportError:
            print("CuPy import failed. SynIM will fall back to CPU use.")
            cp = None
            xp = np
            default_target_device_idx = -1
    else:
        print("Environment variable SYNIM_DISABLE_GPU prevents using the GPU.")
        cp = None
        xp = np
        default_target_device_idx = -1

    # ==================== SET ARRAY LIBRARY ====================
    if default_target_device_idx >= 0:
        xp = cp
        float_dtype_list = [cp.float64, cp.float32]
        complex_dtype_list = [cp.complex128, cp.complex64]
        default_target_device = cp.cuda.Device(default_target_device_idx)
        default_target_device.use()
        print(f'Default device is GPU number {default_target_device_idx}')

        # Try to import cupyx.scipy for GPU acceleration
        try:
            from cupyx.scipy.ndimage import affine_transform as cupy_affine
            from cupyx.scipy.ndimage import binary_dilation as cupy_dilation
            from cupyx.scipy.ndimage import rotate as cupy_rotate
            from cupyx.scipy.ndimage import shift as cupy_shift
            from cupyx.scipy.ndimage import zoom as cupy_zoom
            affine_transform = cupy_affine
            binary_dilation = cupy_dilation
            rotate = cupy_rotate
            shift = cupy_shift
            zoom = cupy_zoom
            print('✓ Using cupyx.scipy.ndimage (GPU-accelerated transforms)')
        except ImportError:
            print('⚠️  cupyx.scipy.ndimage not available, falling back to scipy (CPU)')
            from scipy.ndimage import affine_transform as cpu_affine
            from scipy.ndimage import binary_dilation as cpu_dilation
            from scipy.ndimage import rotate as cpu_rotate
            from scipy.ndimage import shift as cpu_shift
            from scipy.ndimage import zoom as cpu_zoom
            affine_transform = cpu_affine
            binary_dilation = cpu_dilation
            rotate = cpu_rotate
            shift = cpu_shift
            zoom = cpu_zoom
            print('✓ Using scipy.ndimage (CPU)')

    else:
        print('Default device is CPU')
        xp = np
        float_dtype_list = [np.float64, np.float32]
        complex_dtype_list = [np.complex128, np.complex64]

        # *** Use scipy for CPU ***
        from scipy.ndimage import affine_transform as cpu_affine
        from scipy.ndimage import binary_dilation as cpu_dilation
        from scipy.ndimage import rotate as cpu_rotate
        from scipy.ndimage import shift as cpu_shift
        from scipy.ndimage import zoom as cpu_zoom
        affine_transform = cpu_affine
        binary_dilation = cpu_dilation
        rotate = cpu_rotate
        shift = cpu_shift
        zoom = cpu_zoom
        print('✓ Using scipy.ndimage (CPU)')

    float_dtype = float_dtype_list[global_precision]
    complex_dtype = complex_dtype_list[global_precision]

    if precision == 0:
        print('Using double precision (float64)')
    else:
        print('Using single precision (float32)')


def to_xp(target_xp, v, dtype=None, force_copy=False):
    """
    Convert array to target array library (numpy or cupy).
    Optimized to avoid unnecessary copies.
    
    Args:
        target_xp: Target array library (np or cp)
        v: Input array or array-like
        dtype: Optional target dtype
        force_copy: Force copy even if already correct type
    
    Returns:
        Array in target library
    """
    # *** OPTIMIZED: Fast path for already-correct arrays ***
    if not force_copy and dtype is None:
        if target_xp is np and isinstance(v, np.ndarray):
            return v
        if target_xp is cp and cp is not None and isinstance(v, cp.ndarray):
            return v

    # Convert between libraries if needed
    if target_xp is cp:
        if cp is None:
            raise RuntimeError("CuPy not available, cannot convert to GPU array")
        if isinstance(v, cp.ndarray) and not force_copy:
            retval = v
        else:
            retval = cp.array(v)
    else:
        if cp is not None and isinstance(v, cp.ndarray):
            retval = v.get()  # GPU -> CPU
        elif isinstance(v, np.ndarray) and not force_copy:
            retval = v
        else:
            retval = np.array(v)

    # Apply dtype conversion if specified
    if dtype is not None and retval.dtype != dtype:
        return retval.astype(dtype, copy=False)

    return retval


def cpuArray(v, dtype=None, force_copy=False):
    """
    Ensure array is on CPU (numpy).
    
    Args:
        v: Input array
        dtype: Optional target dtype
        force_copy: Force copy
    
    Returns:
        Numpy array
    """
    return to_xp(np, v, dtype=dtype, force_copy=force_copy)


# Initialize with CPU by default, single precision
init(device_idx=-1, precision=1)
