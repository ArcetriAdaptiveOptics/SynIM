# SynIM must be initialized before importing its submodules (see synim.init).
# The tests run in single precision, on CPU by default. To run them on a GPU:
#     SYNIM_TEST_DEVICE_IDX=0 python -m unittest discover -t . -s test
import os

import synim

synim.init(device_idx=int(os.environ.get('SYNIM_TEST_DEVICE_IDX', '-1')), precision=1)


def on_backend(func):
    """
    Wrap a SynIM function for tests written with numpy arrays: the numpy
    array arguments are moved to the SynIM backend (cupy on GPU) and the
    array results (also in tuples, lists and dicts) are returned as numpy
    arrays. On CPU it only calls func.
    """
    import functools

    import numpy as np

    def to_backend(value):
        return synim.xp.asarray(value) if isinstance(value, np.ndarray) else value

    def to_numpy(value):
        if isinstance(value, (tuple, list)):
            return type(value)(to_numpy(v) for v in value)
        if isinstance(value, dict):
            return {k: to_numpy(v) for k, v in value.items()}
        if hasattr(value, '__cuda_array_interface__') or isinstance(value, np.ndarray):
            return synim.cpuArray(value)
        return value

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        args = [to_backend(a) for a in args]
        kwargs = {k: to_backend(v) for k, v in kwargs.items()}
        return to_numpy(func(*args, **kwargs))

    return wrapper
