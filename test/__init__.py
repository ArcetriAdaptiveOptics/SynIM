# SynIM must be initialized before importing its submodules (see synim.init).
# The tests run in single precision, on CPU by default. To run them on a GPU:
#     SYNIM_TEST_DEVICE_IDX=0 python -m unittest discover -t . -s test
import os

import synim

synim.init(device_idx=int(os.environ.get('SYNIM_TEST_DEVICE_IDX', '-1')), precision=1)
