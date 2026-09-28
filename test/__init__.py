# SynIM must be initialized before importing its submodules (see synim.init).
# The tests run on CPU, in single precision.
import synim

synim.init(device_idx=-1, precision=1)
