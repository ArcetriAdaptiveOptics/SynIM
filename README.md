# SynIM

**Synthetic Interaction Matrix generator for Adaptive Optics systems**

[![Documentation Status](https://readthedocs.org/projects/synim/badge/?version=latest)](https://synim.readthedocs.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

SynIM is a Python package for computing synthetic interaction matrices, projection matrices, and covariance matrices for adaptive optics (AO) systems. It supports both Single Conjugate AO (SCAO), Laser Tomography AO (LTAO), Ground Layer AO (GLAO) and Multi-Conjugate AO (MCAO) configurations with Shack-Hartmann sensors.
It also supports GPU acceleration. Some of its functionalities are provided by [SPECULA](https://github.com/ArcetriAdaptiveOptics/SPECULA).

## Key Features

- 🔧 **Interaction matrices** for DM-WFS combinations (SCAO/MCAO)
- 📊 **Projection matrices** for DM-Layer tomography (MCAO/LTAO)
- 📈 **Covariance matrices** for MMSE reconstructor optimization
- ⚙️ **YAML/PRO configuration** with automatic parameter parsing
- 💾 **SPECULA-compatible** FITS format for data exchange
- 🎯 **Smart caching** to minimize redundant computations
- 🚀 **GPU acceleration** via CuPy for high-performance computation
- 🧭 **Mis-registration geometry** (`synim.registration`): a geometric (not pixel-based) model of the relation between global (WFS/DM/guide-star) and local (per-pair) mis-registration in MCAO/GLAO systems, its inversion, and sensitivity/noise analysis - see the [documentation](https://synim.readthedocs.io/en/latest/registration.html)

## Performance

From version 1.2.1, used in the SynIM paper, to version 2.1, the calibration has become significantly faster and uses much less memory, both on CPU and on GPU. For example, the 72 MORFEO interaction matrices of the benchmark in the paper (numerical derivatives, NVIDIA L40S) take 93 s with version 2.1 instead of 265.6 s with version 1.2.1, with the same input and output: influence functions read from a local disk with no file cache, and matrices saved to FITS files. Most of the remaining time is spent reading the influence functions (about 41 GB): with the files already in the file cache and the matrices kept in GPU memory, the computation takes 13.5 s. The main changes are: interaction matrices computed as a sparse linear operator applied to the DM modes, so that the extra memory no longer grows with the number of modes; an MMSE reconstructor that exploits the sparsity of the noise covariance and solves a linear system instead of computing an explicit inverse; the whole reconstruction chain on the GPU; and projection matrices processed one optical source at a time. The measurements are in the [Performance](https://synim.readthedocs.io/en/latest/performance.html) page of the documentation.

## Documentation

Full documentation available at: **[synim.readthedocs.io](https://synim.readthedocs.io)**

## Slope Computation Methods

SynIM supports two slope computation methods for interaction matrices:

- `derivatives` (default): compute numerical derivatives on the transformed phase
- `telsum` (optional): compute telescoping sum from phase differences in each subaperture

The default behavior is `derivatives`. To enable telescoping sum, pass `slope_method='telsum'` when calling the low-level interaction-matrix APIs.

## Requirements

- Python ≥ 3.8
- numpy, scipy, matplotlib
- [specula](https://github.com/ArcetriAdaptiveOptics/SPECULA) - AO simulation framework
- [cupy](https://cupy.dev/) - GPU acceleration (optional, requires CUDA)

## Contributing to SynIM
To contribute to SynIM, follow these steps:

1. Fork this repository.
2. Create a branch: `git checkout -b <branch_name>`
3. Make your changes and **add tests for the new functionality.**
4. Commit your changes: `git commit -m '<commit_message>'`
5. Push to the branch: `git push`
6. Create the pull request.

We require tests for all new features to ensure the stability of the project.

## Citation

If you use SynIM in your research, please cite:

```bibtex
@misc{agapito2026synim,
  title         = {SynIM: a high-performance GPU-accelerated Python library for
                   synthetic interaction and tomographic reconstruction matrices
                   in next-generation adaptive optics},
  author        = {Agapito, Guido and Rossi, Fabio and Puglisi, Alfio},
  year          = {2026},
  eprint        = {2606.07759},
  archivePrefix = {arXiv},
  primaryClass  = {astro-ph.IM},
  url           = {https://arxiv.org/abs/2606.07759}
}
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Authors

- **Guido Agapito** - INAF - Osservatorio Astrofisico di Arcetri

## Acknowledgments

It was built as part of the [MORFEO project](https://elt.eso.org/instrument/MORFEO/).