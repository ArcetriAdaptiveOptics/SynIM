Performance
===========

Since version 1.2.1, used in the SynIM paper, most steps of an MCAO calibration
are faster and use less memory. This page compares version 2.1 with the
benchmark of the paper, and then collects the measurements made on the pull
requests (PRs) that introduced each change. The test case is MORFEO, the MCAO system of the
ELT: 6 LGS Shack-Hartmann WFSs with 68×68 subapertures, 480×480 pupil pixels,
3 DMs with 6564 modes and 9 reconstruction layers with 35042 modes, in single
precision unless noted.

Comparison with the SynIM paper
-------------------------------

The benchmark of the SynIM paper computes the 72 interaction matrices of the
MORFEO LGS WFSs (6 LGS × 3 DMs and 9 layers, numerical derivatives) with
``ParamsManager.compute_interaction_matrices``, on the same machine (NVIDIA L40S)
for both versions:

.. list-table::
   :header-rows: 1
   :widths: 22 33 25 20

   * - Version
     - Influence functions (~41 GB)
     - Interaction matrices
     - Time
   * - 1.2.1 (paper)
     - read from disk, no file cache
     - saved to FITS files
     - 265.6 s
   * - 2.1
     - read from a local disk, no file cache
     - saved to FITS files
     - 93 s
   * - 2.1
     - read from a local disk, no file cache
     - kept in GPU memory
     - 79.5 s
   * - 2.1
     - already in the file cache
     - kept in GPU memory
     - 13.5 s

The first two rows do the same work. With version 2.1 most of the time is
input/output: reading the influence functions takes about 66 s and saving the
matrices about 13 s, while the computation takes about 13.5 s. The peak GPU
memory in use is 14.2 GB with the matrices saved to files and 19.5 GB with the
matrices kept in GPU memory (``in_memory=True``).

Measurements of single pull requests
------------------------------------

The rest of this page collects the measurements made on the pull requests that
introduced each change.

.. note::

   Each measurement compares a PR with the ``main`` branch just before it, on the
   machine used for that PR (NVIDIA L40S, NVIDIA RTX PRO 5000, or a 64-core CPU).
   The numbers of different rows are not measured on the same baseline, so they
   should not be multiplied together. Details are in the release notes and in the
   linked PRs.

Interaction matrices
--------------------

From 2.1, an interaction matrix is a sparse linear operator applied to the cube
of DM modes (#95): the transformation, the derivatives (or the telescoping sum)
and the integration over the subapertures are a single sparse matrix that does
not depend on the number of modes. The explicit step-by-step implementation is
kept as a reference (``method='reference'``).

.. list-table::
   :header-rows: 1
   :widths: 50 20 20 10

   * - Measurement
     - Before
     - After
     - PR
   * - Single IM, GPU, MORFEO size, derivatives / telescoping sum (L40S)
     - 330 / 212 ms
     - 79-84 / 66 ms
     - #95
   * - 72 LGS IMs, GPU, derivatives / telescoping sum (L40S)
     - 19.7 / 13.1 s
     - 6.1 / 5.4 s
     - #95
   * - Extra GPU memory for an IM, 5000 modes (derivatives)
     - 3.25 mode cubes (~15 GB)
     - 0.02 mode cubes
     - #95
   * - Single IM, CPU, 480×480 pixels, 40×40 subapertures, 500 modes
     - 6.2 s
     - 0.11 s
     - #95
   * - Single IM, CPU, derivatives (2.0)
     - 175 s
     - 50 s
     - #90
   * - GPU peak memory, single IM, telescoping sum (2.0)
     - 43.8 GB
     - 30.2 GB
     - #89

Assembly of the full interaction matrices (slopes filter and mode selection on
the compute device):

.. list-table::
   :header-rows: 1
   :widths: 50 20 20 10

   * - Measurement
     - Before
     - After
     - PR
   * - ``assemble_interaction_matrices``, layers (37736 × 35042)
     - 24.4 s
     - 9.4 s
     - #104
   * - ``assemble_interaction_matrices``, layers, kept on the device
     - 5.07 s
     - 0.30 s
     - #111
   * - Compute and assemble without IM files (``in_memory=True``), total
     - 55.6 s
     - 43.3 s
     - #107

Covariances and noise
---------------------

.. list-table::
   :header-rows: 1
   :widths: 50 20 20 10

   * - Measurement
     - Before
     - After
     - PR
   * - ``dm3d_to_2d``, 9 layers (CPU)
     - 159 s
     - 37 s
     - #100
   * - Inversion of the covariance blocks (Cholesky on the device instead of ``pinv``)
     - 38.6 s
     - 3.1 s
     - #114
   * - ``build_noise_covariance``, elongated LGS model (CPU)
     - 21.7 s
     - 1.8 s
     - #102
   * - Extra host memory of the illumination scaling
     - 5.7 GB
     - none
     - #102

MMSE reconstructor
------------------

The MMSE reconstructor uses the sparsity of the noise precision matrix (at most
2 non-zeros per row for the elongated LGS model), a Cholesky or LU solve instead
of an explicit inverse (#105), only the lower tiles of the symmetric matrix
H = Aᵀ Cz⁻¹ A (#116), and copies the covariance inverses to the device only when
they are needed (#117). The whole chain can run on the GPU (#106, #111, #112, #113).

.. list-table::
   :header-rows: 1
   :widths: 50 20 20 10

   * - Measurement
     - Before
     - After
     - PR
   * - ``compute_mmse_reconstructor``, GPU (RTX PRO 5000)
     - 14.6 s
     - 8.2 s
     - #105
   * - ``compute_mmse_reconstructor``, CPU, symmetric Cz⁻¹
     - 169 s
     - 85 s
     - #105
   * - GPU memory above the inputs
     - 20.4 GB
     - 15.5 GB
     - #105
   * - H = Aᵀ Cz⁻¹ A on GPU
     - 2.40 s
     - 1.38 s
     - #116
   * - ``compute_mmse_reconstructor`` in ``compute_tomographic_reconstructor``, CPU vs GPU
     - 73.2 s
     - 10.5 s
     - #113
   * - Host memory peak of the reconstructor chain
     - 34.8 GB
     - 14.2 GB
     - #111, #112

Tomographic projection matrix
-----------------------------

The projection matrices of the optical sources are now read one source at a
time, and only the two weighted sums needed by the projection are accumulated,
in double precision (#110).

.. list-table::
   :header-rows: 1
   :widths: 50 20 20 10

   * - Measurement
     - Before
     - After
     - PR
   * - Peak host memory (17 optical sources)
     - ~100 GB
     - 9.3 GB
     - #110
   * - Relative error of the projection matrix vs a double precision computation, CPU
     - 6.4e-4
     - 2.5e-8
     - #110

The time of this step is now dominated by reading the projection matrix files.
