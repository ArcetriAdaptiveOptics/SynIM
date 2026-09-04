.. _registration:

Mis-registration Geometry
==========================

The ``synim.registration`` subpackage models the geometric relation
between the **global** mis-registration of a wide-field AO (WFAO)
system - the position of each wavefront sensor (WFS) and each
deformable mirror (DM) relative to the pupil, and of each guide star
relative to its assumed position - and the **local** mis-registration
of each WFS-DM pair - the shift, rotation, magnification and
anamorphosis a local estimator measures between one sub-aperture grid
and one actuator grid.

It follows the geometry of Agapito, Plantet & Heritier (2024),
*SPRINT for WFAO systems*, Proc. SPIE 13097
(`arXiv:2406.15336 <https://doi.org/10.48550/arXiv.2406.15336>`_),
generalizing the SCAO registration method of Heritier et al. (2021) to
systems with several WFSs and DMs (MCAO/GLAO).

.. important::
   This subpackage does **not** estimate local mis-registration. It
   takes local measurements as given - produced elsewhere, for example
   by SPRINT (already implemented in `SPECULA
   <https://github.com/ArcetriAdaptiveOptics/SPECULA>`_) - and provides
   the forward map from global parameters to local measurements, its
   inversion, and tools to study its conditioning and noise sensitivity.
   It does not touch interaction matrices or on-sky calibration.

Local and global parameters
----------------------------

Local mis-registration, per WFS-DM pair:

- X and Y shift between sub-apertures and actuators
- rotation between sub-apertures and actuators
- magnification and anamorphosis between sub-apertures and actuators

Global mis-registration, one set per element:

- shift, rotation and magnification of a WFS relative to the pupil
- shift, rotation, magnification and anamorphosis of a DM relative to
  the pupil
- shift and altitude error of a guide star relative to its assumed
  position (altitude error applies to laser guide stars only)

A single DM at the pupil (for example an adaptive secondary mirror)
makes the two levels coincide for that DM: see `System.local_transform`
below.

Building a system
------------------

A system is a set of ``WFS`` objects (each with its own ``GuideStar``)
and a set of ``DM`` objects, collected in a ``System``:

.. code-block:: python

   from synim.registration.model import GuideStar, DM, WFS, System

   system = System(pixel_pitch=0.0802)  # meters per shift unit

   gs = GuideStar(name="lgs1", position=(45.0, 22.5), height=90000.0)
   system.add_wfs(WFS(name="wfs1", guide_star=gs, rotation=6.2))

   system.add_dm(DM(name="dm1", height=600.0))
   system.add_dm(DM(name="dm2", height=6500.0))

Every shift is expressed in the same unit as ``pixel_pitch`` (matching
the ``xShiftPhInPixel``/``yShiftPhInPixel`` convention of a SPECULA
YAML file); angles are in degrees, guide star positions in arcsec.
Nominal (zero mis-registration) values are the default, so a freshly
built system already represents a perfectly aligned instrument.

The forward model: from global to local
------------------------------------------

``System.local_params`` composes a WFS's own mis-registration, the
guide-star parallax at a DM's altitude, and that DM's own
mis-registration into the local shift/rotation/magnification/
anamorphosis of one pair:

.. code-block:: python

   shift, rotation, magnification, anamorphosis = system.local_params("wfs1", "dm1")

``System.pairs()`` lists every WFS-DM pair; restrict it to a subset
(for example the LGS WFSs only) when some sensors carry no useful
shape information - see the :ref:`MORFEO tutorial
<registration_morfeo_tutorial>`.

Loading a system from a configuration file
---------------------------------------------

``System.from_params_manager`` builds a system directly from an already
parsed SPECULA/SynIM YAML configuration, reusing
:func:`synim.params_utils.extract_wfs_list` and related helpers:

.. code-block:: python

   from synim.params_utils import parse_params_file
   from synim.registration.model import System

   config = parse_params_file("params_mcao.yml")
   system = System.from_params_manager(config)

Only the WFS fields already used elsewhere in SynIM (``rotation``,
``xShiftPhInPixel``/``yShiftPhInPixel``, ``magnification``,
``anamorph45``) are read this way; a DM's own shift/magnification and a
guide star's position/altitude error have no dedicated key in that
format and default to zero - set them directly on the returned objects,
or through the reconstruction tools below.

Inverting: global parameters from local measurements
--------------------------------------------------------

Given a set of local measurements (from SPRINT, or - as in the example
above - injected for testing), :mod:`synim.registration.reconstruction`
recovers the global parameters behind them by Gauss-Newton iteration:

.. code-block:: python

   from synim.registration.reconstruction import ParameterSpec, gauss_newton_invert

   specs = [
       ParameterSpec("dm", "dm2", "shift", 0),
       ParameterSpec("dm", "dm2", "rotation"),
       ParameterSpec("wfs", "wfs1", "shift", 0),
   ]
   alpha_hat, system_hat, history = gauss_newton_invert(
       system, specs, system.pairs(), local_measurements)

A ``ParameterSpec`` names one scalar unknown (an object, a field, and -
for a shift - which component); ``history`` is the residual norm at
each iteration, for a convergence check. At each step, a sensitivity
matrix (the same role as an interaction matrix, mapping global
parameters to local measurements) is built by finite differences and
inverted by pseudo-inverse.

Sensitivity and degeneracy
----------------------------

Not every choice of unknowns is well posed. Two cases the paper points
out, and that :mod:`synim.registration.analysis` can detect directly
from the sensitivity matrix built above:

- a common shift applied to every WFS together with every DM
  (each scaled by its own parallax factor) leaves every local
  measurement unchanged - an exact gauge freedom, since only the
  *relative* WFS-DM alignment is observable;
- a guide star's position error and that WFS's own shift affect a
  single WFS-DM pair identically, and can only be told apart using
  several DMs at different altitudes.

.. code-block:: python

   from synim.registration.analysis import analyze

   report = analyze(system, specs, system.pairs(), sigma=0.02)
   print(report["condition_number"], report["modes"][-1])

``report["modes"]`` lists, from best to worst determined, which
``ParameterSpec`` entries dominate each direction; passing ``sigma``
(the local-measurement noise, one value or one per degree of freedom)
also returns the resulting global-parameter uncertainty, both from the
pseudo-inverse directly (``analytic_std``) and from a Monte Carlo
cross-check (``montecarlo_std``, via `monte_carlo_noise_propagation`).

.. important::
   That Monte Carlo check reconstructs with a single linear step,
   linearized exactly at the true values used to generate ``report``'s
   sensitivity matrix - the same linear model behind ``analytic_std``.
   Its expectation therefore equals those true values by construction,
   so it can validate the analytic UNCERTAINTY, but it is unbiased by
   construction and cannot reveal an estimation BIAS. To check for a
   real bias, run `monte_carlo_gauss_newton` instead: it repeats the
   actual iterative estimator (`gauss_newton_invert`) from a given
   starting system - typically the nominal one, as an operational
   estimator would, not knowing the true values in advance. Compare its
   returned mean to the true values, in units of ``std / sqrt(n_trials)``
   (the standard error of that mean), to tell a real bias from Monte
   Carlo sampling noise.

Visualization
--------------

:mod:`synim.registration.viz` provides an altitude schematic, a
mis-registration table, and a top-down footprint view per DM, all
colour-matched and combinable through ``plot_system_overview``. See the
:ref:`MORFEO tutorial <registration_morfeo_tutorial>` for example
output.

See Also
--------

- :ref:`registration_morfeo_tutorial` - a worked example on a MORFEO-like geometry
- :doc:`api/registration` - full API reference
- Agapito, Plantet & Heritier (2024), *SPRINT for WFAO systems*,
  `arXiv:2406.15336 <https://doi.org/10.48550/arXiv.2406.15336>`_
