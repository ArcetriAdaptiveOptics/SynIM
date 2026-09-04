.. _registration_morfeo_tutorial:

=========================================================
Mis-registration Geometry Tutorial: A MORFEO LGS Example
=========================================================

This tutorial builds a MORFEO-like geometry with ``synim.registration``
(see :ref:`registration` for the concepts) and walks through the forward
model, its inversion, a degeneracy check, and noise propagation. It only
needs ``synim`` and runs on its own - no SPECULA configuration file
required.

.. note::
   For a companion example built from a real MORFEO SPECULA
   configuration - including the extra visualizations available through
   SPECULA ``plot_utils`` - see
   ``SPECULA_scripts/morfeo/mis_registration_geometry/`` alongside this
   repository.

1. Building the geometry
=========================

MORFEO senses six laser guide stars (LGS) on an evenly spaced ring at
45 arcsec, and corrects with three DMs at 600 m, 6500 m and 17500 m.
This tutorial keeps to the LGS path: a 1x1 low-order sensor, like
MORFEO natural guide star channel, has no sub-aperture grid to
register, so it carries no useful shift/rotation/magnification of its
own.

.. code-block:: python

   import numpy as np
   from synim.registration.model import GuideStar, DM, WFS, System

   pixel_pitch = 0.0802  # m per shift unit, matching a 480-pixel, ~38.5 m pupil
   pupil_diameter = 38.5  # m

   system = System(pixel_pitch=pixel_pitch)
   for i, angle_deg in enumerate([30, 90, 150, 210, 270, 330], start=1):
       angle = np.radians(angle_deg)
       position = (45.0 * np.cos(angle), 45.0 * np.sin(angle))
       gs = GuideStar(name=f"gs_lgs{i}", position=position, height=90000.0)
       system.add_wfs(WFS(name=f"lgs{i}", guide_star=gs))

   for name, height in [("dm1", 600.0), ("dm2", 6500.0), ("dm3", 17500.0)]:
       system.add_dm(DM(name=name, height=height))

   pairs = system.pairs()  # all 18 WFS-DM pairs

Every field just built is at its default (zero shift, zero rotation,
unit magnification): this ``system`` is a perfectly aligned MORFEO.

2. A first look
================

.. code-block:: python

   from synim.registration.viz import plot_system_overview

   fig, axes = plot_system_overview(system, pupil_diameter=pupil_diameter)
   fig.savefig("morfeo_lgs_overview.png", bbox_inches="tight")

This produces an altitude schematic, a parameter table, and one
top-down footprint panel per DM - at 600 m the six LGS footprints
nearly coincide, and they separate increasingly at 6500 m and 17500 m
as the cone effect narrows each footprint and the parallax shifts it.

3. Forward model and inversion
================================

Inject a mis-registration on two DMs and one WFS, compute the local
parameters it produces (standing in for a local estimator
measurement), and recover it by Gauss-Newton inversion:

.. code-block:: python

   from synim.registration.reconstruction import (
       ParameterSpec, apply_alpha, local_params_vector, gauss_newton_invert, LOCAL_DOF,
   )

   specs = [
       ParameterSpec("dm", "dm2", "shift", 0),
       ParameterSpec("dm", "dm2", "shift", 1),
       ParameterSpec("dm", "dm3", "rotation"),
       ParameterSpec("wfs", "lgs1", "shift", 0),
   ]
   true_alpha = np.array([0.08, -0.05, -0.10, 0.06]) / [pixel_pitch, pixel_pitch, 1.0, pixel_pitch]

   true_system = apply_alpha(system, specs, true_alpha)
   measurements = local_params_vector(true_system, pairs, dof=LOCAL_DOF)

   alpha_hat, system_hat, history = gauss_newton_invert(system, specs, pairs, measurements)
   print(np.abs(alpha_hat - true_alpha).max())  # near machine precision

``history`` is the residual norm at each iteration; it should drop to
machine precision within a handful of steps, since this is a noise-free
inversion (Sec. 4 of the paper does the same check for MAVIS).

4. A degeneracy to watch for
==============================

Not every choice of unknowns is observable. Shifting every WFS and
every DM *together* is a gauge freedom: only the relative alignment
between a WFS and a DM matters, so a common shift of the whole system
leaves every local measurement unchanged.

.. code-block:: python

   from synim.registration.reconstruction import jacobian
   from synim.registration.analysis import condition_number, describe_modes, svd_of_jacobian

   wide_specs = (
       [ParameterSpec("wfs", f"lgs{i}", "shift", axis) for i in range(1, 7) for axis in (0, 1)]
       + [ParameterSpec("dm", name, "shift", axis) for name in ("dm1", "dm2", "dm3") for axis in (0, 1)]
   )
   Lambda = jacobian(system, wide_specs, pairs, dof=("shift_x", "shift_y"))
   print(condition_number(Lambda))  # ~1e11: an exact degeneracy

   _, _, Vt = svd_of_jacobian(Lambda)
   for label, coeff in describe_modes(Vt, wide_specs, top_n=9)[-1]:
       print(f"{label:20s} {coeff:+.3f}")

The worst direction assigns nearly the same coefficient to every WFS
shift and to every DM shift: exactly that common-shift combination. The
scenario in Step 3 avoids it by leaving most WFSs and one DM fixed - a
real estimator must restrict its unknowns the same way, or accept that
an absolute offset of the whole system is unrecoverable (which is
harmless: only the relative registration affects AO correction).

5. Propagating measurement noise
===================================

Given a per-degree-of-freedom measurement noise, ``analyze`` returns
the resulting uncertainty on the global parameters, both from the
pseudo-inverse directly and from a Monte Carlo cross-check:

.. code-block:: python

   from synim.registration.analysis import analyze

   sigma_per_dof = [0.02 / pixel_pitch, 0.02 / pixel_pitch, 0.03, 0.005, 0.005]
   sigma = np.tile(sigma_per_dof, len(pairs))

   report = analyze(true_system, specs, pairs, sigma=sigma, n_trials=2000)
   print(report["analytic_std"])
   print(report["montecarlo_std"])

The two should agree to a few percent. Try adding a guide-star altitude
error (``ParameterSpec("gs", "gs_lgs1", "height_shift")``) to ``specs``:
its recovered uncertainty is orders of magnitude larger than the
others', because at these DM altitudes the parallax factor is only
weakly sensitive to a realistic sodium-layer altitude error - a real
system needs a dedicated technique (for example focus sensing) for
that parameter, not geometric registration alone.

Both estimates above are UNCERTAINTY, not BIAS: `analyze` reconstructs
with a single linear step linearized exactly at the truth, so its
expectation equals the truth by construction. To check for a real bias,
run the actual iterative estimator, starting from the nominal system,
over many noisy realizations:

.. code-block:: python

   from synim.registration.analysis import monte_carlo_gauss_newton

   _, mean, covariance = monte_carlo_gauss_newton(
       system, specs, pairs, measurements, sigma, n_trials=300)
   bias = mean - true_alpha
   standard_error = np.sqrt(np.diag(covariance) / 300)
   print(bias / standard_error)  # values of order 1-2 are just sampling noise

See :ref:`registration` for why the two checks are not interchangeable.

Summary
=========

- A ``System`` of ``WFS``/``DM``/``GuideStar`` objects models the
  global mis-registration; ``local_params`` gives the local
  shift/rotation/magnification/anamorphosis a real estimator would see.
- ``gauss_newton_invert`` inverts that map given local measurements.
- Not every set of unknowns is observable: check ``condition_number``
  and ``describe_modes`` before trusting a reconstruction.
- ``analyze`` propagates measurement noise to the global estimate.
