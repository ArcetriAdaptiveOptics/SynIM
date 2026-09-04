"""
Numerical analysis of the global <-> local sensitivity matrix.

`reconstruction.jacobian` builds Lambda = dD/dalpha, i.e. the exact
analogue of an AO interaction matrix (IM), mapping the global
mis-registration unknowns `alpha` to the local measurements `D` a
SPRINT-like estimator provides. This module treats Lambda accordingly:

  - its SVD gives the "reconstruction matrix" (pseudo-inverse) exactly as
    for a wavefront reconstructor, including truncation of poorly
    determined modes (small singular values);
  - the right singular vectors associated with small singular values are
    the degenerate/poorly-observed combinations of global parameters
    (Sec. 2/5 of the SPIE paper, e.g. GS position error vs WFS-DM shift
    when not all DMs are used) - `describe_mode` labels them in terms of
    the `ParameterSpec`s that dominate them;
  - noise on the local measurements is propagated to the global estimate
    both analytically (through the truncated pseudo-inverse) and by Monte
    Carlo, so the two can be cross-checked against each other.

Nothing here touches interaction matrices, probe signals or on-sky
calibration - see the package docstring in `__init__.py`.
"""
import numpy as np

from .reconstruction import get_alpha, apply_alpha, local_params_vector, LOCAL_DOF

__all__ = [
    "svd_of_jacobian", "normalize_jacobian", "condition_number", "describe_mode",
    "describe_modes", "reconstruction_matrix", "covariance_from_noise",
    "monte_carlo_noise_propagation", "monte_carlo_gauss_newton", "analyze",
]


def svd_of_jacobian(Lambda):
    """SVD of the sensitivity matrix Lambda (n_dof x n_specs), singular
    values sorted descending (`numpy.linalg.svd` default)."""
    return np.linalg.svd(Lambda, full_matrices=False)


def normalize_jacobian(Lambda, specs, scale):
    """
    Rescale each column (unknown) of `Lambda` by a characteristic/typical
    physical size for that parameter, so its SVD's coefficients become
    comparable across parameter TYPES with very different natural units -
    a shift in `pixel_pitch` units, a rotation in degrees, a
    magnification/anamorphosis expressed as a near-1 dimensionless ratio,
    and so on. Without this, a fixed coefficient threshold (as used by
    `describe_mode`) is unfair: it takes a much larger PHYSICAL change in
    magnification/anamorphosis than in shift to move their raw
    coefficient by the same amount, so a uniform threshold silently
    over-weights shift-like parameters and under-weights the others.

    Rescaling columns by positive factors never changes WHICH directions
    are exactly degenerate (a rank-deficient matrix stays rank-deficient
    under any nonzero column scaling - only the ORTHOGONAL BASIS chosen
    within a many-dimensional null space, and hence how fairly its
    coefficients can be compared, depends on this scaling). Compute the
    SVD of the returned, rescaled matrix (not of the original `Lambda`)
    for `describe_mode`'s labels to be meaningful across parameter types.

    Parameters
    ----------
    Lambda : ndarray, shape (n_dof, n_specs)
    specs : list of reconstruction.ParameterSpec
    scale : dict or array-like
        Either an array of length `len(specs)` (one characteristic size
        per spec, in that spec's own unit), or a dict mapping
        `spec.field` (e.g. `"shift"`, `"rotation"`, `"magnification"`,
        `"position_shift"`, ...) to a characteristic size.

    Returns
    -------
    Lambda_scaled : ndarray, same shape as `Lambda`
    """
    if isinstance(scale, dict):
        try:
            scale = np.array([scale[spec.field] for spec in specs], dtype=float)
        except KeyError as exc:
            raise KeyError(f"no characteristic scale given for field {exc}") from exc
    else:
        scale = np.asarray(scale, dtype=float)
    return Lambda * scale[np.newaxis, :]


def condition_number(Lambda):
    """``sigma_max / sigma_min`` of Lambda (``inf`` if Lambda is rank
    deficient, i.e. a global combination of parameters has zero effect on
    every local measurement - a hard/exact degeneracy)."""
    s = np.linalg.svd(Lambda, compute_uv=False)
    if s[-1] == 0:
        return np.inf
    return float(s[0] / s[-1])


def describe_mode(v, specs, top_n=3):
    """
    Human-readable description of one right singular vector `v` (a linear
    combination of `specs`, i.e. one column of Lambda's V / row of Vt):
    the `top_n` specs with the largest absolute coefficient, as
    ``[(spec.get_label(), coefficient), ...]``, sorted by absolute
    coefficient descending.
    """
    order = np.argsort(-np.abs(v))[:top_n]
    return [(specs[i].get_label(), float(v[i])) for i in order]


def describe_modes(Vt, specs, top_n=3):
    """`describe_mode` for every row of `Vt` (i.e. every singular
    vector), in the same order as the corresponding singular values."""
    return [describe_mode(row, specs, top_n=top_n) for row in Vt]


def reconstruction_matrix(Lambda, n_modes=None, rcond=None):
    """
    Build the (pseudo-)inverse "reconstruction matrix" R such that
    ``alpha_hat = R @ D``, from the SVD of `Lambda`, with explicit modal
    truncation:

      - `n_modes`: keep only the `n_modes` best-determined modes (largest
        singular values), discarding the rest - the direct analogue of
        modal filtering in AO reconstruction;
      - `rcond`: instead keep only singular values >= `rcond * sigma_max`
        (same convention as `numpy.linalg.pinv`);
      - neither given: keep all nonzero singular values (`rcond=1e-10`).

    Returns
    -------
    R : ndarray, shape (n_specs, n_dof)
    n_kept : int
        Number of modes actually kept.
    singular_values : ndarray
        All singular values of `Lambda`, for reference/plotting.
    """
    U, S, Vt = svd_of_jacobian(Lambda)
    if n_modes is not None:
        n_kept = min(n_modes, len(S))
    else:
        threshold = (rcond if rcond is not None else 1e-10) * S[0]
        n_kept = int(np.sum(S >= threshold))

    S_inv = np.zeros_like(S)
    S_inv[:n_kept] = 1.0 / S[:n_kept]
    R = (Vt.T * S_inv) @ U.T
    return R, n_kept, S


def covariance_from_noise(Lambda, sigma, n_modes=None, rcond=None):
    """
    Analytic propagation of local-measurement noise to the global
    parameter estimate, through the (possibly truncated)
    `reconstruction_matrix`:

        Cov(alpha_hat) = R @ diag(sigma**2) @ R.T

    `sigma` is either a scalar (same std for every local measurement) or
    an array of length `Lambda.shape[0]` (one std per local d.o.f.).

    Returns
    -------
    covariance : ndarray, shape (n_specs, n_specs)
    std : ndarray, shape (n_specs,)
        ``sqrt(diag(covariance))``, directly comparable to the paper's
        Table 1 style accuracy numbers.
    """
    R, n_kept, S = reconstruction_matrix(Lambda, n_modes=n_modes, rcond=rcond)
    sigma = np.broadcast_to(np.asarray(sigma, dtype=float), (Lambda.shape[0],))
    covariance = (R * sigma ** 2) @ R.T
    return covariance, np.sqrt(np.diag(covariance))


def monte_carlo_noise_propagation(system0, specs, pairs, D_true, sigma,
                                   dof=LOCAL_DOF, n_trials=2000, n_modes=None,
                                   rcond=None, rng=None):
    """
    Monte Carlo cross-check of `covariance_from_noise`: repeatedly draw a
    noisy local measurement ``D_meas = D_true + noise`` and reconstruct
    `alpha` with a single linear step around `system0`'s current
    `specs` values (i.e. using a fixed Lambda, exactly matching the
    linear model behind the analytic covariance formula - not a full
    `gauss_newton_invert`, which would additionally re-linearize at every
    iteration).

    Because this reconstruction is linear and linearized exactly at
    `system0`'s own values, its expectation over noise equals those
    values by construction: it can validate the analytic COVARIANCE, but
    it is unbiased by construction and therefore cannot reveal a real
    estimation BIAS. For that, use `monte_carlo_gauss_newton`, which runs
    the actual iterative estimator instead - typically starting from a
    different (e.g. nominal) system.

    Returns
    -------
    alpha_hat_samples : ndarray, shape (n_trials, n_specs)
    mean : ndarray, shape (n_specs,)
    covariance : ndarray, shape (n_specs, n_specs)
    """
    from .reconstruction import jacobian  # local import: avoids a cycle at module load

    if rng is None:
        rng = np.random.default_rng()
    sigma = np.broadcast_to(np.asarray(sigma, dtype=float), (len(pairs) * len(dof),))

    alpha0 = get_alpha(system0, specs)
    Lambda = jacobian(system0, specs, pairs, dof=dof)
    R, n_kept, S = reconstruction_matrix(Lambda, n_modes=n_modes, rcond=rcond)
    D0 = local_params_vector(system0, pairs, dof=dof)

    noise = rng.normal(scale=sigma, size=(n_trials, len(sigma)))
    D_meas = D_true + noise
    alpha_hat_samples = alpha0 + (D_meas - D0) @ R.T

    mean = alpha_hat_samples.mean(axis=0)
    # `numpy.cov` returns a bare 0-d scalar (not a 1x1 array) for a single
    # variable, which breaks callers expecting a (n_specs, n_specs) matrix
    # (e.g. `numpy.diag`) - reshape defensively.
    covariance = np.cov(alpha_hat_samples, rowvar=False).reshape(len(specs), len(specs))
    return alpha_hat_samples, mean, covariance


def monte_carlo_gauss_newton(system_start, specs, pairs, D_true, sigma,
                              dof=LOCAL_DOF, n_trials=500, n_iter=6,
                              rcond=1e-10, rng=None):
    """
    Monte Carlo study of the REAL iterative estimator
    (`reconstruction.gauss_newton_invert`): for each of `n_trials` noisy
    local measurements ``D_meas = D_true + noise``, run the full
    Gauss-Newton iteration from `system_start` (typically the NOMINAL
    system, not the true one - i.e. what an operational estimator
    actually starts from) to convergence, re-linearizing at every step.

    Unlike `monte_carlo_noise_propagation` (a single linear step
    linearized exactly at the truth, unbiased by construction - see its
    docstring), this can reveal a genuine estimation BIAS, from:

      - `system_start` not being the true system (a real estimator does
        not know the truth in advance);
      - the nonlinearity of the local <-> global map itself (the same
        nonlinearity that makes `gauss_newton_invert` take several
        iterations even without noise);
      - noise interacting with a poorly conditioned direction (large
        noise-driven excursions can leave the region where the
        linearization is a good approximation).

    This is `n_iter` times more expensive per trial than
    `monte_carlo_noise_propagation` (one `jacobian` evaluation per
    iteration instead of one total), so `n_trials` defaults lower.

    Returns
    -------
    alpha_hat_samples : ndarray, shape (n_trials, n_specs)
    mean : ndarray, shape (n_specs,)
        Compare to the true `alpha` (e.g. via `apply_alpha`/`get_alpha`
        on the system `D_true` was generated from): a nonzero
        ``mean - true_alpha`` is the estimator's bias.
    covariance : ndarray, shape (n_specs, n_specs)
    """
    from .reconstruction import gauss_newton_invert  # local import: avoids a cycle at module load

    if rng is None:
        rng = np.random.default_rng()
    sigma = np.broadcast_to(np.asarray(sigma, dtype=float), (len(pairs) * len(dof),))

    alpha_hat_samples = np.empty((n_trials, len(specs)))
    for trial in range(n_trials):
        noise = rng.normal(scale=sigma)
        D_meas = D_true + noise
        alpha_hat, _, _ = gauss_newton_invert(
            system_start, specs, pairs, D_meas, dof=dof, n_iter=n_iter, rcond=rcond)
        alpha_hat_samples[trial] = alpha_hat

    mean = alpha_hat_samples.mean(axis=0)
    # `numpy.cov` returns a bare 0-d scalar (not a 1x1 array) for a single
    # variable, which breaks callers expecting a (n_specs, n_specs) matrix
    # (e.g. `numpy.diag`) - reshape defensively.
    covariance = np.cov(alpha_hat_samples, rowvar=False).reshape(len(specs), len(specs))
    return alpha_hat_samples, mean, covariance


def analyze(system0, specs, pairs, dof=LOCAL_DOF, sigma=None, n_modes=None,
            rcond=None, n_trials=2000, rng=None, top_n=3):
    """
    Convenience report tying together sensitivity, degeneracy and (if
    `sigma` is given) noise propagation for one `(system0, specs, pairs)`
    setup, evaluated at `system0`'s current values of `specs`.

    Returns a dict with: `Lambda`, `singular_values`, `condition_number`,
    `modes` (`describe_modes` output, ordered like `singular_values`,
    i.e. best-determined first), and - only if `sigma` is given -
    `analytic_std`, `analytic_covariance`, `montecarlo_std`,
    `montecarlo_covariance`, `montecarlo_mean`.

    The Monte Carlo entries come from `monte_carlo_noise_propagation`
    (a single linear step linearized at `system0`'s own values): good to
    cross-check `analytic_std`/`analytic_covariance`, but unbiased by
    construction, so `montecarlo_mean` is not a real bias estimate - use
    `monte_carlo_gauss_newton` directly for that.
    """
    from .reconstruction import jacobian  # local import: avoids a cycle at module load

    Lambda = jacobian(system0, specs, pairs, dof=dof)
    U, S, Vt = svd_of_jacobian(Lambda)

    report = {
        "Lambda": Lambda,
        "singular_values": S,
        "condition_number": condition_number(Lambda),
        "modes": describe_modes(Vt, specs, top_n=top_n),
    }

    if sigma is not None:
        analytic_cov, analytic_std = covariance_from_noise(
            Lambda, sigma, n_modes=n_modes, rcond=rcond)
        D_true = local_params_vector(system0, pairs, dof=dof)
        _, mc_mean, mc_cov = monte_carlo_noise_propagation(
            system0, specs, pairs, D_true, sigma, dof=dof, n_trials=n_trials,
            n_modes=n_modes, rcond=rcond, rng=rng)
        report.update({
            "analytic_covariance": analytic_cov,
            "analytic_std": analytic_std,
            "montecarlo_covariance": mc_cov,
            "montecarlo_std": np.sqrt(np.diag(mc_cov)),
            # mc_mean == get_alpha(system0, specs) up to Monte Carlo noise,
            # by construction (see monte_carlo_noise_propagation) - exposed
            # mainly as a sanity check, not as a bias estimate; for a real
            # bias estimate see `monte_carlo_gauss_newton`.
            "montecarlo_mean": mc_mean,
        })

    return report
