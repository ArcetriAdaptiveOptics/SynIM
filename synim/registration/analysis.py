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
    "svd_of_jacobian", "condition_number", "describe_mode", "describe_modes",
    "reconstruction_matrix", "covariance_from_noise",
    "monte_carlo_noise_propagation", "analyze",
]


def svd_of_jacobian(Lambda):
    """SVD of the sensitivity matrix Lambda (n_dof x n_specs), singular
    values sorted descending (`numpy.linalg.svd` default)."""
    return np.linalg.svd(Lambda, full_matrices=False)


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
    the `top_n` specs with the largest |coefficient|, as
    ``[(spec.get_label(), coefficient), ...]``, sorted by |coefficient|
    descending.
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
    covariance = np.cov(alpha_hat_samples, rowvar=False)
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
    `montecarlo_covariance`.
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
        })

    return report
