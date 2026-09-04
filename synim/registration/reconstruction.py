"""
Global mis-registration reconstruction from LOCAL measurements.

Given a `System` (the Sec. 2 geometry of the SPIE paper, implemented in
`model.py`) and a chosen set of global unknowns (a list of
`ParameterSpec`), this module linearizes the forward map
``alpha -> local parameters`` (`System.local_params` for every WFS-DM
pair) by finite differences - the sensitivity matrix Lambda of Eq. 2 of
the paper, generalized here from a full interaction matrix to the local
geometric parameters (shift, rotation, magnification, anamorphosis)
themselves - and inverts it by Gauss-Newton + pseudo-inverse least
squares, given a vector of LOCAL measurements assumed already available.

This module does NOT estimate local mis-registrations itself and knows
nothing about interaction matrices, probe signals or on-sky calibration
(that is SPRINT's job, already implemented separately in SPECULA - see the
package docstring in `__init__.py`): its only input is a vector of local
(shift, rotation, magnification, anamorphosis) values, one per WFS-DM
pair, with whatever error they carry.
"""
import copy
from dataclasses import dataclass

import numpy as np

__all__ = [
    "LOCAL_DOF", "ParameterSpec", "get_alpha", "apply_alpha",
    "local_params_vector", "jacobian", "gauss_newton_invert",
]

# Which local degrees of freedom to stack per pair, and in which order.
LOCAL_DOF = ("shift_x", "shift_y", "rotation", "magnification",
             "anamorphosis_45", "anamorphosis_90")

_DEFAULT_STEP = {"rotation": 1e-4, "shift": 1e-4, "position_shift": 1e-4,
                  "height_shift": 1e-2, "magnification": 1e-5,
                  "anamorphosis_45": 1e-5, "anamorphosis_90": 1e-5}


@dataclass(frozen=True)
class ParameterSpec:
    """
    One scalar global unknown, i.e. one column of the sensitivity matrix
    Lambda: which object and field of the `System` it controls.

    kind : 'wfs', 'dm' or 'gs'
    name : the WFS/DM/guide-star name in the System (for 'gs', the `.name`
        of the guide star, looked up across all WFSs)
    field : attribute name ('shift', 'rotation', 'magnification',
        'anamorphosis_45', 'anamorphosis_90' for kind='wfs'/'dm';
        'position_shift', 'height_shift' for kind='gs')
    component : vector component (0=x, 1=y) for a tuple field, None for a
        scalar field
    """
    kind: str
    name: str
    field: str
    component: int = None
    label: str = None

    def get_label(self):
        if self.label:
            return self.label
        suffix = "" if self.component is None else ("_x" if self.component == 0 else "_y")
        return f"{self.kind}:{self.name}:{self.field}{suffix}"


def _find_guide_star(system, gs_name):
    for wfs in system.wfss.values():
        if wfs.guide_star.name == gs_name:
            return wfs.guide_star
    raise KeyError(f"No guide star named {gs_name!r} in this system.")


def _get_object(system, kind, name):
    if kind == "wfs":
        return system.wfss[name]
    if kind == "dm":
        return system.dms[name]
    if kind == "gs":
        return _find_guide_star(system, name)
    raise ValueError(f"Unknown ParameterSpec kind: {kind!r}")


def get_alpha(system, specs):
    """Read the current value of every `spec` from `system`."""
    alpha = np.zeros(len(specs))
    for i, spec in enumerate(specs):
        obj = _get_object(system, spec.kind, spec.name)
        value = getattr(obj, spec.field)
        alpha[i] = value[spec.component] if spec.component is not None else value
    return alpha


def apply_alpha(system, specs, alpha):
    """Return a deep copy of `system` with every `spec` set to `alpha[i]`."""
    new_system = copy.deepcopy(system)
    for value, spec in zip(alpha, specs):
        obj = _get_object(new_system, spec.kind, spec.name)
        if spec.component is not None:
            current = list(getattr(obj, spec.field))
            current[spec.component] = float(value)
            setattr(obj, spec.field, tuple(current))
        else:
            setattr(obj, spec.field, float(value))
    return new_system


def local_params_vector(system, pairs, dof=LOCAL_DOF):
    """
    Flatten ``system.local_params(*pair)`` for every `pair` in `pairs`,
    keeping only the requested local degrees of freedom (any subset/order
    of `LOCAL_DOF`). This plays the role of the paper's local IM vector
    D, but expressed directly in geometric parameters.
    """
    values = []
    for wfs_name, dm_name in pairs:
        shift, rotation, magnification, anam45, anam90 = system.local_params(wfs_name, dm_name)
        full = {"shift_x": shift[0], "shift_y": shift[1], "rotation": rotation,
                "magnification": magnification, "anamorphosis_45": anam45,
                "anamorphosis_90": anam90}
        values.extend(full[d] for d in dof)
    return np.array(values, dtype=float)


def jacobian(system, specs, pairs, dof=LOCAL_DOF, step=None):
    """
    Central finite-difference sensitivity matrix ``Lambda = dD/dalpha``
    (Eq. 2 of the SPIE paper, generalized from a full IM to
    `local_params_vector`).

    `step` is either a scalar, a per-spec array, or None (default: a
    small step per field type from `_DEFAULT_STEP`, e.g. 1e-4 deg/shift
    unit and 1e-5 for dimensionless magnification/anamorphosis near 1.0).
    """
    alpha0 = get_alpha(system, specs)
    if step is None:
        step = np.array([_DEFAULT_STEP.get(spec.field, 1e-4) for spec in specs])
    elif np.isscalar(step):
        step = np.full(len(specs), float(step))
    else:
        step = np.asarray(step, dtype=float)

    n_dof = len(pairs) * len(dof)
    Lambda = np.zeros((n_dof, len(specs)))
    for k in range(len(specs)):
        alpha_plus = alpha0.copy(); alpha_plus[k] += step[k]
        alpha_minus = alpha0.copy(); alpha_minus[k] -= step[k]
        D_plus = local_params_vector(apply_alpha(system, specs, alpha_plus), pairs, dof)
        D_minus = local_params_vector(apply_alpha(system, specs, alpha_minus), pairs, dof)
        Lambda[:, k] = (D_plus - D_minus) / (2.0 * step[k])
    return Lambda


def gauss_newton_invert(system0, specs, pairs, D_meas, dof=LOCAL_DOF,
                         n_iter=10, rcond=1e-10, step=None, verbose=False):
    """
    Reconstruct the global parameter vector alpha from a set of (possibly
    noisy) local measurements `D_meas`, by Gauss-Newton iteration starting
    from `system0`'s current values of `specs`:

        alpha_(k+1) = alpha_k + pinv(Lambda_k) @ (D_meas - D(alpha_k))

    This is Eq. 4-7 of the SPIE paper without the IM-amplitude gain matrix
    G of Eq. 6, which has no counterpart once the local measurements are
    already expressed as geometric parameters rather than raw slopes.

    Returns
    -------
    alpha_hat : ndarray
    system_hat : System
        `system0` with `specs` set to `alpha_hat`.
    history : list of float
        ``||D_meas - D(alpha_k)||`` at each iteration, for convergence
        diagnostics.
    """
    alpha = get_alpha(system0, specs)
    history = []
    for _ in range(n_iter):
        system = apply_alpha(system0, specs, alpha)
        D = local_params_vector(system, pairs, dof)
        residual = D_meas - D
        history.append(float(np.linalg.norm(residual)))
        Lambda = jacobian(system, specs, pairs, dof, step=step)
        d_alpha = np.linalg.pinv(Lambda, rcond=rcond) @ residual
        alpha = alpha + d_alpha
        if verbose:
            print(f"iter {len(history) - 1}: |residual|={history[-1]:.3e} "
                  f"|d_alpha|={np.linalg.norm(d_alpha):.3e}")
    system_hat = apply_alpha(system0, specs, alpha)
    return alpha, system_hat, history
