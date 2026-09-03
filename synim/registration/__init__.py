"""
SPRINT for WFAO systems - analytic (parameter-space) mis-registration model.

This subpackage implements the geometry and identification approach described
in Agapito, Plantet & Heritier, "SPRINT for WFAO systems", Proc. SPIE 13097,
130975P (2024), generalizing Heritier et al. 2021 (MNRAS 504, 4274) to
multi-WFS/multi-DM (MCAO/GLAO) systems.

Unlike `synim.synim` / `synim.utils`, which synthesize full pixel-based
interaction matrices, this module works purely in mis-registration parameter
space (shift, rotation, magnification, anamorphosis) so that the global <->
local geometry relationship can be composed, inverted and studied
(conditioning, degeneracies, noise propagation) cheaply and repeatedly.
"""
