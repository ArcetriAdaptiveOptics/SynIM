"""
Complete example of tomographic reconstructor computation using ParamsManager.

This script demonstrates:
1. Computing interaction matrices on-the-fly (not saved to disk)
2. Computing/loading atmospheric covariance matrices (cached to disk)
3. Computing MMSE reconstructor with noise modeling
4. Saving the final reconstructor

The key advantage: only covariance matrices are saved permanently,
while IMs are computed transiently as needed.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import specula
specula.init(device_idx=0, precision=1)

import synim
synim.init(device_idx=0, precision=1)
from synim.params_manager import ParamsManager

def v_min_max(arr):
    arr = np.abs(arr)
    arr = arr[arr > 0]
    if arr.size == 0:
        vmin = 1e-12
        vmax = 1
    else:
        vmin = max(np.percentile(arr, 1), 1e-12)
        vmax = np.percentile(arr, 99)
        if vmin >= vmax:
            vmax = vmin * 10
    return vmin, vmax

# ===================================================================
# Configuration
# ===================================================================
# Get paths
synim_init_path = os.path.dirname(__file__)
specula_init_path = specula.__file__
specula_package_dir = os.path.dirname(specula_init_path)
specula_repo_path = os.path.dirname(specula_package_dir)

# PRO file and directories
yaml_file = "/raid1/guido/pythonLib/SPECULA_scripts/params_aof/params_aof.yml"
root_dir = "/raid1/guido/PASSATA/AOF/"

print(f"\n{'='*70}")
print(f"TOMOGRAPHIC RECONSTRUCTOR COMPUTATION")
print(f"{'='*70}")
print(f"YAML file: {yaml_file}")
print(f"Root dir: {root_dir}")
print(f"{'='*70}\n")

# Output directories
output_rec_dir = os.path.join(root_dir, "synrec/")
output_im_dir = os.path.join(root_dir, "synim/")
output_cov_dir = os.path.join(root_dir, "covariance/")

# ===================================================================
# Initialize Parameters Manager
# ===================================================================
params_mgr = ParamsManager(yaml_file, root_dir=root_dir, verbose=True)

# Display configuration info
out = params_mgr.count_mcao_stars()
print(f"\nConfiguration summary:")
print(f"  LGS: {out['n_lgs']}")
print(f"  NGS: {out['n_ngs']}")
print(f"  DMs: {out['n_dm']}")
print(f"  Layers: {out['n_rec_layer']}")
print(f"  Optical sources: {out['n_opt']}")
print(f"  Science sources: {out['n_star']}")

# ===================================================================
# Compute Tomographic Reconstructor
# ===================================================================
# This method does everything:
# - Computes IMs on-the-fly (not saved)
# - Computes/loads covariances (saved to disk)
# - Computes MMSE reconstructor
# - Saves reconstructor (optional)

print(f"\n{'='*70}")
print(f"COMPUTING TOMOGRAPHIC RECONSTRUCTOR")
print(f"{'='*70}\n")

# Atmospheric parameters
r0 = 0.2  # Fried parameter [m]
L0 = 25.0  # Outer scale [m]

# Reconstruction parameters
wfs_type = 'lgs'  # Which WFSs to use ('lgs', 'ngs', 'ref')
component_type = 'layer'  # Reconstruct on 'layer' or 'dm'

# Noise variance per WFS [rad^2]
# If None, will be computed from magnitude and detector parameters
# Can also be a list with one value per WFS
noise_variance = None  # Let the method compute it

# Alternative: provide full noise covariance matrix
C_noise = None  # If provided, overrides noise_variance

result = params_mgr.compute_tomographic_reconstructor(
    r0=r0,
    L0=L0,
    wfs_type=wfs_type,
    component_type=component_type,
    noise_variance=noise_variance,
    C_noise=C_noise,
    output_dir=output_rec_dir,
    save=True,  # Save the reconstructor
    verbose=True
)

# ===================================================================
# Save DMs interaction matrix
# ===================================================================

path = params_mgr.save_assembled_interaction_matrix(
    wfs_type=wfs_type,
    component_type='dm',
    output_dir=output_im_dir,
    overwrite=True,
    apply_filter=True,
    verbose=True)

# ===================================================================
# Extract Results
# ===================================================================
reconstructor = result['reconstructor']
im_full = result['im_full']
C_atm_full = result['C_atm_full']
C_noise = result['C_noise']
mode_indices = result['mode_indices']
component_indices = result['component_indices']

print(f"\n{'='*70}")
print(f"RESULTS SUMMARY")
print(f"{'='*70}")
print(f"Reconstructor shape: {reconstructor.shape}")
print(f"  (n_modes, n_slopes) = ({reconstructor.shape[0]}, {reconstructor.shape[1]})")
print(f"\nInteraction matrix shape: {im_full.shape}")
print(f"  (n_modes, n_slopes) = ({im_full.shape[0]}, {im_full.shape[1]})")
print(f"\nAtmospheric covariance shape: {C_atm_full.shape}")
print(f"Noise covariance shape: {C_noise.shape}")
print(f"\nComponents used: {component_indices}")
print(f"Modes per component: {[len(mi) for mi in mode_indices]}")
print(f"Total modes: {sum(len(mi) for mi in mode_indices)}")

if result['rec_filename']:
    print(f"\nReconstructor saved to:")
    print(f"  {result['rec_filename']}")

print(f"{'='*70}\n")

# ===================================================================
# Visualizations
# ===================================================================
print(f"Creating visualizations...")

# 1. Reconstructor matrix (log scale, normalized)
vmin, vmax = v_min_max(reconstructor)
plt.figure(figsize=(12, 8))
plt.imshow(np.abs(reconstructor), cmap='seismic', aspect='auto',
           norm=LogNorm(vmin=vmin,
                        vmax=vmax))
plt.colorbar(label='|Reconstructor| (log scale)')
plt.title(f"Tomographic Reconstructor ({wfs_type.upper()}→{component_type})\n"
          f"r0={r0}m, L0={L0}m")
plt.xlabel("Slope Index")
plt.ylabel("Mode Index")
plt.tight_layout()
plt.savefig(os.path.join(output_rec_dir, "reconstructor_matrix_log.png"), dpi=150)
print(f"  ✓ Saved reconstructor_matrix_log.png")

# 2. Interaction matrix (log scale)
vmin, vmax = v_min_max(im_full)
plt.figure(figsize=(12, 8))
plt.imshow(np.abs(im_full), cmap='viridis', aspect='auto',
           norm=LogNorm(vmin=vmin,
                        vmax=vmax))
plt.colorbar(label='|IM| (log scale)')
plt.title(f"Interaction Matrix ({wfs_type.upper()}→{component_type})")
plt.xlabel("Slope Index")
plt.ylabel("Mode Index")
plt.tight_layout()
plt.savefig(os.path.join(output_rec_dir, "interaction_matrix_log.png"), dpi=150)
print(f"  ✓ Saved interaction_matrix_log.png")

# 3. Atmospheric PRECISION matrix (log scale)
vmin, vmax = v_min_max(C_atm_full)
plt.figure(figsize=(10, 8))
plt.imshow(np.abs(C_atm_full), cmap='viridis',
           norm=LogNorm(vmin=vmin, vmax=vmax))
plt.colorbar(label='|Precision| [rad⁻²] (log scale)')
plt.title(f"Atmospheric Precision Matrix (Inverse Covariance)\n"
          f"r0={r0}m, L0={L0}m")
plt.xlabel("Mode Index")
plt.ylabel("Mode Index")
plt.tight_layout()
plt.savefig(os.path.join(output_rec_dir, "atmospheric_precision_log.png"), dpi=150)
print(f"  ✓ Saved atmospheric_precision_log.png")

# 4. Noise PRECISION matrix (zoomed on first WFS, log scale)
n_slopes_per_wfs = result['n_slopes_per_wfs']
vmin, vmax = v_min_max(C_noise)
plt.figure(figsize=(10, 8))
plt.imshow(np.abs(C_noise[:2*n_slopes_per_wfs, :2*n_slopes_per_wfs]), 
           cmap='viridis',
           norm=LogNorm(vmin=vmin, vmax=vmax))
plt.colorbar(label='|Noise precision| [rad⁻²] (log scale)')
plt.title(f"Noise Precision Matrix (Inverse Covariance) - First WFS")
plt.xlabel("Slope Index")
plt.ylabel("Slope Index")
plt.tight_layout()
plt.savefig(os.path.join(output_rec_dir, "noise_precision_log.png"), dpi=150)
print(f"  ✓ Saved noise_precision_log.png")

# 5. Diagonal of IM @ Reconstructor
plt.figure(figsize=(10, 6))
diag_rec_im = np.diag(reconstructor @ im_full)
plt.plot(diag_rec_im, '.-')
plt.xlabel('Mode Index')
plt.ylabel('Diagonal value')
plt.title('Diagonal of Reconstructor @ IM')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(output_rec_dir, "diag_reconstructor_im.png"), dpi=150)
print(f"  ✓ Saved diag_reconstructor_im.png")

# 6. Diagonal of the inverse atmospheric covariance
plt.figure(figsize=(10, 6))
diag_cov = np.diag(C_atm_full)
plt.plot(diag_cov, '.-')
plt.xlabel('Mode Index')
plt.ylabel('Precision [rad⁻²]')
plt.title('Diagonal of the Inverse Atmospheric Covariance Matrix')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(output_rec_dir, "diag_precision_atmospheric.png"), dpi=150)
print(f"  ✓ Saved diag_precision_atmospheric.png")

# Show all plots
plt.show()

print(f"\n{'='*70}")
print(f"All visualizations saved to: {output_rec_dir}")
print(f"{'='*70}\n")

# ===================================================================
# Optional: Analysis
# ===================================================================
print(f"{'='*70}")
print(f"PERFORMANCE ANALYSIS")
print(f"{'='*70}")

# Condition number analysis
from numpy.linalg import svd, cond

print(f"\nInteraction Matrix:")
print(f"  Condition number: {cond(im_full):.2e}")
u, s, vh = svd(im_full, full_matrices=False)
print(f"  Singular values: min={s.min():.2e}, max={s.max():.2e}")
print(f"  Rank: {np.sum(s > s.max() * 1e-10)}/{len(s)}")

print(f"\nAtmospheric Precision Matrix (Cx^(-1)):")
print(f"  Condition number: {cond(C_atm_full):.2e}")
eigenvalues = np.linalg.eigvalsh(C_atm_full)
print(f"  Eigenvalues: min={eigenvalues.min():.2e}, max={eigenvalues.max():.2e}")

print(f"\nNoise Precision Matrix (Cz^(-1)):")
noise_diag = np.diag(C_noise)
print(f"  Mean precision: {noise_diag.mean():.2e} rad⁻²")
print(f"  Std precision: {noise_diag.std():.2e} rad⁻²")
print(f"  Equivalent noise variance: {1.0/noise_diag.mean():.2e} rad²")

print(f"\nReconstructor Statistics:")
print(f"  Mean: {reconstructor.mean():.2e}")
print(f"  Std: {reconstructor.std():.2e}")
print(f"  Max abs: {np.abs(reconstructor).max():.2e}")
print(f"  Sparsity: {100 * np.sum(reconstructor == 0) / reconstructor.size:.1f}%")

print(f"{'='*70}\n")