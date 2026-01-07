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

os.environ['SPECULA_DISABLE_GPU'] = 'TRUE' # Force CPU usage for Specula
import specula
specula.init(device_idx=-1, precision=1)

import synim
synim.init(device_idx=0, precision=1)
from synim.params_manager import ParamsManager

# ===================================================================
# Configuration
# ===================================================================
# Get paths
synim_init_path = os.path.dirname(__file__)
specula_init_path = specula.__file__
specula_package_dir = os.path.dirname(specula_init_path)
specula_repo_path = os.path.dirname(specula_package_dir)

# PRO file and directories
pro_file = os.path.join(synim_init_path, "params_maoryPhaseC_newpup_focus_filt.pro")
root_dir = "/raid1/guido/PASSATA/MAORYC/"

print(f"\n{'='*70}")
print(f"TOMOGRAPHIC RECONSTRUCTOR COMPUTATION")
print(f"{'='*70}")
print(f"PRO file: {pro_file}")
print(f"Root dir: {root_dir}")
print(f"{'='*70}\n")

# Output directories
output_rec_dir = os.path.join(root_dir, "synrec/")
output_cov_dir = os.path.join(root_dir, "covariance/")

# ===================================================================
# Initialize Parameters Manager
# ===================================================================
params_mgr = ParamsManager(pro_file, root_dir=root_dir, verbose=True)

# Display configuration info
out = params_mgr.count_mcao_stars()
print(f"\nConfiguration summary:")
print(f"  LGS: {out['n_lgs']}")
print(f"  NGS: {out['n_ngs']}")
print(f"  REF: {out['n_ref']}")
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
    skip_gpu_covariance=True,  # Use CPU for covariance computation
    verbose=True
)

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

# 1. Reconstructor matrix
plt.figure(figsize=(12, 8))
plt.imshow(reconstructor, cmap='seismic', aspect='auto', 
           vmin=-np.percentile(np.abs(reconstructor), 99),
           vmax=np.percentile(np.abs(reconstructor), 99))
plt.colorbar(label='Reconstructor coefficient')
plt.title(f"Tomographic Reconstructor ({wfs_type.upper()}→{component_type})\n"
          f"r0={r0}m, L0={L0}m, reg={reg_factor:.0e}")
plt.xlabel("Slope Index")
plt.ylabel("Mode Index")
plt.tight_layout()
plt.savefig(os.path.join(output_rec_dir, "reconstructor_matrix.png"), dpi=150)
print(f"  ✓ Saved reconstructor_matrix.png")

# 2. Interaction matrix
plt.figure(figsize=(12, 8))
plt.imshow(im_full, cmap='viridis', aspect='auto')
plt.colorbar(label='IM coefficient')
plt.title(f"Interaction Matrix ({wfs_type.upper()}→{component_type})")
plt.xlabel("Slope Index")
plt.ylabel("Mode Index")
plt.tight_layout()
plt.savefig(os.path.join(output_rec_dir, "interaction_matrix.png"), dpi=150)
print(f"  ✓ Saved interaction_matrix.png")

# 3. Atmospheric covariance
plt.figure(figsize=(10, 8))
plt.imshow(C_atm_full, cmap='viridis')
plt.colorbar(label='Covariance [rad²]')
plt.title(f"Atmospheric Covariance Matrix\n"
          f"r0={r0}m, L0={L0}m, weights={weights}")
plt.xlabel("Mode Index")
plt.ylabel("Mode Index")
plt.tight_layout()
plt.savefig(os.path.join(output_rec_dir, "covariance_atmospheric.png"), dpi=150)
print(f"  ✓ Saved covariance_atmospheric.png")

# 4. Noise covariance (zoomed on first WFS)
n_slopes_per_wfs = result['n_slopes_per_wfs']
plt.figure(figsize=(10, 8))
plt.imshow(C_noise[:2*n_slopes_per_wfs, :2*n_slopes_per_wfs], cmap='viridis')
plt.colorbar(label='Noise variance [rad²]')
plt.title(f"Noise Covariance Matrix (first WFS)")
plt.xlabel("Slope Index")
plt.ylabel("Slope Index")
plt.tight_layout()
plt.savefig(os.path.join(output_rec_dir, "covariance_noise.png"), dpi=150)
print(f"  ✓ Saved covariance_noise.png")

# 5. Reconstructor histogram
plt.figure(figsize=(10, 6))
plt.hist(reconstructor.flatten(), bins=100, edgecolor='black', alpha=0.7)
plt.xlabel('Reconstructor Coefficient')
plt.ylabel('Frequency')
plt.title(f'Distribution of Reconstructor Coefficients\n'
          f'Total: {reconstructor.size}, Non-zero: {np.sum(reconstructor != 0)}')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(output_rec_dir, "reconstructor_histogram.png"), dpi=150)
print(f"  ✓ Saved reconstructor_histogram.png")

# 6. Mode distribution per component
plt.figure(figsize=(12, 6))
for i, (comp_idx, modes) in enumerate(zip(component_indices, mode_indices)):
    plt.subplot(1, len(component_indices), i+1)
    
    # Extract reconstructor for this component
    start_mode = sum(len(mi) for mi in mode_indices[:i])
    end_mode = start_mode + len(modes)
    rec_component = reconstructor[start_mode:end_mode, :]
    
    # Plot RMS per mode
    rms_per_mode = np.sqrt(np.mean(rec_component**2, axis=1))
    plt.plot(modes, rms_per_mode, 'o-', markersize=3)
    plt.xlabel('Mode Index')
    plt.ylabel('RMS Coefficient')
    plt.title(f'{component_type.capitalize()}{comp_idx}\n'
              f'({len(modes)} modes)')
    plt.grid(True, alpha=0.3)

plt.suptitle(f'Reconstructor RMS per Mode and Component', fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(output_rec_dir, "reconstructor_per_component.png"), dpi=150)
print(f"  ✓ Saved reconstructor_per_component.png")

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

print(f"\nAtmospheric Covariance:")
print(f"  Condition number: {cond(C_atm_full):.2e}")
eigenvalues = np.linalg.eigvalsh(C_atm_full)
print(f"  Eigenvalues: min={eigenvalues.min():.2e}, max={eigenvalues.max():.2e}")

print(f"\nNoise Covariance:")
noise_diag = np.diag(C_noise)
print(f"  Mean variance: {noise_diag.mean():.2e} rad²")
print(f"  Std variance: {noise_diag.std():.2e} rad²")

print(f"\nReconstructor Statistics:")
print(f"  Mean: {reconstructor.mean():.2e}")
print(f"  Std: {reconstructor.std():.2e}")
print(f"  Max abs: {np.abs(reconstructor).max():.2e}")
print(f"  Sparsity: {100 * np.sum(reconstructor == 0) / reconstructor.size:.1f}%")

print(f"{'='*70}\n")