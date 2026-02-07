import os
import numpy as np
import matplotlib.pyplot as plt
from synim.params_manager import ParamsManager
import specula
specula.init(device_idx=-1, precision=1)

# -------------------------------------------------------------------
# Get paths
specula_init_path = specula.__file__
synim_init_path = os.path.dirname(__file__)
specula_package_dir = os.path.dirname(specula_init_path)
specula_repo_path = os.path.dirname(specula_package_dir)

# Set up file paths
#pro_file = os.path.join(synim_init_path, "params_maoryPhaseC_newpup_100m.pro")
pro_file = os.path.join(synim_init_path, "params_maoryPhaseC_newpup_focus_filt.pro")
root_dir = "/raid1/guido/PASSATA/MAORYC/"
print(f"PRO file path: {pro_file}")

output_pm_dir = os.path.join(root_dir, "synpro/")
print(f"Output PM directory: {output_pm_dir}")

# ===================================================================
# Initialize Parameters Manager
# ===================================================================
params_mgr = ParamsManager(pro_file, root_dir=root_dir, verbose=True)

# ===================================================================
# Step 1: Compute/Load Projection Matrices
# ===================================================================
print(f"\n{'='*70}")
print(f"STEP 1: Computing/Loading Projection Matrices")
print(f"{'='*70}\n")

pm_paths = params_mgr.compute_projection_matrices(
    output_dir=output_pm_dir,
    overwrite=False
)

print(f"\n✓ Computed/loaded {len(pm_paths)} projection matrices")

# ===================================================================
# Step 2: Compute Tomographic Projection Matrix
# ===================================================================
print(f"\n{'='*70}")
print(f"STEP 2: Computing Tomographic Projection Matrix")
print(f"{'='*70}\n")

p_opt, pm_full_dm, pm_full_layer, info = params_mgr.compute_tomographic_projection_matrix(
    output_dir=output_pm_dir,
    save=True,
    verbose=True
)

# ===================================================================
# Step 3: Display Results
# ===================================================================
print(f"\n{'='*70}")
print(f"RESULTS SUMMARY")
print(f"{'='*70}")
print(f"Tomographic projection matrix (p_opt):")
print(f"  Shape: {p_opt.shape}")
print(f"  (n_dm_modes, n_layer_modes) = ({p_opt.shape[0]}, {p_opt.shape[1]})")
print(f"  RMS: {np.sqrt(np.mean(p_opt**2)):.4f}")
print(f"  Min/Max: {p_opt.min():.4f} / {p_opt.max():.4f}")

print(f"\nFull DM projection matrix (pm_full_dm):")
print(f"  Shape: {pm_full_dm.shape}")
print(f"  (n_opt_sources, n_dm_modes, n_pupil_modes)")

print(f"\nFull Layer projection matrix (pm_full_layer):")
print(f"  Shape: {pm_full_layer.shape}")
print(f"  (n_opt_sources, n_layer_modes, n_pupil_modes)")

print(f"\nOptical sources:")
print(f"  Count: {info['n_opt_sources']}")
print(f"  Weights: {info['weights']}")

print(f"\nRegularization:")
print(f"  reg_factor: {info['reg_factor']}")
print(f"  rcond: {info['rcond']}")

print(f"\n{'='*70}\n")

# ===================================================================
# Step 4: Visualizations
# ===================================================================
# Visualize p_opt (tomographic projection)
plt.figure(figsize=(12, 8))
plt.imshow(p_opt, cmap='seismic', aspect='auto', vmin=-0.1, vmax=0.1)
plt.colorbar(label='Projection coefficient')
plt.title(f"Tomographic Projection Matrix (p_opt)\n"
          f"DM modes → Layer modes\n")
plt.xlabel("Layer Mode Index")
plt.ylabel("DM Mode Index")
plt.tight_layout()
plt.savefig(os.path.join(output_pm_dir, "p_opt_tomographic.png"), dpi=150)
print(f"✓ Saved p_opt visualization")

# Visualize DM projection for first optical source
if pm_full_dm is not None:
    plt.figure(figsize=(10, 6))
    plt.imshow(pm_full_dm[0, :, :], cmap='viridis', aspect='auto')
    plt.colorbar(label='Projection value')
    plt.title("DM Projection Matrix - First Optical Source (opt1)")
    plt.xlabel("Pupil Mode Index")
    plt.ylabel("DM Mode Index")
    plt.tight_layout()
    plt.savefig(os.path.join(output_pm_dir, "pm_dm_opt1.png"), dpi=150)
    print(f"✓ Saved DM projection visualization")

# Visualize Layer projection for first optical source
if pm_full_layer is not None:
    plt.figure(figsize=(10, 6))
    plt.imshow(pm_full_layer[0, :, :], cmap='viridis', aspect='auto')
    plt.colorbar(label='Projection value')
    plt.title("Layer Projection Matrix - First Optical Source (opt1)")
    plt.xlabel("Pupil Mode Index")
    plt.ylabel("Layer Mode Index")
    plt.tight_layout()
    plt.savefig(os.path.join(output_pm_dir, "pm_layer_opt1.png"), dpi=150)
    print(f"✓ Saved Layer projection visualization")

# Visualize all optical source weights in DM projection
if pm_full_dm is not None and info['n_opt_sources'] > 1:
    plt.figure(figsize=(12, 6))

    # Plot DM mode 0 for all optical sources
    for i in range(min(info['n_opt_sources'], 17)):  # Max 17 sources
        plt.subplot(3, 6, i+1)
        plt.imshow(pm_full_dm[i, :100, :100], cmap='viridis', aspect='auto')
        plt.title(f"opt{i+1}\n(w={info['weights'][i]:.2f})", fontsize=8)
        plt.axis('off')

    plt.suptitle("DM Projection Matrices - All Optical Sources (first 100 modes)",
                 fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(output_pm_dir, "pm_dm_all_sources.png"), dpi=150)
    print(f"✓ Saved all sources visualization")

# Histogram of p_opt values
plt.figure(figsize=(10, 6))
plt.hist(p_opt.flatten(), bins=100, edgecolor='black', alpha=0.7)
plt.xlabel('Projection Coefficient')
plt.ylabel('Frequency')
plt.title(f'Distribution of Tomographic Projection Coefficients\n'
          f'Total elements: {p_opt.size}, Non-zero: {np.sum(p_opt != 0)}')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(output_pm_dir, "p_opt_histogram.png"), dpi=150)
print(f"✓ Saved histogram")

# Show all plots
plt.show()

print(f"\n{'='*70}")
print(f"All visualizations saved to: {output_pm_dir}")
print(f"{'='*70}\n")