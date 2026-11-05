import os
import numpy as np
import matplotlib.pyplot as plt
import synim
synim.init(device_idx=-1, precision=1)  # CPU, single precision
from synim.params_manager import ParamsManager
import specula
specula.init(device_idx=-1, precision=1)

# -------------------------------------------------------------------
# Configurazione file e directory
yaml_file = "params_aof.yml"
root_dir = "calib_aof/"
output_pm_dir = os.path.join(root_dir, "synpm/")
print(f"YAML file: {yaml_file}")
print(f"Output PM directory: {output_pm_dir}")

# ===================================================================
# Inizializza ParamsManager
# ===================================================================
params_mgr = ParamsManager(yaml_file, root_dir=root_dir, verbose=True)

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

# Non serve più passare reg_factor: viene letto dalla sezione projection del YAML
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

print(f"\nDM projection matrix of the first optical source (pm_full_dm):")
print(f"  Shape: {pm_full_dm.shape}")
print(f"  (1, n_dm_modes, n_pupil_modes)")

print(f"\nLayer projection matrix of the first optical source (pm_full_layer):")
print(f"  Shape: {pm_full_layer.shape}")
print(f"  (1, n_layer_modes, n_pupil_modes)")

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
# compute_tomographic_projection_matrix returns the PMs of the first optical
# source only (the sources are processed one at a time)
if pm_full_dm is not None:
    plt.figure(figsize=(6, 5))
    plt.imshow(pm_full_dm[0, :100, :100], cmap='viridis', aspect='auto')
    plt.colorbar()
    plt.title(f"DM projection matrix, opt1 (w={info['weights'][0]:.2f}), first 100 modes")
    plt.xlabel('Pupil modes')
    plt.ylabel('DM modes')
    plt.tight_layout()
    plt.savefig(os.path.join(output_pm_dir, "pm_dm_opt1.png"), dpi=150)
    print(f"✓ Saved first source visualization")

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

plt.show()

print(f"\n{'='*70}")
print(f"All visualizations saved to: {output_pm_dir}")
print(f"{'='*70}\n")