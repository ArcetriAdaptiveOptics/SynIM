import numpy as np
import matplotlib.pyplot as plt
import synim
synim.init(device_idx=-1, precision=1)  # CPU, single precision
import synim.synim as synim

import specula
specula.init(-1)  # CPU
from specula.data_objects.ifunc import IFunc
from specula.lib.make_mask import make_mask

# pupil parameters
pup_m = 38.5
pup_npoints = 480
pup_obsratio = 0.283
pixel_pitch = pup_m / pup_npoints

# DM parameters
dm_npoints = 600
dm_obsratio = 0.0
dm_height = 10000
n_modes = 100
dm_rotation = 22.5

# WFS parameters
# number of subapertures of the WFS (sub-multiple of pup_npoints)
nsubaps = 20
wfs_rotation = 0
wfs_translation = (0,0)
wfs_magnification = (1,1)
wfs_fov_arcsec = 2.4

# LGS parameters
LGS_pol_coo = [45,30] #arcsec, deg
LGS_height = 90000*1/np.cos(30/180*np.pi)

print('DM height [m]', dm_height)
print('LGS height [m]', LGS_height)
print('LGS star pos. [arcsec, deg]', LGS_pol_coo)

# Generate Zernike influence functions
ifunc = IFunc(
    type_str='zernike',
    nmodes=n_modes,
    npixels=dm_npoints,
    obsratio=0.0,
    diaratio=1.0
)
dm_mask = np.asarray(ifunc.mask_inf_func, dtype=np.float32)
dm_array = ifunc.ifunc_2d_to_3d(normalize=False)

# Pupil masks
pup_mask = make_mask(
    pup_npoints, obsratio=pup_obsratio, diaratio=1.0, xc=0.0, yc=0.0, square=False, inverse=False, centeronpixel=False)

# estimate the valid sub-apertures indices
pup_mask_sa = synim.rebin(pup_mask, (nsubaps,nsubaps), method='sum')
idx_valid_sa = np.ravel(np.array(np.where(pup_mask_sa.flat > (0.5*np.max(pup_mask_sa)))).astype(np.int32))

print('Masks and Zernike computed.')

fig, axs = plt.subplots(1,2)
im1 = axs[0].imshow(pup_mask)
im1 = axs[1].imshow(dm_mask)
axs[0].set_title('Pupil masks')
axs[1].set_title('DM masks')
fig.colorbar(im1, ax=axs.ravel().tolist(),fraction=0.02)

fig, axs = plt.subplots(1,2)
im2 = axs[0].imshow(dm_array[:,:,0])
im2 = axs[1].imshow(dm_array[:,:,4])
axs[0].set_title('DM shapes (idx 0 and 4)')
fig.colorbar(im2, ax=axs.ravel().tolist(),fraction=0.02)
plt.show()

intmat = synim.interaction_matrix(pup_m,pup_mask,dm_array,
                                  dm_mask,dm_height,dm_rotation,
                                  nsubaps,wfs_rotation,wfs_translation,
                                  wfs_magnification,wfs_fov_arcsec,LGS_pol_coo,LGS_height,
                                  idx_valid_sa=idx_valid_sa,verbose=True,display=True,
                                  specula_convention=False)

fig, _ = plt.subplots()
plt.imshow(intmat)
plt.title('Interaction Matrix')
plt.colorbar()
plt.show()

plt.ioff()
plt.show()  # Show the plots