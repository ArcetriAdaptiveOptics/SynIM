

{main,
  root_dir: '/home/guido/IDLWorkspace/PASSATA/MAORYC/',       ; Root directory for calibration manager
  instrument_name: 'maory'
  verbose: 0b,
  pixel_pupil: 480,                                             ; Linear dimension of pupil phase array
  pixel_pitch: 38.5/480,                                         ; [m] Pitch of the pupil phase array (size of 1 element)
  total_time: 2.000,                                            ; [s] Total simulation running time
  time_step : 0.002,                                            ; [s] Simulation time step
  zenithAngleInDeg : 30
  precision: 0,                                                 ; Calculation precision, 0=float 1=double
}

{dm1,
  ifunc_tag : 'M4_eso_zonal_486p_5352a_slaved',                 ; DM command-to-phase influence matrix
  m2c_tag   : 'M4_eso_39.5m_486p_20230215'
  nmodes    : 4519L                                             ; number of modes
  height    : 600.                                              ; DM height [m]
}

{dm2,
  ifunc_tag : 'MAORY_542pix_37nacts_cir_0.000obs_zonal_ifs'                   
  m2c_tag   : 'MAORY_542pix_37nacts_cir_0.000obs_5zern_10000.0cn'
  nmodes    : 1026L
  ;ifunc_tag : 'MAORY_542pix_31nacts_cir_0.000obs_zonal_ifs'
  ;m2c_tag   : 'MAORY_542pix_31nacts_cir_0.000obs_5zern_10000.0cn'
  ;nmodes    : 720L
  start_mode: 5                                                 
  height    : 6500.  
}

{dm3,
  ifunc_tag : 'MAORY_648pix_37nacts_cir_0.000obs_zonal_ifs'                   
  m2c_tag   : 'MAORY_648pix_37nacts_cir_0.000obs_5zern_10000.0cn'
  nmodes    : 1026L
  start_mode: 2                                                 
  height    : 17500.  
}

{launcherPos1,
    position  : [ 15.2, 15.2, 0]
    sizeArcsec: 1.5                                              
}

{launcherPos2,
    position  : [ 15.2, 15.2, 0]
    sizeArcsec: 1.5                                              
}

{launcherPos3,
    position  : [-15.2, 15.2, 0]
    sizeArcsec: 1.5
}

{launcherPos4,
    position  : [-15.2,-15.2, 0]
    sizeArcsec: 1.5
}

{launcherPos5,
    position  : [-15.2,-15.2, 0]
    sizeArcsec: 1.5
}

{launcherPos6,
    position  : [ 15.2,-15.2, 0]
    sizeArcsec: 1.5
}

{sh_lgs1,
   wavelengthInNm   : 589
   sensor_fov       : 16.1,
   sensor_pxscale   : 16.1/14,
   sensor_npx       : 14
   subap_on_diameter: 68
   energy_th        : 0.5                                       ; energy threshold for including subapertures 
   subapdata_tag    : 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_rot6.2'
   convolGaussSpotSize: 1.8
   kernel_application : 'FFT'
   fov_ovs_coeff: 1.52
   rotAnglePhInDeg: 6.2 ; this is associated with launcher @ [ 13.5,  16.8]
   xShiftPhInPixel: 0.
   yShiftPhInPixel: 0.
}

{sh_lgs2,
   wavelengthInNm   : 589
   sensor_fov       : 16.1,
   sensor_pxscale   : 16.1/14,
   sensor_npx       : 14
   subap_on_diameter: 68
   energy_th        : 0.5                                       ; energy threshold for including subapertures 
   subapdata_tag    : 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_rot14.2'
   convolGaussSpotSize: 1.8
   kernel_application : 'FFT'
   fov_ovs_coeff: 1.52
   rotAnglePhInDeg: 14.2 ; this is associated with launcher @ [ 10.9,  18.3]
   xShiftPhInPixel: 0.
   yShiftPhInPixel: 0.
}

{sh_lgs3,
   wavelengthInNm   : 589
   sensor_fov       : 16.1,
   sensor_pxscale   : 16.1/14,
   sensor_npx       : 14
   subap_on_diameter: 68
   energy_th        : 0.5                                       ; energy threshold for including subapertures 
   subapdata_tag    : 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_rot-6.2'
   convolGaussSpotSize: 1.8
   kernel_application : 'FFT'
   fov_ovs_coeff: 1.52
   rotAnglePhInDeg: -6.2 ; this is associated with launcher @ [-13.5,  16.8]
   xShiftPhInPixel: 0.
   yShiftPhInPixel: 0.
}

{sh_lgs4,
   wavelengthInNm   : 589
   sensor_fov       : 16.1,
   sensor_pxscale   : 16.1/14,
   sensor_npx       : 14
   subap_on_diameter: 68
   energy_th        : 0.5                                       ; energy threshold for including subapertures 
   subapdata_tag    : 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_rot6.2'
   convolGaussSpotSize: 1.8
   kernel_application : 'FFT'
   fov_ovs_coeff: 1.52
   rotAnglePhInDeg: 6.2 ; this is associated with launcher @ [-13.5, -16.8]
   xShiftPhInPixel: 0.
   yShiftPhInPixel: 0.
}

{sh_lgs5,
   wavelengthInNm   : 589
   sensor_fov       : 16.1,
   sensor_pxscale   : 16.1/14,
   sensor_npx       : 14
   subap_on_diameter: 68
   energy_th        : 0.5                                       ; energy threshold for including subapertures 
   subapdata_tag    : 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_rot14.2'
   convolGaussSpotSize: 1.8
   kernel_application : 'FFT'
   fov_ovs_coeff: 1.52
   rotAnglePhInDeg: 14.2 ; this is associated with launcher @ [-10.9, -18.3]
   xShiftPhInPixel: 0.
   yShiftPhInPixel: 0.

}

{sh_lgs6,
   wavelengthInNm   : 589
   sensor_fov       : 16.1,
   sensor_pxscale   : 16.1/14,
   sensor_npx       : 14
   subap_on_diameter: 68
   energy_th        : 0.5                                       ; energy threshold for including subapertures 
   subapdata_tag    : 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_rot-6.2'
   convolGaussSpotSize: 1.8
   kernel_application : 'FFT'
   fov_ovs_coeff: 1.52
   rotAnglePhInDeg: -6.2 ; this is associated with launcher @ [ 13.5, -16.8]
   xShiftPhInPixel: 0.
   yShiftPhInPixel: 0.
}

{sh_ngs1,
   wavelengthInNm   : 1650,
   sensor_fov       : 1.8,
   sensor_pxscale   : 1.8/120,
   sensor_npx       : 120,
   subap_on_diameter: 2
   energy_th        : 0.5                                       ; energy threshold for including subapertures 
   subapdata_tag    : 'auto'                                    ; (compute by give_me_the_tag_lngs)
}

{sh_ngs2,
   wavelengthInNm   : 1650,
   sensor_fov       : 1.8,
   sensor_pxscale   : 1.8/120,
   sensor_npx       : 120,
   subap_on_diameter: 2
   energy_th        : 0.5                                       ; energy threshold for including subapertures 
   subapdata_tag    : 'auto'                                    ; (compute by give_me_the_tag_lngs)
}

{sh_ngs3,
   wavelengthInNm   : 1650,
   sensor_fov       : 1.8,
   sensor_pxscale   : 1.8/120,
   sensor_npx       : 120,
   subap_on_diameter: 2
   energy_th        : 0.5                                       ; energy threshold for including subapertures 
   subapdata_tag    : 'auto'                                    ; (compute by give_me_the_tag_lngs)
}

{detector1,
  name: 'CCD220',                                               ; name of the detector (used by ccd auto_params_management method)
  size: [952,952]      ; Detector size in pixels
  dt  : 0.002d         ; Detector integration time
  bandw: 20            ; [nm] Sensor bandwidth
  binning: 1           ; Detector binning (1,2,3)

  photon_noise      : 1b    ; activate photon noise
  excess_noise      : 0b,                                             ; activate Excess noise of sqrt(2.) characteristic of EMCCDs.
  background_noise  : 0b,   ; activate sky background noise
  darkcurrent_noise : 0b,   ; activate dark current noise
  readout_noise     : 1b    ; activate detector readout noise
  readout_level     : 3.0   ; readout noise in e-/pix/frame
  darkcurrent_level : 0
  background_level  : 0
  sky_bg_norm       : 0
  quantum_eff       : 0.196
}

{detector2,
  name: 'CCD220',                                               ; name of the detector (used by ccd auto_params_management method)
  size: [952,952]      ; Detector size in pixels
  dt  : 0.002d         ; Detector integration time
  bandw: 20            ; [nm] Sensor bandwidth
  binning: 1           ; Detector binning (1,2,3)

  photon_noise      : 1b    ; activate photon noise
  excess_noise      : 0b,                                             ; activate Excess noise of sqrt(2.) characteristic of EMCCDs.
  background_noise  : 0b,   ; activate sky background noise
  darkcurrent_noise : 0b,   ; activate dark current noise
  readout_noise     : 1b    ; activate detector readout noise
  readout_level     : 3.0   ; readout noise in e-/pix/frame
  darkcurrent_level : 0
  background_level  : 0
  sky_bg_norm       : 0
  quantum_eff       : 0.196  
}

{detector3,
  name: 'CCD220',                                               ; name of the detector (used by ccd auto_params_management method)
  size: [952,952]      ; Detector size in pixels
  dt  : 0.002d         ; Detector integration time
  bandw: 20            ; [nm] Sensor bandwidth
  binning: 1           ; Detector binning (1,2,3)

  photon_noise      : 1b    ; activate photon noise
  excess_noise      : 0b,                                             ; activate Excess noise of sqrt(2.) characteristic of EMCCDs.
  background_noise  : 0b,   ; activate sky background noise
  darkcurrent_noise : 0b,   ; activate dark current noise
  readout_noise     : 1b    ; activate detector readout noise
  readout_level     : 3.0   ; readout noise in e-/pix/frame
  darkcurrent_level : 0
  background_level  : 0
  sky_bg_norm       : 0
  quantum_eff       : 0.196   
}

{detector4,
  name: 'CCD220',                                               ; name of the detector (used by ccd auto_params_management method)
  size: [952,952]      ; Detector size in pixels
  dt  : 0.002d         ; Detector integration time
  bandw: 20            ; [nm] Sensor bandwidth
  binning: 1           ; Detector binning (1,2,3)

  photon_noise      : 1b    ; activate photon noise
  excess_noise      : 0b,                                             ; activate Excess noise of sqrt(2.) characteristic of EMCCDs.
  background_noise  : 0b,   ; activate sky background noise
  darkcurrent_noise : 0b,   ; activate dark current noise
  readout_noise     : 1b    ; activate detector readout noise
  readout_level     : 3.0   ; readout noise in e-/pix/frame
  darkcurrent_level : 0
  background_level  : 0
  sky_bg_norm       : 0
  quantum_eff       : 0.196  
}

{detector5,
  name: 'CCD220',                                               ; name of the detector (used by ccd auto_params_management method)
  size: [952,952]      ; Detector size in pixels
  dt  : 0.002d         ; Detector integration time
  bandw: 20            ; [nm] Sensor bandwidth
  binning: 1           ; Detector binning (1,2,3)

  photon_noise      : 1b    ; activate photon noise
  excess_noise      : 0b,                                             ; activate Excess noise of sqrt(2.) characteristic of EMCCDs.
  background_noise  : 0b,   ; activate sky background noise
  darkcurrent_noise : 0b,   ; activate dark current noise
  readout_noise     : 1b    ; activate detector readout noise
  readout_level     : 3.0   ; readout noise in e-/pix/frame
  darkcurrent_level : 0
  background_level  : 0
  sky_bg_norm       : 0
  quantum_eff       : 0.196     
}

{detector6,
  name: 'CCD220',                                               ; name of the detector (used by ccd auto_params_management method)
  size: [952,952]      ; Detector size in pixels
  dt  : 0.002d         ; Detector integration time
  bandw: 20            ; [nm] Sensor bandwidth
  binning: 1           ; Detector binning (1,2,3)

  photon_noise      : 1b    ; activate photon noise
  excess_noise      : 0b,                                             ; activate Excess noise of sqrt(2.) characteristic of EMCCDs.
  background_noise  : 0b,   ; activate sky background noise
  darkcurrent_noise : 0b,   ; activate dark current noise
  readout_noise     : 1b    ; activate detector readout noise
  readout_level     : 3.0   ; readout noise in e-/pix/frame
  darkcurrent_level : 0
  background_level  : 0
  sky_bg_norm       : 0
  quantum_eff       : 0.196   
}

{detector_ngs1,
  name: 'C-RED',                                                ; name of the detector (used by ccd auto_params_management method)
  size: [240,240],                                              ; Detector size in pixels
  dt: 0.002d,                                                   ; [s] Detector integration time
  bandw: 330,                                                   ; [nm] Sensor bandwidth
  binning: 1,                                                   ; detector binning (1,2,3...)
  photon_noise: 1b,                                             ; activate photon noise
  background_noise: 1b,                                         ; activate sky background noise
  darkcurrent_noise: 1b,                                        ; activate dark current noise
  excess_noise: 1b,                                             ; activate Excess noise of sqrt(2.) characteristic of EMCCDs.
  excess_delta: 3.33
  readout_noise: 1b,                                            ; activate readout noise
  readout_level: 0.5,                                           ; readout noise in e-/pix/frame (computed by vlt_lgs_sh_cloop)
  darkcurrent_level: 'auto',                                        ; dark current value in e-/pix/frame (computed by ccd auto_params_management method)
  background_level: 'auto',                                         ; sky backgroundt value in e-/pix/frame (computed by ccd auto_params_management method)
  sky_bg_norm: 2035,                                               ; [e-/m^2/arcsec^2/s] sky background value (used by ccd auto_params_management method)
  quantum_eff: 0.314,                                           ; detector quantum efficiency
  photon_seed: 1,                                               ; photon noise seed
  readout_seed: 1,                                              ; RON seed
  excess_seed: 1                                                ; excess noise seed
}

{detector_ngs2,
  name: 'C-RED',                                                ; name of the detector (used by ccd auto_params_management method)
  size: [240,240],                                              ; Detector size in pixels
  dt: 0.002d,                                                   ; [s] Detector integration time
  bandw: 330,                                                   ; [nm] Sensor bandwidth
  binning: 1,                                                   ; detector binning (1,2,3...)
  photon_noise: 1b,                                             ; activate photon noise
  background_noise: 1b,                                         ; activate sky background noise
  darkcurrent_noise: 1b,                                        ; activate dark current noise
  excess_noise: 1b,                                             ; activate Excess noise of sqrt(2.) characteristic of EMCCDs.
  excess_delta: 3.33
  readout_noise: 1b,                                            ; activate readout noise
  readout_level: 0.5,                                           ; readout noise in e-/pix/frame (computed by vlt_lgs_sh_cloop)
  darkcurrent_level: 'auto',                                        ; dark current value in e-/pix/frame (computed by ccd auto_params_management method)
  background_level: 'auto',                                         ; sky backgroundt value in e-/pix/frame (computed by ccd auto_params_management method)
  sky_bg_norm: 2035,                                               ; [e-/m^2/arcsec^2/s] sky background value (used by ccd auto_params_management method)
  quantum_eff: 0.314,                                           ; detector quantum efficiency
  photon_seed: 1,                                               ; photon noise seed
  readout_seed: 1,                                              ; RON seed
  excess_seed: 1                                                ; excess noise seed
}

{detector_ngs3,
  name: 'C-RED',                                                ; name of the detector (used by ccd auto_params_management method)
  size: [240,240],                                              ; Detector size in pixels
  dt: 0.002d,                                                   ; [s] Detector integration time
  bandw: 330,                                                   ; [nm] Sensor bandwidth
  binning: 1,                                                   ; detector binning (1,2,3...)
  photon_noise: 1b,                                             ; activate photon noise
  background_noise: 1b,                                         ; activate sky background noise
  darkcurrent_noise: 1b,                                        ; activate dark current noise
  excess_noise: 1b,                                             ; activate Excess noise of sqrt(2.) characteristic of EMCCDs.
  excess_delta: 3.33
  readout_noise: 1b,                                            ; activate readout noise
  readout_level: 0.5,                                           ; readout noise in e-/pix/frame (computed by vlt_lgs_sh_cloop)
  darkcurrent_level: 'auto',                                        ; dark current value in e-/pix/frame (computed by ccd auto_params_management method)
  background_level: 'auto',                                         ; sky backgroundt value in e-/pix/frame (computed by ccd auto_params_management method)
  sky_bg_norm: 2035,                                               ; [e-/m^2/arcsec^2/s] sky background value (used by ccd auto_params_management method)
  quantum_eff: 0.314,                                           ; detector quantum efficiency
  photon_seed: 1,                                               ; photon noise seed
  readout_seed: 1,                                              ; RON seed
  excess_seed: 1                                                ; excess noise seed
}

{slopec1,
   subapdata_tag = 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_rot6.2'
   sn_tag        = 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_ce_rot6.2'
   filtmat_tag   = 'maory_np_filtmat_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_mn1000_3_ce_rot6.2'
   use_sn        = 0B
   thr_value     = 0.0
}

{slopec2,
   subapdata_tag = 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_rot14.2'
   sn_tag        = 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_ce_rot14.2'
   filtmat_tag   = 'maory_np_filtmat_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_mn1000_3_ce_rot14.2'
   use_sn        = 0B
   thr_value     = 0.0
}

{slopec3,
   subapdata_tag = 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_rot-6.2'
   sn_tag        = 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_ce_rot-6.2'
   filtmat_tag   = 'maory_np_filtmat_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_mn1000_3_ce_rot-6.2'
   use_sn        = 0B
   thr_value     = 0.0
}

{slopec4,
   subapdata_tag = 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_rot6.2'
   sn_tag        = 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_ce_rot6.2'
   filtmat_tag   = 'maory_np_filtmat_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_mn1000_3_ce_rot6.2'
   use_sn        = 0B
   thr_value     = 0.0
}

{slopec5,
   subapdata_tag = 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_rot14.2'
   sn_tag        = 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_ce_rot14.2'
   filtmat_tag   = 'maory_np_filtmat_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_mn1000_3_ce_rot14.2'
   use_sn        = 0B
   thr_value     = 0.0
}

{slopec6,
   subapdata_tag = 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_rot-6.2'
   sn_tag        = 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_ce_rot-6.2'
   filtmat_tag   = 'maory_np_filtmat_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_mn1000_3_ce_rot-6.2'
   use_sn        = 0B
   thr_value     = 0.0
}

{slopec_ngs1,
   subapdata_tag = 'auto'
   sn_tag = 'auto'
   thr_value = 0.0                                         ; [ph] threshold value used in the slope computation
   weightedPixRad = 1.0
}

{slopec_ngs2,
   subapdata_tag = 'auto'
   sn_tag = 'auto'
   thr_value = 0.0                                         ; [ph] threshold value used in the slope computation
   weightedPixRad = 1.0
}

{slopec_ngs3,
   subapdata_tag = 'auto'
   sn_tag = 'auto'
   thr_value = 0.0                                         ; [ph] threshold value used in the slope computation
   weightedPixRad = 1.0
}

{modalrec1,
   intmat_tag: 'autoTN'
   recmat_tag: 'autoTN'
   projmat_tag: 'autoTN'
   nmodes: 4519L
   nmodes_popt: 6869L
   lgsTT: 1B
   sigma2innm2: 2e4
   noiseCovScaleVector_tag: 'maory_np_ps480p0.080_shs68x68_wl589_fv16.1_np14_th0.50_fluxPerSa'
   turb_var_from_dm: 1B
   polc: 1B
   polcRev: 1B
   proj_regFactor: 1e-4
   tag_ifunc4proj: 'MAORY_480pix_99nacts_cir_spider2023_5zern_10000.0cn'
   singlePrecision: 1B
   POLCwoLOmodes: 1B
   doNotPutOnGpu: 1B
   doNotReduceModesInProj: 1B
   noiseCovTallon: 1B
   naThicknessInM: 10e3;8e3; 
   tGparameter: 0.0
   onlyTallonDiag: 0B
}

{modalrec2,
   intmat_tag: 'autoTN'
   recmat_tag: 'autoTN'
   projmat_tag: 'autoTN'
   nmodes: 1021L
   lgsTT: 1B
   polc: 1B
   doNotPutOnGpu: 1B
}

{modalrec3,
   intmat_tag: 'autoTN'
   recmat_tag: 'autoTN'
   projmat_tag: 'autoTN'
   nmodes: 1024
   lgsTT: 1B
   polc: 1B
   doNotPutOnGpu: 1B
}

{modalrec_ngs1,
   intmat_tag: 'autoTN'
   recmat_tag: 'autoTN'
   projmat_tag: 'autoTN'
   nmodes: 2
   ;nmodesIM: 5
   lgsTT: 1B
   nphCoeff: -1
   directLO: 1B 
   polc: 1B
   proj_regFactor: 4e-3
   tag_ifunc4proj: 'MAORY_480pix_99nacts_cir_spider2023_5zern_10000.0cn'
   ;doNotReduceModesInProj = 1B
   calibPlateScaleModes = 1B
   plateScale = 1B
}

{modalrec_ngs2,
   intmat_tag: 'autoTN'
   recmat_tag: 'autoTN'
   projmat_tag: 'autoTN'
   nmodes: 0
   lgsTT: 1B
   polc: 1B
}

{modalrec_ngs3,
   intmat_tag: 'autoTN'
   recmat_tag: 'autoTN'
   projmat_tag: 'autoTN'
   nmodes: 3
   lgsTT: 1B
   polc: 1B
}

{control1,
   delay :   2                                                  ; Total temporal delay (relative to CCD integration time)
   type  : 'INT'                                                ; type of control ('INT', 'OMGI')
   int_gain: [0]         ; Integrator gain (for 'INT' control)
}

{control2,
   delay :   2                                                  ; Total temporal delay (relative to CCD integration time)
   type  : 'INT'                                                ; type of control ('INT', 'OMGI')
   int_gain: [0]          ; Integrator gain (for 'INT' control)
}

{control3,
   delay :   2                                                  ; Total temporal delay (relative to CCD integration time)
   type  : 'INT'                                                ; type of control ('INT', 'OMGI')
   int_gain: [0]          ; Integrator gain (for 'INT' control)
}

{source_lgs1,

  polar_coordinate: [45.0, 30], ; Source coordinates in [distance, angle] [arcsec, degrees]
                                 ; Zero distance is the telescope axis
  height: 90000.,                ; Source height [m]
  magnitude = 6.05,               ; Source magnitude
  wavelengthInNm: 589            ; Source wavelength SUPPORTED ?
}

{source_lgs2,

  polar_coordinate: [45.0, 90],   ; Source coordinates in [distance, angle] [arcsec, degrees]
                                    ; Zero distance is the telescope axis
  height: 90000.,                ; Source height [m]
  magnitude = 6.05,               ; Source magnitude
  wavelengthInNm: 589            ; Source wavelength SUPPORTED ?
}

{source_lgs3,

  polar_coordinate: [45.0, 150],   ; Source coordinates in [distance, angle] [arcsec, degrees]
                                    ; Zero distance is the telescope axis
  height: 90000.,                ; Source height [m]
  magnitude = 6.05,               ; Source magnitude
  wavelengthInNm: 589            ; Source wavelength SUPPORTED ?
}

{source_lgs4,

  polar_coordinate: [45.0, 210],   ; Source coordinates in [distance, angle] [arcsec, degrees]
                                    ; Zero distance is the telescope axis
  height: 90000.,                ; Source height [m]
  magnitude = 6.05,               ; Source magnitude
  wavelengthInNm: 589            ; Source wavelength SUPPORTED ?
}

{source_lgs5,

  polar_coordinate: [45.0, 270],   ; Source coordinates in [distance, angle] [arcsec, degrees]
                                    ; Zero distance is the telescope axis
  height: 90000.,                ; Source height [m]
  magnitude = 6.05,               ; Source magnitude
  wavelengthInNm: 589            ; Source wavelength SUPPORTED ?
}

{source_lgs6,

  polar_coordinate: [45.0, 330],   ; Source coordinates in [distance, angle] [arcsec, degrees]
                                    ; Zero distance is the telescope axis
  height: 90000.,                ; Source height [m]
  magnitude = 6.05,               ; Source magnitude
  wavelengthInNm: 589            ; Source wavelength SUPPORTED ?
}

{source_axis,

  polar_coordinate: [0.0, 45.0],    ; Source coordinates in [distance, angle] [arcsec, degrees]
                                    ; Zero distance is the telescope axis
  height: !VALUES.F_INFINITY,       ; Source height [m]
  magnitude = 8,                    ; Source magnitude
  wavelengthInNm: 750               ; Source wavelength SUPPORTED ?
}

{science_source1,

  polar_coordinate: [40.0, 45.0],    ; Source coordinates in [distance, angle] [arcsec, degrees]
                                    ; Zero distance is the telescope axis
  height: !VALUES.F_INFINITY,       ; Source height [m]
  magnitude = 8,                    ; Source magnitude
  wavelengthInNm: 750               ; Source wavelength SUPPORTED ?
}

{science_source2,

  polar_coordinate: [80.0, 45.0],    ; Source coordinates in [distance, angle] [arcsec, degrees]
                                    ; Zero distance is the telescope axis
  height: !VALUES.F_INFINITY,       ; Source height [m]
  magnitude = 8,                    ; Source magnitude
  wavelengthInNm: 750               ; Source wavelength SUPPORTED ?
}


{source_ngs1,
  polar_coordinate: [55.0, 0.0],    ; Source coordinates in [distance, angle] [arcsec, degrees]
                                    ; Zero distance is the telescope axis
  height: !VALUES.F_INFINITY,       ; Source height [m]
  magnitude       : 10,                 ; Source magnitude
  wavelengthInNm  : 1650                ; Source wavelength (must be the same as WFS)
  zeroPoint       : 1.14d-9
}
  
{source_ngs2,
  polar_coordinate: [55.0, 120.0],    ; Source coordinates in [distance, angle] [arcsec, degrees]
                                    ; Zero distance is the telescope axis
  height: !VALUES.F_INFINITY,       ; Source height [m]
  magnitude       : 10,                 ; Source magnitude
  wavelengthInNm  : 1650                ; Source wavelength (must be the same as WFS)
  zeroPoint       : 1.14d-9
}

{source_ngs3,
  polar_coordinate: [55.0, 240.0],    ; Source coordinates in [distance, angle] [arcsec, degrees]
                                    ; Zero distance is the telescope axis
  height: !VALUES.F_INFINITY,       ; Source height [m]
  magnitude       : 10,                 ; Source magnitude
  wavelengthInNm  : 1650                ; Source wavelength (must be the same as WFS)
  zeroPoint       : 1.14d-9
}

{camera,
  wavelengthInNm : 2200              ; [nm] Imaging wavelength
  nd             : 4,
}

{atmo,
  L0: 25,                         ; [m] Outer scale
  wavelengthInNm: 500,
  heights: [30.0000, 90.0000, 150.000, 200.000, 245.000, 300.000, 390.000, 600.000, 1130.00, 1880.00, 2630.00, 3500.00, 4500.00, 5500.00, 6500.00, 7500.00, 8500.00, 9500.00, 10500.0, 11500.0, 12500.0, 13500.0, 14500.0, 15500.0, 16500.0, 17500.0, 18500.0, 19500.0, 20500.0, 21500.0, 22500.0, 23500.0, 24500.0, 25500.0, 26500.0]      
  Cn2    : [0.241954, 0.119977, 0.0968817, 0.0589889, 0.0472911, 0.0472911, 0.0472911, 0.0472911, 0.0398925, 0.0323939, 0.0161969, 0.0260951, 0.0155971, 0.0103980, 0.00999811, 0.0119977, 0.00400924, 0.0139974, 0.0129975, 0.00700868, 0.0159970, 0.0258951, 0.0190964, 0.00986813, 0.00616883, 0.00400924, 0.00246953, 0.00215959, 0.00184965, 0.00135974, 0.00110979, 0.000616883, 0.000925825, 0.000493907, 0.000431918]
  mcao_fov : 160.
  seed: 1
}

{seeing,
 func_type = 'SIN'
 constant =  0.644
}

{wind_speed,
 func_type = 'SIN'
 constant = [5.5, 5.5, 5.1, 5.5, 5.6, 5.7, 5.8, 6.0, 6.5, 7.0, 7.5, 8.5, 9.5, 11.5, 17.5, 23.0, 26.0, 29.0, 32.0, 27.0, 22.0, 14.5, 9.5, 6.3, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.0] ; [m/s] Wind speed
}

{wind_direction,
 func_type = 'SIN'
 constant = [0, -180, 0, 0, 90, 180, 0, 0, 0, -180, 0, 0, -90, 0, 90, -180, 90, 0, -90, -90, 0, -90, 0, 0, 180, 180, 0, -180, 90, 0, 0, 180, -90, 90, -90]                 ; [degrees] Wind direction
}

{pupil_stop,
  pupil_mask_tag: 'EELT480pp0.0803m_obs0.283_spider2023'
}

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;for optimization of specific directions

{source_opt1,
polar_coordinate: [0.0, 0.0]
weight: 0.5
}

{source_opt2,
polar_coordinate: [30.0, 0.0]
weight: 1.0
}

{source_opt3,
polar_coordinate: [30.0, 45.0]
weight: 1.0
}

{source_opt4,
polar_coordinate: [30.0, 90.0]
weight: 1.0
}

{source_opt5,
polar_coordinate: [30.0, 135.0]
weight: 1.0
}

{source_opt6,
polar_coordinate: [30.0, 180.0]
weight: 1.0
}

{source_opt7,
polar_coordinate: [30.0, 225.0]
weight: 1.0
}

{source_opt8,
polar_coordinate: [30.0, 270.0]
weight: 1.0
}

{source_opt9,
polar_coordinate: [30.0, 315.0]
weight: 1.0
}

{source_opt10,
polar_coordinate: [80.0, 0.0]
weight: 0.05
}

{source_opt11,
polar_coordinate: [80.0, 45.0]
weight: 0.05
}

{source_opt12,
polar_coordinate: [80.0, 90.0]
weight: 0.05
}

{source_opt13,
polar_coordinate: [80.0, 135.0]
weight: 0.05
}

{source_opt14,
polar_coordinate: [80.0, 180.0]
weight: 0.05
}

{source_opt15,
polar_coordinate: [80.0, 225.0]
weight: 0.05
}

{source_opt16,
polar_coordinate: [80.0, 270.0]
weight: 0.05
}

{source_opt17,
polar_coordinate: [80.0, 315.0]
weight: 0.05
}

{layer1,
  ifunc_tag : 'MAORY_480pix_99nacts_cir_spider2023_zonal_ifs',
  m2c_tag   : 'MAORY_480pix_99nacts_cir_spider2023_5zern_10000.0cn',
  nmodes    : 5000L
  start_mode: 3                                             
  height    : 0   
  doNotPutOnGpu: 1B
}

{layer2,
  ifunc_tag : 'MAORY_486pix_81nacts_cir_0.000obs_zonal_ifs',
  m2c_tag   : 'MAORY_486pix_81nacts_cir_0.000obs_5zern_10000.0cn'
  nmodes    : 4920
  start_mode: 3                                              
  height    : 600                                            ; takes into account airmass
  doNotPutOnGpu: 1B
}

{layer3,
  ifunc_tag : 'MAORY_500pix_83nacts_cir_0.000obs_zonal_ifs'
  m2c_tag   : 'MAORY_500pix_83nacts_cir_0.000obs_5zern_10000.0cn'
  nmodes    : 5166
  start_mode: 3                                              
  height    : 2000                                           ; takes into account airmass
  doNotPutOnGpu: 1B
}

{layer4,
  ifunc_tag : 'MAORY_522pix_73nacts_cir_0.000obs_zonal_ifs'
  m2c_tag   : 'MAORY_522pix_73nacts_cir_0.000obs_5zern_10000.0cn'
  nmodes    : 3996
  start_mode: 3                                              
  height    : 4500                                           ; takes into account airmass
  doNotPutOnGpu: 1B
}

{layer5,
  ifunc_tag : 'MAORY_552pix_61nacts_cir_0.000obs_zonal_ifs' 
  m2c_tag   : 'MAORY_552pix_61nacts_cir_0.000obs_5zern_10000.0cn'
  nmodes    : 2790
  start_mode: 3                                              
  height    : 7500                                           ; takes into account airmass   
  doNotPutOnGpu: 1B
}

{layer6,
  ifunc_tag : 'MAORY_586pix_65nacts_cir_0.000obs_zonal_ifs'
  m2c_tag   : 'MAORY_586pix_65nacts_cir_0.000obs_5zern_10000.0cn'
  nmodes    : 3168
  start_mode: 3                                              
  height    : 11000                                           ; takes into account airmass
  doNotPutOnGpu: 1B
}

{layer7,
  ifunc_tag : 'MAORY_624pix_69nacts_cir_0.000obs_zonal_ifs'
  m2c_tag   : 'MAORY_624pix_69nacts_cir_0.000obs_5zern_10000.0cn'
  nmodes    : 3570
  start_mode: 3                                              
  height    : 15000                                           ; takes into account airmass
  doNotPutOnGpu: 1B
}

{layer8,
  ifunc_tag : 'MAORY_652pix_73nacts_cir_0.000obs_zonal_ifs'
  m2c_tag   : 'MAORY_652pix_73nacts_cir_0.000obs_5zern_10000.0cn'
  nmodes    : 3996
  start_mode: 3                                              
  height    : 18000                                           ; takes into account airmass
  doNotPutOnGpu: 1B
}

{layer9,
  ifunc_tag : 'MAORY_690pix_57nacts_cir_0.000obs_zonal_ifs'
  m2c_tag   : 'MAORY_690pix_57nacts_cir_0.000obs_5zern_10000.0cn'
  nmodes    : 2436
  start_mode: 3                                              
  height    : 22000                                           ; takes into account airmass
  doNotPutOnGpu: 1B
}

