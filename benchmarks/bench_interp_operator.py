"""
Benchmark: rotshiftzoom_array with interp='affine' (affine_transform) vs
interp='sparse' (precomputed bilinear sparse operator).

For each case it reports execution time, approximate peak memory and the
difference between the two results. It also checks the equivalence on the
selected device and exits with status 1 if the check fails.

Examples:
    # GPU 0, default sizes
    python benchmarks/bench_interp_operator.py --device 0

    # MORFEO-like M4 cube and end-to-end interaction matrix
    python benchmarks/bench_interp_operator.py --device 0 --sizes 480 \\
        --nmodes 4519 --im --im-nmodes 4519 --json bench_l40s.json

    # CPU
    python benchmarks/bench_interp_operator.py --device -1 --sizes 240 --nmodes 300
"""
import argparse
import json
import os
import platform
import sys
import time
import tracemalloc

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--device', type=int, default=0,
                        help='GPU index, -1 for CPU (default: 0)')
    parser.add_argument('--precision', type=int, default=1,
                        help='0 double, 1 single (default: 1)')
    parser.add_argument('--sizes', type=int, nargs='+', default=[480, 690],
                        help='grid sizes in pixels (default: 480 690)')
    parser.add_argument('--nmodes', type=int, default=1000,
                        help='number of modes of the test cube (default: 1000)')
    parser.add_argument('--repeat', type=int, default=3,
                        help='repetitions, the best time is reported (default: 3)')
    parser.add_argument('--im', action='store_true',
                        help='also time an end-to-end interaction matrix')
    parser.add_argument('--im-npix', type=int, default=480)
    parser.add_argument('--im-nmodes', type=int, default=1000)
    parser.add_argument('--im-nsubaps', type=int, default=68)
    parser.add_argument('--json', type=str, default=None,
                        help='save the results to this JSON file')
    return parser.parse_args()


args = parse_args()

# synim.init must be called before importing the synim submodules
import synim  # noqa: E402
synim.init(device_idx=args.device, precision=args.precision)
import synim.utils as synim_utils  # noqa: E402
import synim.synim as synim_core  # noqa: E402
from synim import cpuArray  # noqa: E402

sys.path.insert(0, os.path.join(HERE, '..', 'test'))
from reference_rotshiftzoom import rotshiftzoom_array_reference  # noqa: E402

xp = synim.xp
ON_GPU = xp is not np
if ON_GPU:
    import cupy as cp

# LGS-like geometry: DM shift and magnification, rotated and shifted WFS
GEOMETRY = dict(dm_translation=(3.7, -1.2), dm_rotation=0.5,
                dm_magnification=(1.003, 1.003), wfs_translation=(0.4, -0.3),
                wfs_rotation=30.0, wfs_magnification=(1.0, 1.0))

# Equivalence tolerance relative to max(|input|). On GPU the affine kernel may
# compute in single precision and may round coordinates differently at the
# grid border, so pixels whose input position is within BORDER_EPS of the
# border are reported separately.
RTOL = 1e-5 if args.precision == 1 else 1e-10
BORDER_EPS = 1e-9


def synchronize():
    if ON_GPU:
        cp.cuda.Device().synchronize()


def free_memory():
    if ON_GPU:
        cp.get_default_memory_pool().free_all_blocks()


def measure(func):
    """Run func once, return (result, seconds, peak additional bytes)."""
    free_memory()
    synchronize()
    if ON_GPU:
        pool = cp.get_default_memory_pool()
        used_before = pool.used_bytes()
        t0 = time.perf_counter()
        result = func()
        synchronize()
        elapsed = time.perf_counter() - t0
        # The pool keeps freed blocks, so its total size is the high-water mark
        peak = pool.total_bytes() - used_before
    else:
        tracemalloc.start()
        tracemalloc.reset_peak()
        t0 = time.perf_counter()
        result = func()
        elapsed = time.perf_counter() - t0
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    return result, elapsed, peak


def best_of(func, repeat):
    times, peaks, result = [], [], None
    for _ in range(repeat):
        result = None
        result, elapsed, peak = measure(func)
        times.append(elapsed)
        peaks.append(peak)
    return result, min(times), max(peaks)


class CaptureOperatorGeometry:
    """Record the (matrix, offset) passed to build_bilinear_operator."""

    def __init__(self):
        self.calls = []
        self._original = synim_utils.build_bilinear_operator

    def __enter__(self):
        def wrapper(input_shape, output_shape, matrix, offset, **kwargs):
            self.calls.append((input_shape, output_shape,
                               np.asarray(cpuArray(matrix), dtype=np.float64),
                               np.asarray(cpuArray(offset), dtype=np.float64)))
            return self._original(input_shape, output_shape, matrix, offset, **kwargs)
        synim_utils.build_bilinear_operator = wrapper
        return self

    def __exit__(self, *exc):
        synim_utils.build_bilinear_operator = self._original


def border_ambiguous_mask(input_shape, output_shape, matrix, offset):
    ii = np.arange(output_shape[0], dtype=np.float64)[:, None]
    jj = np.arange(output_shape[1], dtype=np.float64)[None, :]
    y = (offset[0] + ii * matrix[0, 0]) + jj * matrix[0, 1]
    x = (offset[1] + ii * matrix[1, 0]) + jj * matrix[1, 1]
    near = np.zeros(y.shape, dtype=bool)
    for coord, n in ((y, input_shape[0]), (x, input_shape[1])):
        near |= np.abs(coord) < BORDER_EPS
        near |= np.abs(coord - (n - 1)) < BORDER_EPS
    return near


def compare(result, reference, scale, ambiguous=None):
    diff = np.abs(cpuArray(result).astype(np.float64) - cpuArray(reference).astype(np.float64))
    info = {'max_rel_diff': float(diff.max() / scale)}
    if ambiguous is not None and diff.ndim == 3:
        info['n_border_ambiguous_pixels'] = int(ambiguous.sum())
        info['max_rel_diff_excluding_border'] = float(diff[~ambiguous].max() / scale) \
            if (~ambiguous).any() else 0.0
    else:
        info['max_rel_diff_excluding_border'] = info['max_rel_diff']
    return info


def environment():
    env = {'python': platform.python_version(), 'numpy': np.__version__,
           'platform': platform.platform(), 'device': args.device,
           'precision': args.precision}
    import scipy
    env['scipy'] = scipy.__version__
    if ON_GPU:
        env['cupy'] = cp.__version__
        props = cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)
        name = props['name']
        env['gpu'] = name.decode() if isinstance(name, bytes) else str(name)
        env['gpu_memory_gb'] = cp.cuda.Device().mem_info[1] / 1e9
    return env


def bench_rotshiftzoom(npix, nmodes, layout):
    rng = np.random.default_rng(0)
    host = rng.standard_normal((npix, npix, nmodes)).astype(synim.float_dtype)
    cube = xp.asarray(host)
    del host
    if layout == 'transposed':
        # As in the pipeline with the SPECULA convention
        cube = xp.transpose(cube, (1, 0, 2))
    cube_bytes = cube.nbytes
    scale = float(xp.max(xp.abs(cube)))

    out_affine, t_affine, m_affine = best_of(
        lambda: synim_utils.rotshiftzoom_array(cube, interp='affine', **GEOMETRY),
        args.repeat)
    out_affine = cpuArray(out_affine)
    free_memory()

    with CaptureOperatorGeometry() as capture:
        out_sparse, t_sparse, m_sparse = best_of(
            lambda: synim_utils.rotshiftzoom_array(cube, interp='sparse', **GEOMETRY),
            args.repeat)
    out_sparse = cpuArray(out_sparse)
    free_memory()

    # Operator construction alone (CPU), to separate it from the product
    in_shape, out_shape, matrix, offset = capture.calls[-1]
    t0 = time.perf_counter()
    operator = synim_utils.build_bilinear_operator(in_shape, out_shape, matrix, offset,
                                                   dtype=synim.float_dtype)
    t_build = time.perf_counter() - t0
    operator_mb = (operator.data.nbytes + operator.indices.nbytes
                   + operator.indptr.nbytes) / 1e6

    ambiguous = border_ambiguous_mask(in_shape, out_shape, matrix, offset)
    equivalence = compare(out_sparse, out_affine, scale, ambiguous)

    # On the first (smallest) slices also compare with the frozen reference
    n_ref = min(nmodes, 8)
    ref = cpuArray(rotshiftzoom_array_reference(cube[:, :, :n_ref], **GEOMETRY))
    equivalence['max_rel_diff_vs_frozen_reference'] = compare(
        out_sparse[:, :, :n_ref], ref, scale, ambiguous)['max_rel_diff_excluding_border']
    del cube, out_affine, out_sparse, ref
    free_memory()

    return {
        'case': 'rotshiftzoom', 'npix': npix, 'nmodes': nmodes, 'layout': layout,
        'cube_gb': cube_bytes / 1e9,
        'time_affine_s': t_affine, 'time_sparse_s': t_sparse,
        'time_sparse_build_s': t_build, 'speedup': t_affine / t_sparse,
        'peak_affine_gb': m_affine / 1e9, 'peak_sparse_gb': m_sparse / 1e9,
        'peak_affine_cubes': m_affine / cube_bytes, 'peak_sparse_cubes': m_sparse / cube_bytes,
        'operator_mb': operator_mb,
        **equivalence,
    }


def bench_interaction_matrix(npix, nmodes, nsubaps):
    yy, xx = np.mgrid[:npix, :npix]
    r = np.hypot(xx - npix / 2 + 0.5, yy - npix / 2 + 0.5)
    pup_mask = ((r < npix / 2 - 2) & (r > 0.28 * npix / 2)).astype(np.float32)
    dm_mask = (r < npix / 2).astype(np.float32)
    rng = np.random.default_rng(1)
    dm_array = xp.asarray(rng.standard_normal((npix, npix, nmodes)).astype(synim.float_dtype))
    cube_bytes = dm_array.nbytes
    kwargs = dict(pup_diam_m=38.5, pup_mask=pup_mask, dm_array=dm_array, dm_mask=dm_mask,
                  dm_height=0.0, dm_rotation=0.0, wfs_nsubaps=nsubaps, wfs_fov_arcsec=16.0,
                  gs_pol_coo=(45.0, 30.0), gs_height=90e3, wfs_rotation=30.0,
                  wfs_translation=(0.2, -0.1), wfs_mag_global=1.0)
    results = {}
    for method in ('affine', 'sparse'):
        synim_utils.set_interp_method(method)
        try:
            im, t, m = best_of(lambda: synim_core.interaction_matrix(**kwargs), args.repeat)
        finally:
            synim_utils.set_interp_method('affine')
        results[method] = (cpuArray(im), t, m)
        del im
        free_memory()
    del dm_array
    free_memory()
    im_a, t_a, m_a = results['affine']
    im_s, t_s, m_s = results['sparse']
    return {
        'case': 'interaction_matrix', 'npix': npix, 'nmodes': nmodes, 'nsubaps': nsubaps,
        'cube_gb': cube_bytes / 1e9,
        'time_affine_s': t_a, 'time_sparse_s': t_s, 'speedup': t_a / t_s,
        'peak_affine_gb': m_a / 1e9, 'peak_sparse_gb': m_s / 1e9,
        'peak_affine_cubes': m_a / cube_bytes, 'peak_sparse_cubes': m_s / cube_bytes,
        'max_rel_diff': float(np.max(np.abs(im_a - im_s)) / np.max(np.abs(im_a))),
    }


def main():
    env = environment()
    print('Environment:', json.dumps(env))
    results = []
    for npix in args.sizes:
        for layout in ('contiguous', 'transposed'):
            print(f'\nrotshiftzoom: {npix}x{npix}x{args.nmodes}, {layout} input ...', flush=True)
            res = bench_rotshiftzoom(npix, args.nmodes, layout)
            results.append(res)
            print(f"  affine {res['time_affine_s']:.3f} s, peak {res['peak_affine_gb']:.2f} GB"
                  f" ({res['peak_affine_cubes']:.1f} cubes)")
            print(f"  sparse {res['time_sparse_s']:.3f} s (build {res['time_sparse_build_s']:.3f} s),"
                  f" peak {res['peak_sparse_gb']:.2f} GB ({res['peak_sparse_cubes']:.1f} cubes),"
                  f" operator {res['operator_mb']:.1f} MB")
            print(f"  speed-up {res['speedup']:.2f}x, max rel diff {res['max_rel_diff']:.2e}"
                  f" (excluding {res['n_border_ambiguous_pixels']} border-ambiguous pixels:"
                  f" {res['max_rel_diff_excluding_border']:.2e}),"
                  f" vs frozen reference {res['max_rel_diff_vs_frozen_reference']:.2e}")
    if args.im:
        print(f'\ninteraction_matrix: {args.im_npix}px, {args.im_nmodes} modes,'
              f' {args.im_nsubaps}x{args.im_nsubaps} subapertures ...', flush=True)
        res = bench_interaction_matrix(args.im_npix, args.im_nmodes, args.im_nsubaps)
        results.append(res)
        print(f"  affine {res['time_affine_s']:.3f} s, peak {res['peak_affine_gb']:.2f} GB"
              f" ({res['peak_affine_cubes']:.1f} cubes)")
        print(f"  sparse {res['time_sparse_s']:.3f} s, peak {res['peak_sparse_gb']:.2f} GB"
              f" ({res['peak_sparse_cubes']:.1f} cubes)")
        print(f"  speed-up {res['speedup']:.2f}x, IM max rel diff {res['max_rel_diff']:.2e}")

    failed = [r for r in results if r['max_rel_diff_excluding_border' if r['case'] ==
                                     'rotshiftzoom' else 'max_rel_diff'] > RTOL]
    if args.json:
        with open(args.json, 'w') as f:
            json.dump({'environment': env, 'rtol': RTOL, 'results': results}, f, indent=2)
        print(f'\nResults saved to {args.json}')
    if failed:
        print(f'\nEQUIVALENCE CHECK FAILED (tolerance {RTOL:.0e}) for {len(failed)} case(s)')
        sys.exit(1)
    print(f'\nEquivalence check passed (tolerance {RTOL:.0e})')


if __name__ == '__main__':
    main()
