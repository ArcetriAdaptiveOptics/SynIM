"""
SynIM performance benchmark on the MORFEO configuration.

Tests (GPU and CPU run in separate processes, each initialized before
importing the SynIM submodules):
  1. single LGS interaction matrix (lgs1 x dm1) on GPU
  2. the same on CPU
  3. compute_interaction_matrices for all the LGS WFS and all the
     components (DMs and layers) on GPU, I/O included

For the single IM the first call (loading of the influence functions and
computation) and the second call (computation only: the component is then
in the ParamsManager cache) are timed separately, on both devices.
compute_interaction_matrices is run --multi-repeat times: the first run also
reads the FITS files from disk (cold file cache), the last one is reported.

GPU memory: peak of the memory in use in the cupy pool (MemoryHook on every
malloc/free) and peak of the device memory (sampled every 20 ms, as
nvidia-smi; includes the cupy pool blocks kept for reuse and the CUDA context).

Example (run it on each SynIM branch with the same arguments):
    python synim_benchmark.py --device 0 --json bench_perfMemory.json
"""
import argparse
import json
import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

DEFAULT_YAML = '/home/guido/pythonLib/SPECULA_scripts/morfeo/params_morfeo_calib.yml'
DEFAULT_ROOT = '/raid1/guido/PASSATA/MAORYC/'


def _init_libs(device_idx, precision):
    # SynIM and SPECULA bind the backend at import time: each device needs
    # its own process, initialized before importing the submodules
    import specula
    specula.init(device_idx, precision=precision)
    import synim
    synim.init(device_idx, precision=precision)
    from synim.params_manager import ParamsManager
    return ParamsManager


def _synim_version():
    import synim
    path = os.path.dirname(os.path.abspath(synim.__file__))
    try:
        commit = subprocess.run(['git', '-C', path, 'rev-parse', '--short', 'HEAD'],
                                capture_output=True, text=True, timeout=10).stdout.strip()
        branch = subprocess.run(['git', '-C', path, 'rev-parse', '--abbrev-ref', 'HEAD'],
                                capture_output=True, text=True, timeout=10).stdout.strip()
        dirty = subprocess.run(['git', '-C', path, 'status', '--porcelain', '--untracked-files=no'],
                               capture_output=True, text=True, timeout=10).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        commit, branch, dirty = '', '', ''
    return {'path': path, 'commit': commit, 'branch': branch, 'modified': bool(dirty)}


class GpuMemoryMonitor:
    """Peak GPU memory: in use in the cupy pool (hook) and on the device (sampled)."""

    def __init__(self, device_idx, interval=0.02):
        import cupy as cp
        self.cp = cp
        self.device_idx = device_idx
        self.interval = interval
        self.used = 0
        self.peak_in_use = 0
        self.peak_device = 0
        self._stop = threading.Event()
        monitor = self

        class _Hook(cp.cuda.MemoryHook):
            name = 'BenchmarkMemoryHook'

            def malloc_postprocess(self, **kwargs):
                monitor.used += kwargs['mem_size']
                if monitor.used > monitor.peak_in_use:
                    monitor.peak_in_use = monitor.used

            def free_postprocess(self, **kwargs):
                monitor.used -= kwargs['mem_size']

        self._hook = _Hook()

    def _sample(self):
        with self.cp.cuda.Device(self.device_idx):
            while not self._stop.is_set():
                free, total = self.cp.cuda.runtime.memGetInfo()
                self.peak_device = max(self.peak_device, total - free)
                time.sleep(self.interval)

    def __enter__(self):
        self.cp.cuda.Device(self.device_idx).synchronize()
        self.used = self.cp.get_default_memory_pool().used_bytes()
        self.peak_in_use = self.used
        self.peak_device = 0
        self._stop.clear()
        self._thread = threading.Thread(target=self._sample, daemon=True)
        self._thread.start()
        self._hook.__enter__()
        return self

    def __exit__(self, *exc):
        self.cp.cuda.Device(self.device_idx).synchronize()
        self._hook.__exit__(*exc)
        self._stop.set()
        self._thread.join()
        return False

    def result(self):
        return {'gpu_peak_in_use_gb': self.peak_in_use / 1e9,
                'gpu_peak_device_gb': self.peak_device / 1e9}


def _sync(device_idx):
    if device_idx >= 0:
        import cupy as cp
        cp.cuda.Device(device_idx).synchronize()


def _timed(func, device_idx):
    _sync(device_idx)
    t0 = time.perf_counter()
    result = func()
    _sync(device_idx)
    return result, time.perf_counter() - t0


def _single_im(pm, slope_method, device_idx, label):
    """First call (loading + computation) and second call (computation only)."""
    def call():
        return pm.compute_interaction_matrix(wfs_type='lgs', wfs_index=1, dm_index=1,
                                             slope_method=slope_method)
    im, t_first = _timed(call, device_idx)
    im, t_second = _timed(call, device_idx)
    print(f"  {label} single IM: first call (loading + computation) {t_first:.2f} s,"
          f" second call (computation) {t_second:.2f} s, shape {im.shape}, {im.dtype}")
    sys.stdout.flush()
    return {'t_single_first': t_first, 't_single_second': t_second,
            'im_shape': list(im.shape)}


def _run_gpu(yaml_file, root_dir, slope_method, precision, device_idx, multi_repeat):
    """Tests 1 and 3 on GPU."""
    ParamsManager = _init_libs(device_idx, precision)
    import synim
    result = {'synim': _synim_version()}

    print("\n[1/3] Single IM on GPU...")
    pm = ParamsManager(yaml_file, root_dir=root_dir, verbose=False)
    with GpuMemoryMonitor(device_idx) as monitor:
        result.update(_single_im(pm, slope_method, device_idx, 'GPU'))
    result['single'] = monitor.result()
    result['n_lgs'] = pm._count_wfs('lgs')
    del pm

    print(f"\n[3/3] compute_interaction_matrices, all LGS x all components"
          f" ({multi_repeat} run(s), the last one is reported)...")
    out_dir = os.path.join(root_dir, 'benchmark_temp')
    os.makedirs(out_dir, exist_ok=True)
    runs = []
    for i in range(multi_repeat):
        import cupy as cp
        cp.get_default_memory_pool().free_all_blocks()
        pm = ParamsManager(yaml_file, root_dir=root_dir, verbose=False)
        with GpuMemoryMonitor(device_idx) as monitor:
            paths, t_multi = _timed(lambda: pm.compute_interaction_matrices(
                output_im_dir=out_dir, output_rec_dir=out_dir, wfs_type='lgs',
                slope_method=slope_method, overwrite=True, verbose=False), device_idx)
        runs.append(dict(t_multi=t_multi, n_im=len(paths), **monitor.result()))
        print(f"  run {i + 1}: {t_multi:.1f} s for {len(paths)} IMs"
              f" ({t_multi / max(len(paths), 1):.2f} s/IM), GPU peak in use"
              f" {runs[-1]['gpu_peak_in_use_gb']:.1f} GB, device peak"
              f" {runs[-1]['gpu_peak_device_gb']:.1f} GB")
        sys.stdout.flush()
        del pm
    result['multi_runs'] = runs
    result['cupy'] = synim.cp.__version__ if synim.cp is not None else None
    return result


def _run_cpu(yaml_file, root_dir, slope_method, precision):
    """Test 2 on CPU."""
    ParamsManager = _init_libs(-1, precision)
    print("\n[2/3] Single IM on CPU...")
    pm = ParamsManager(yaml_file, root_dir=root_dir, verbose=False)
    return _single_im(pm, slope_method, -1, 'CPU')


def _in_subprocess(func, *args):
    # 'spawn': a fresh interpreter, no SynIM/SPECULA state inherited
    with ProcessPoolExecutor(max_workers=1, mp_context=get_context('spawn')) as ex:
        return ex.submit(func, *args).result()


def run_benchmark(yaml_file, root_dir, device_idx=0, slope_method='telsum', precision=1,
                  multi_repeat=2, skip_cpu=False, json_path=None):
    print("=" * 60)
    print("SynIM PERFORMANCE BENCHMARK (MORFEO ELT CONFIGURATION)")
    print("=" * 60)
    print(f"slope method: {slope_method}, precision: {precision}, GPU: {device_idx}")

    gpu = _in_subprocess(_run_gpu, yaml_file, root_dir, slope_method, precision,
                         device_idx, multi_repeat)
    cpu = None if skip_cpu else _in_subprocess(_run_cpu, yaml_file, root_dir,
                                               slope_method, precision)

    last = gpu['multi_runs'][-1]
    print("\n" + "=" * 60)
    print(f"SynIM: {gpu['synim']['branch']} {gpu['synim']['commit']}"
          f"{' (modified)' if gpu['synim']['modified'] else ''}  [{gpu['synim']['path']}]")
    print(f"Single IM {gpu['im_shape']}, computation only (second call):")
    line = f"  GPU {gpu['t_single_second']:.2f} s"
    if cpu is not None:
        line += (f" | CPU {cpu['t_single_second']:.2f} s"
                 f" | speed-up {cpu['t_single_second'] / gpu['t_single_second']:.0f}x")
    print(line)
    print("Single IM, loading + computation (first call):")
    line = f"  GPU {gpu['t_single_first']:.2f} s"
    if cpu is not None:
        line += f" | CPU {cpu['t_single_first']:.2f} s"
    print(line)
    print(f"  GPU memory: peak in use {gpu['single']['gpu_peak_in_use_gb']:.1f} GB,"
          f" device peak {gpu['single']['gpu_peak_device_gb']:.1f} GB")
    print(f"compute_interaction_matrices ({last['n_im']} IMs, {gpu['n_lgs']} LGS, I/O included):")
    for i, run in enumerate(gpu['multi_runs']):
        tag = ' (may include reading the FITS files from disk)' if i == 0 and len(gpu['multi_runs']) > 1 else ''
        print(f"  run {i + 1}: {run['t_multi']:.1f} s{tag}")
    print(f"  reported: {last['t_multi']:.1f} s, {last['t_multi'] / max(last['n_im'], 1):.2f} s/IM,"
          f" GPU peak in use {last['gpu_peak_in_use_gb']:.1f} GB,"
          f" device peak {last['gpu_peak_device_gb']:.1f} GB")
    print("=" * 60)

    if json_path:
        with open(json_path, 'w') as f:
            json.dump({'gpu': gpu, 'cpu': cpu, 'slope_method': slope_method,
                       'precision': precision, 'yaml': yaml_file}, f, indent=2)
        print(f"Results saved to {json_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--yaml', default=DEFAULT_YAML)
    parser.add_argument('--root-dir', default=DEFAULT_ROOT)
    parser.add_argument('--device', type=int, default=0)
    parser.add_argument('--slope-method', default='telsum')
    parser.add_argument('--precision', type=int, default=1, help='0 double, 1 single')
    parser.add_argument('--multi-repeat', type=int, default=2,
                        help='runs of compute_interaction_matrices (the last is reported)')
    parser.add_argument('--skip-cpu', action='store_true')
    parser.add_argument('--json', default=None, help='save the results to this file')
    args = parser.parse_args()
    run_benchmark(args.yaml, args.root_dir, device_idx=args.device,
                  slope_method=args.slope_method, precision=args.precision,
                  multi_repeat=args.multi_repeat, skip_cpu=args.skip_cpu,
                  json_path=args.json)
