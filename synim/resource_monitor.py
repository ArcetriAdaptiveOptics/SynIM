import os
import time
import threading


class ResourceMonitor:
    """
    Tracks peak CPU %, RAM and GPU memory usage per named section.

    A background thread samples at regular intervals; bracket operations with
    start_section() / end_section() or with the section() context manager.
    Call report() at any point to print a summary table.

    Usage::

        monitor = ResourceMonitor()
        with monitor.section("compute_interaction_matrices"):
            ...
        monitor.report()

        # Or integrated via ParamsManager:
        mgr = ParamsManager(yaml_file, monitor=True)
        mgr.compute_tomographic_reconstructor(...)
        mgr.resource_report()
    """

    def __init__(self, interval_ms=200):
        try:
            import psutil
        except ImportError:
            raise ImportError(
                "ResourceMonitor requires psutil.  Install with: pip install psutil"
            )
        self._psutil = psutil
        self._interval = interval_ms / 1000.0
        self._sections = []         # completed sections
        self._current = None        # section being sampled right now
        self._lock = threading.Lock()
        self._stop_event = threading.Event()
        self._process = psutil.Process(os.getpid())
        self._process.cpu_percent(interval=None)  # warm-up; first call always returns 0
        self._has_gpu = self._detect_gpu()
        self._thread = threading.Thread(target=self._sample_loop, daemon=True)
        self._thread.start()

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _detect_gpu(self):
        try:
            import cupy as cp
            cp.get_default_memory_pool().used_bytes()
            return True
        except Exception:
            return False

    def _sample_loop(self):
        while not self._stop_event.is_set():
            with self._lock:
                if self._current is not None:
                    try:
                        cpu = self._process.cpu_percent(interval=None)
                        self._current['peak_cpu'] = max(self._current['peak_cpu'], cpu)
                        ram = self._process.memory_info().rss
                        self._current['peak_ram'] = max(self._current['peak_ram'], ram)
                    except Exception:
                        pass
                    if self._has_gpu:
                        try:
                            import cupy as cp
                            used = cp.get_default_memory_pool().used_bytes()
                            self._current['peak_gpu'] = max(self._current['peak_gpu'], used)
                        except Exception:
                            pass
            time.sleep(self._interval)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def start_section(self, name):
        """Begin tracking a named section."""
        with self._lock:
            self._current = {
                'name': name,
                't_start': time.perf_counter(),
                'peak_cpu': 0.0,
                'peak_ram': 0,
                'peak_gpu': 0,
            }

    def end_section(self):
        """End the current section and store its peak statistics."""
        with self._lock:
            if self._current is None:
                return
            self._current['elapsed'] = time.perf_counter() - self._current['t_start']
            self._sections.append(self._current)
            self._current = None

    def section(self, name):
        """Return a context manager that brackets a named section."""
        return _SectionContext(self, name)

    def report(self):
        """Print a formatted table of peak resource usage per completed section."""
        with self._lock:
            sections = list(self._sections)

        if not sections:
            print("ResourceMonitor: no sections recorded yet.")
            return

        name_w = max(len(s['name']) for s in sections)
        name_w = max(name_w, 42)
        gpu_col = self._has_gpu

        header = (f"{'Section':<{name_w}}  {'Time (s)':>9}  {'Peak RAM':>9}")
        if gpu_col:
            header += f"  {'Peak GPU':>9}"
        sep = "─" * len(header)

        print(f"\n{'═' * len(header)}")
        print("  Resource usage summary")
        print(f"{'═' * len(header)}")
        print(header)
        print(sep)
        for s in sections:
            ram_gb = s['peak_ram'] / 1024 ** 3
            line = f"{s['name']:<{name_w}}  {s['elapsed']:>9.1f}  {ram_gb:>8.2f}G"
            if gpu_col:
                gpu_gb = s['peak_gpu'] / 1024 ** 3
                line += f"  {gpu_gb:>8.2f}G"
            print(line)
        print(f"{'═' * len(header)}\n")

    def reset(self):
        """Clear all recorded sections."""
        with self._lock:
            self._sections.clear()
            self._current = None

    def stop(self):
        """Stop the background sampling thread (optional; daemon thread dies with process)."""
        self._stop_event.set()
        self._thread.join()


# ---------------------------------------------------------------------------
# Context manager helper
# ---------------------------------------------------------------------------

class _SectionContext:
    def __init__(self, monitor, name):
        self._monitor = monitor
        self._name = name

    def __enter__(self):
        self._monitor.start_section(self._name)
        return self

    def __exit__(self, *_):
        self._monitor.end_section()


# ---------------------------------------------------------------------------
# No-op monitor used when monitoring is disabled (null-object pattern)
# ---------------------------------------------------------------------------

class _NullMonitor:
    """Drop-in replacement for ResourceMonitor when monitoring is disabled."""

    def start_section(self, name):
        pass

    def end_section(self):
        pass

    def section(self, name):
        return _NullSectionContext()

    def report(self):
        print("Resource monitoring is not enabled."
              " Pass monitor=True to ParamsManager to activate it.")

    def reset(self):
        pass


class _NullSectionContext:
    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass
