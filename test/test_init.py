"""
Tests of synim.init() and of the import-order checks.

The configuration is global to the Python process, so every case runs in a
fresh interpreter.
"""
import os
import subprocess
import sys
import textwrap
import unittest


def run_python(code):
    env = dict(os.environ)
    env['PYTHONPATH'] = os.pathsep.join(p for p in sys.path if p)
    env['SYNIM_DISABLE_GPU'] = 'TRUE'   # deterministic: CPU only
    env['MPLBACKEND'] = 'Agg'
    return subprocess.run([sys.executable, '-c', textwrap.dedent(code)],
                          capture_output=True, text=True, env=env, timeout=300)


class TestInit(unittest.TestCase):

    def assert_ok(self, result):
        self.assertEqual(result.returncode, 0, msg=result.stderr[-2000:])

    def test_import_without_init_raises(self):
        for module in ('synim.synim', 'synim.synpm', 'synim.utils',
                       'synim.params_utils', 'synim.params_manager'):
            with self.subTest(module=module):
                result = run_python(f"""
                    try:
                        import {module}
                    except RuntimeError as err:
                        assert 'synim.init' in str(err), err
                        print('raised')
                """)
                self.assert_ok(result)
                self.assertIn('raised', result.stdout)

    def test_package_import_does_not_initialize(self):
        result = run_python("""
            import synim
            assert synim.xp is None and synim.float_dtype is None
            # helpers that do not depend on the configuration still work
            import numpy as np
            assert synim.cpuArray(np.ones(3)).sum() == 3
            # the registration subpackage does not depend on the configuration
            import synim.registration.geometry
        """)
        self.assert_ok(result)

    def test_init_then_import(self):
        result = run_python("""
            import numpy as np
            import synim
            synim.init(device_idx=-1, precision=0)
            import synim.synim, synim.utils
            assert synim.synim.float_dtype is np.float64
            assert synim.utils.float_dtype is np.float64
            assert synim.synim.xp is np
        """)
        self.assert_ok(result)

    def test_reinit_before_import_is_allowed(self):
        result = run_python("""
            import numpy as np
            import synim
            synim.init(device_idx=-1, precision=1)
            synim.init(device_idx=-1, precision=0)
            import synim.utils
            assert synim.utils.float_dtype is np.float64
        """)
        self.assert_ok(result)

    def test_reinit_after_import_same_config(self):
        result = run_python("""
            import numpy as np
            import synim
            synim.init(device_idx=-1, precision=1)
            import synim.synim
            synim.init(device_idx=-1, precision=1)      # no effect, no error
            # GPU requested but not available: same effective configuration
            synim.init(device_idx=0, precision=1)
            assert synim.default_target_device_idx == -1
            assert synim.synim.float_dtype is np.float32
        """)
        self.assert_ok(result)

    def test_reinit_after_import_different_config_raises(self):
        result = run_python("""
            import numpy as np
            import synim
            synim.init(device_idx=-1, precision=1)
            import synim.synim
            try:
                synim.init(device_idx=-1, precision=0)
            except RuntimeError as err:
                assert 'synim.synim' in str(err), err
                print('raised')
            # the configuration is unchanged
            assert synim.float_dtype is np.float32 and synim.global_precision == 1
        """)
        self.assert_ok(result)
        self.assertIn('raised', result.stdout)

    def test_invalid_precision(self):
        result = run_python("""
            import synim
            try:
                synim.init(device_idx=-1, precision=2)
            except ValueError:
                print('raised')
        """)
        self.assert_ok(result)
        self.assertIn('raised', result.stdout)

    def test_specula_configuration_is_respected(self):
        # Importing params_utils must not re-initialize an already initialized SPECULA
        result = run_python("""
            import specula
            specula.init(device_idx=-1, precision=0)
            import synim
            synim.init(device_idx=-1, precision=1)
            import synim.params_utils
            assert specula.global_precision == 0, specula.global_precision
        """)
        self.assert_ok(result)


if __name__ == '__main__':
    unittest.main()
