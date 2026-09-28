"""
The messages printed by SynIM must be ASCII only: on Windows, with stdout
redirected (e.g. cp1252 encoding), non-ASCII characters raise
UnicodeEncodeError.
"""
import ast
import glob
import os
import unittest

PACKAGE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'synim')


class TestAsciiOutput(unittest.TestCase):

    def test_print_calls_are_ascii(self):
        offending = []
        for path in sorted(glob.glob(os.path.join(PACKAGE_DIR, '**', '*.py'), recursive=True)):
            with open(path, encoding='utf-8') as f:
                source = f.read()
            for node in ast.walk(ast.parse(source)):
                if isinstance(node, ast.Call) and getattr(node.func, 'id', None) == 'print':
                    segment = ast.get_source_segment(source, node)
                    if any(ord(c) > 127 for c in segment):
                        offending.append(f'{os.path.relpath(path, PACKAGE_DIR)}:{node.lineno}')
        self.assertEqual(offending, [], msg='print calls with non-ASCII characters')


if __name__ == '__main__':
    unittest.main()
