"""Regressions for serial/parallel module ordering, without a Fortran compiler."""
from pathlib import Path
import importlib.util
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('dependencies', ROOT/'scripts/fortran_dependencies.py')
dependencies = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dependencies)


class ModuleDependencies(unittest.TestCase):
    def test_statements_comments_strings_and_continuations(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            providers = []
            for name in ('first', 'second', 'third', 'fourth', 'fake', 'ieee_arithmetic'):
                path = root/(name+'.f90')
                path.write_text(f'module {name} ! definition comment\nend module\n')
                providers.append(path)
            consumer = root/'consumer.f90'
            consumer.write_text('''program consumer
use first; USE second; use :: third ! ; use fake
use, non_intrinsic :: & ! continued declaration
  & fourth, only: value
use, intrinsic :: ieee_arithmetic
print *, "text; use fake !", 'it''s a string; use fake'
end program
''')
            lines = list(dependencies.dependency_lines(providers+[consumer]))
            self.assertEqual(lines[-1], 'build/obj/consumer.o: build/obj/first.o build/obj/fourth.o build/obj/second.o build/obj/third.o')
            self.assertTrue(all(line.endswith(': ') for line in lines[:-1]))

    def test_actual_program_requires_all_semicolon_modules(self):
        sources = sorted((ROOT/'src').glob('*.f90'))+sorted((ROOT/'src').glob('*.f'))
        line = next(line for line in dependencies.dependency_lines(sources)
                    if line.startswith('build/obj/jetvheto.o:'))
        for name in ('expansion', 'resummation', 'reader', 'pdfs_tools', 'rad_tools'):
            self.assertIn(f'build/obj/{name}.o', line.split())


if __name__ == '__main__':
    unittest.main()
