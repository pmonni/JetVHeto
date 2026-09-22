"""Remove only known, generated products inside this source tree."""
from pathlib import Path
import shutil
root=Path(__file__).resolve().parents[1]
build=root/'build'
if build.is_dir() and not build.is_symlink():shutil.rmtree(build)
exe=root/'jetvheto'
if exe.is_file():exe.unlink()
for source in (root/'tests').glob('test_*.f90'):
    binary=source.with_suffix('')
    if binary.is_file():binary.unlink()
