"""Prevent default-real rounding in the N3LL reference coefficient kernels."""
from pathlib import Path
import re

source = Path(__file__).resolve().parents[1] / 'src/coefficient_functions_n3ll.f90'
# Consume strings/comments before checking complete numeric tokens.
tokens = re.compile(
    r"'(?:(?:'')|[^'])*'|\"(?:(?:\"\")|[^\"])*\"|![^\n]*"
    r'|(?<![\w.])(?P<real>(?:\d+\.\d*|\.\d+)(?:[eEdD][+-]?\d+)?'
    r'|\d+[eEdD][+-]?\d+)(?:_(?P<kind>\w+))?(?![\w.])'
)

def untyped(text):
    return [m.group() for m in tokens.finditer(text)
            if m.group('real') and not m.group('kind')
            and 'd' not in m.group('real').lower()]

assert untyped('7.324074074074074 + 1e-11') == ['7.324074074074074', '1e-11']
assert not untyped("1._dp + .5_dp + 1d-11 + 2.0D0 ! 8.574074074074074\n'1.23'")
bad = untyped(source.read_text())
assert not bad, f'Default-real literals in {source.name}: {bad}'
print('N3LL reference coefficient literal-precision check passed')
