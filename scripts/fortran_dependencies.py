"""Generate local module dependencies; external modules belong to HOPPET."""
from pathlib import Path
import re
import sys


def statements(text):
    """Handle free-form statements, not just the first USE on a physical line.

    Mask strings before removing comments/splitting semicolons, so literal
    text cannot create dependencies. Fortran doubled quote escapes are kept
    within the masked string. Join ampersand continuations after comments
    have been removed (including `& ! comment`).
    """
    text = re.sub(r"'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\"|![^\n]*", " ", text)
    text = re.sub(r"&[ \t]*\r?\n[ \t]*&?", " ", text)
    return re.split(r"[;\r\n]", text)


def dependency_lines(sources):
    parsed = {p: statements(p.read_text()) for p in sources}
    modules = {}
    for path, lines in parsed.items():
        for line in lines:
            match = re.fullmatch(r"\s*module\s+(\w+)\s*", line, re.I)
            if match:
                modules[match[1].lower()] = path
    for path, lines in parsed.items():
        deps = set()
        for line in lines:
            match = re.match(r"\s*use\b\s*(?:,\s*(intrinsic|non_intrinsic)\s*)?(?:::)?\s*(\w+)\b", line, re.I)
            if not match or (match[1] or '').lower() == 'intrinsic':
                continue
            provider = modules.get(match[2].lower())
            if provider is not None and provider != path:
                deps.add(f'build/obj/{provider.stem}.o')
        yield f'build/obj/{path.stem}.o: ' + ' '.join(sorted(deps))


if __name__ == '__main__':
    print('\n'.join(dependency_lines([Path(p) for p in sys.argv[1:]])))
