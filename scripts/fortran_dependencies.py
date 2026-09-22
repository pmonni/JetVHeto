"""Generate local module dependencies; external modules belong to HOPPET."""
from pathlib import Path
import re, sys
sources=[Path(p) for p in sys.argv[1:]]
texts={p:p.read_text() for p in sources}
modules={}
for p,text in texts.items():
    for name in re.findall(r'^\s*module\s+(\w+)\s*$',text,re.M|re.I):
        modules[name.lower()]=p
for p,text in texts.items():
    used=re.findall(r'^\s*use\s+(?:[^:!\n]*::\s*)?(\w+)',text,re.M|re.I)
    deps=sorted({f'build/obj/{modules[n.lower()].stem}.o' for n in used
                 if n.lower() in modules and modules[n.lower()]!=p})
    print(f'build/obj/{p.stem}.o: '+ ' '.join(deps))
