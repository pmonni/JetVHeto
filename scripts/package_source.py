"""Create a deterministic source-review archive; never upload or publish."""
from pathlib import Path
import hashlib,io,tarfile

ROOT=Path(__file__).resolve().parents[1]
DIST=ROOT/'dist';DIST.mkdir(exist_ok=True)
TOP='JetVHeto-n3ll-source-review'
files=[ROOT/n for n in ['README.md','README_SHORT.md','CHANGELOG.md','Makefile','COPYING','AUTHORS','.gitignore']]
for directory,suffixes in [('src',{'.f90','.f','.inc'}),('tests',{'.py','.f90','.fxd'}),
                          ('scripts',{'.py'}),('data',{'.txt','.wl'}),
                          ('docs',{'.md','.json'})]:
    files.extend(p for p in (ROOT/directory).rglob('*') if p.is_file() and p.suffix in suffixes)
# Examples may contain private prediction runs. Only ship the named public
# upstream fixture, when present; a test copy always lives under tests/fixtures.
example=ROOT/'examples/H125-LHC13-N3LO.fxd'
if example.is_file(): files.append(example)
files.extend(p for p in (ROOT/'licenses').iterdir() if p.is_file())
files=sorted(set(files))
for p in files:
    assert not p.is_symlink(),f'Symlink in source archive: {p}'
    if p.suffix in {'.f90','.f','.py','.md','.txt','.fxd'}:
        text=p.read_text()
        personal_roots=[chr(47)+name+chr(47) for name in ('Users','home')]
        assert not any(prefix in text for prefix in personal_roots), f'Personal absolute path: {p}'
inventory=''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT)}\n' for p in files)
(DIST/'SOURCE_SHA256.txt').write_text(inventory)
archive=DIST/(TOP+'.tar.gz')
# gzip mtime is fixed as well as tar metadata so repeated builds are identical.
import gzip
with archive.open('wb') as raw, gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as zipped:
    with tarfile.open(fileobj=zipped,mode='w') as tar:
        for name,data in [(str(p.relative_to(ROOT)),p.read_bytes()) for p in files]+[('SOURCE_SHA256.txt',inventory.encode())]:
            entry=tarfile.TarInfo(TOP+'/'+name)
            entry.size=len(data);entry.mode=0o644;entry.mtime=0
            tar.addfile(entry,io.BytesIO(data))
(DIST/'ARCHIVE_SHA256.txt').write_text(hashlib.sha256(archive.read_bytes()).hexdigest()+'  '+archive.name+'\n')
print(f'{archive.name}: {len(files)} source/data/documentation files; publication checklist remains in docs/PROVENANCE.md')
