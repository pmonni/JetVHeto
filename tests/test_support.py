"""Self-contained CLI fixture, not a physical prediction with a changed PDF."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parents[1]
PDF=os.environ.get('PDF_SET','NNPDF40_nnlo_as_01180')
BUILD=ROOT/'build/test-results'
BUILD.mkdir(parents=True,exist_ok=True)
CARD=BUILD/'H125-LHC13-test.fxd'
source=ROOT/'tests/fixtures/H125-LHC13-N3LO.fxd'
# Only the test fixture uses a caller-selected PDF. The published example is
# preserved unchanged: its fixed-order coefficients belong to PDF4LHC15.
lines=source.read_text().splitlines()
CARD.write_text('# Test-only PDF substitution; not for physics predictions\n'+'\n'.join(
    'pdf_name = '+PDF if line.startswith('pdf_name =') else line for line in lines)+'\n')
NROWS=sum(1 for line in lines if line.strip() and not line.startswith('#') and '=' not in line)
