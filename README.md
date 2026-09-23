# JetVHeto 4.0.0

A Fortran code for the calculation of jet-veto cross sections and 
efficiencies in colour-singlet production, extending JetVHeto 3.0.0 
to N3LL resummation with NLO, NNLO or N3LO fixed-order matching.

## Example

Build with `make -j` after installing GNU Fortran, HOPPET, LHAPDF 6 and
CHAPLIN. Install the `NNPDF40_nnlo_as_01180` PDF set for this example. From the
source directory, run:

```sh
./jetvheto -in inputs/ptj1_xmuR_0.5_xmuF_0.5_xQ_0.5.fxd \
  -order 3 -fixed_order 3 -scheme anew -cross-section \
  -new_modlog -xM 0.5 -out matched.dat
```

This produces an N3LL+N3LO Higgs jet-veto prediction at 13.6 TeV, with
M=125 GeV and R=0.4, using the supplied HEFT fixed-order input. Output columns
are veto momentum [GeV], matched, resummed and fixed-order cross sections
[nb]. Input-card cross sections must also be converted into nb. For other predictions,
supply fixed-order inputs with the appropriate process, PDF, scales, radius
and jet algorithm; do not change these independently of the input calculation.

Running the command line

python3 scripts/reproduce_figures.py --input-dir inputs --output-dir outputs

creates a directory outputs/ containing the plots of Ref. XXXX.XXXX.

N3LL requires `-loop-mass none` (the default). At NNLL (`-order 2`), optional
small-R resummation is enabled by `-small-r -R0 1.0`; it is not available at
N3LL.

## References to cite

For public scientific work using this code, cite the original references 
relevant to the ingredients and options used.

### For N3LL predictions implemented starting from version 4.0.0

- S. Abreu, J. R. Gaunt, P. F. Monni, L. Rottoli, R. Szafron
  *Third-Order Logarithmic Resummation for Jet-Vetoed Higgs Production*:
  **final citation to be supplied by the authors**.
- S. Abreu, J. R. Gaunt, P. F. Monni, R. Szafron, 
  *The analytic two-loop soft function for leading-jet pT*,
  [arXiv:2204.02987](https://arxiv.org/abs/2204.02987).
- S. Abreu, J. R. Gaunt, P. F. Monni, L. Rottoli, R. Szafron, 
  *Quark and gluon two-loop beam functions for leading-jet pT and slicing at NNLO*,
  [arXiv:2207.07037](https://arxiv.org/abs/2207.07037).
- S. Abreu, J. R. Gaunt, P. F. Monni, L. Rottoli, R. Szafron, 
  *Three-loop rapidity anomalous dimension for jet-veto cross sections*:
  **citation and associated data-deposit DOI/version to be supplied by the authors**.


### For NNLL predictions, cite original work on the JetVHeto code

- A. Banfi, P. F. Monni, G. P. Salam and G. Zanderighi,
  *Higgs and Z-boson production with a jet veto*,
  [arXiv:1206.4998](https://arxiv.org/abs/1206.4998).
- A. Banfi, P. F. Monni and G. Zanderighi,
  *Quark masses in Higgs production with a jet veto*,
  [arXiv:1308.4634](https://arxiv.org/abs/1308.4634).
- A. Banfi et al., 
  *Jet-vetoed Higgs cross section in gluon fusion at N3LO+NNLL with small-R resummation*,
  [arXiv:1511.02886](https://arxiv.org/abs/1511.02886).


### Small-R resummation, when enabled

- M. Dasgupta, F. A. Dreyer, G. P. Salam and G. Soyez,
  *Small-radius jets to all orders in QCD*,
  [arXiv:1411.5182](https://arxiv.org/abs/1411.5182).


### Software, PDFs and fixed-order inputs

- **NNLOJET** A. Huss et al., 
  [arXiv:2503.22804](https://arxiv.org/abs/2503.22804).
- **HOPPET:** G. P. Salam and J. Rojo,
  [arXiv:0804.3755](https://arxiv.org/abs/0804.3755).
- **LHAPDF 6:** A. Buckley et al.,
  [arXiv:1412.7420](https://arxiv.org/abs/1412.7420).
- **Harmonic polylogarithms and HPL:** E. Remiddi and J. A. M. Vermaseren,
  [hep-ph/9905237](https://arxiv.org/abs/hep-ph/9905237); T. Gehrmann and
  E. Remiddi, [hep-ph/0107173](https://arxiv.org/abs/hep-ph/0107173).
- **CHAPLIN:** S. Buehler and C. Duhr,
  [arXiv:1106.5739](https://arxiv.org/abs/1106.5739).
