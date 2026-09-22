! User-facing help for the unified JetVHeto interface.
module banner
  implicit none
contains
  subroutine print_banner
    write(0,'(a)') 'JetVHeto - LL through N3LL jet-veto resummation'
    write(0,'(a)') 'Based on JetVHeto 3.0.0; see AUTHORS, COPYING and docs/PROVENANCE.md.'
  end subroutine
  subroutine print_help
    write(0,'(a)') 'Usage: jetvheto [options]'
    write(0,'(a)') '  -in FILE            Fixed-order input card (cross sections in nb)'
    write(0,'(a)') '  -out FILE           Output file (default /dev/stdout)'
    write(0,'(a)') '  -order 0|1|2|3      LL, NLL, NNLL, N3LL (default 2)'
    write(0,'(a)') '  -fixed_order 1|2|3  NLO, NNLO, N3LO matching (default 2)'
    write(0,'(a)') '  -scheme anew        Only supported matching scheme (default)'
    write(0,'(a)') '  -cross-section      Cross sections in nb instead of efficiencies'
    write(0,'(a)') '  -differential       Finite-difference spectrum instead of cumulant'
    write(0,'(a)') '  -new_modlog         Profile Q and switch off logarithms above xM*M'
    write(0,'(a)') '  -xM VALUE           Profile matching-scale ratio in [0.5,1] (default 0.5)'
    write(0,'(a)') '  -xQ VALUE           Initial resummation scale / M (default 0.5)'
    write(0,'(a)') '  -p VALUE            Original fixed-Q logarithm power (default 5)'
    write(0,'(a)') '  -truncate-rapidity-rge-at-n3ll  Strict radius exponent; order 3 only'
    write(0,'(a)') '                      Default: RadISH inclusive products beyond N3LL'
    write(0,'(a)') '  -jet-algorithm antikt|kt|ca   N3LL algorithm (default antikt)'
    write(0,'(a)') '  -fixed-order-algorithm NAME  Input algorithm (default antikt)'
    write(0,'(a)') '  -n3ll-data DIR      Boundary tables (default data/TMDs_ptj relative to cwd)'
    write(0,'(a)') '  -loop-mass none|t|t+b        Heavy-quark correction (default none)'
    write(0,'(a)') '  -resum-only         No matching; parameters supplied on the command line'
    write(0,'(a)') '  -proc H|DY -collider pp -sqrts ENERGY -M MASS -pdf_name SET'
    write(0,'(a)') '  -muR SCALE -muF SCALE -Q SCALE -R RADIUS -pdf_set MEMBER'
    write(0,'(a)') '  -nbins N -ptmin PT -ptmax PT  Resum-only grid; see README.md'
    write(0,'(a)') '  -h, -help           Show this help'
    write(0,'(a)') 'N3LL: unprimed, SU(3), nf=5, H in HEFT or DY, 0.05 <= R <= 1.'
    write(0,'(a)') 'Only ptj/anew is supported. -small-r [-R0 1] is available only at NNLL (antikt).'
    write(0,'(a)') '-small-r-ln2z and c1-squared options are rejected.'
    write(0,'(a)') 'At N3LL, -resum-only requires -cross-section; unused matched/FO columns are zero.'
    write(0,'(a)') 'No experimental acknowledgement flag is required. See docs/VALIDATION.md.'
  end subroutine
end module banner
