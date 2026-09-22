#!/usr/bin/env python
# script to attempt to gather all the uncertainties and put them into an uncertainty band
# according to different prescriptions
from __future__ import division # always do floating point division

usage="""
Usage: python/gather-uncertainties.py [-h] [... options ...]

This script reads files carried out with different run configurations
(e.g. scale choices, schemes, etc.) and assembles the information to
produce resummed and/or fixed order central values and uncertainty
bands, both for the jet-veto efficiency and for the jet-veto cross
section.

NB: this script is a research tool and not 100% robust. Do not rely on
it blindly, but instead check the output carefully.

Options
-------

Files are assumed to be named according to a convention

  DIR/HEAD-xmurXXX-xmufXXX-xQXXX-ORDERSCHEME[.R0-YYY].res

where XXX indicates numerical values (100 == 1.0*mH, 050=0.5*mH,
etc.), and DIR, HEAD and ORDER can be set through explicit options.


  -dir DIR         directory in which individual results files are stored
                   [defaults to benchmarks-n3lo/]. The name must include 
                   the trailing slash. 
  
  -filehead HEAD   common part (head) of filename across all runs
                   [defaults to H125-LHC13-R04-heft]

                   (Alternatively, set DIR to the combination of
                   directory and filehead, and set filehead to "")
  
  -order ORDER     the order of the calculation that is being used
                   e.g. N3LO+NNLL [default] or NNLO+NNLL. Note that
                   even if you want fixed-order results

  -type {matched|fixed|resummed}

                   indicates which column should be taken from the
                   output files. Note that to get fixed-order results,
                   you'll want to specify, e.g. 

                       "-order N3LO+NNLL -type fixed"

Various options control the scale, scheme and R0 options that are used:

  -schemes SCHEMEDEF 

                   a SCHEMEDEF of "a-RF-Q-R0,b" (the default) means
                   that the result is the envelope of the following:

                   - renormalisation and factorisation variations
                     ("RF"), Q variations("Q") and R0 variations
                     ("R0") in scheme a,
                   - the central scale and R0 choice from scheme b.

                   The schemes that are potentially available are
                   a,b,c,cunx (called c' in 1511.xxxxx) and d.

                   If is also possible to use schemes such as b:a, as
                   decribsed in appendix A of 1511.xxxxx

  -central {mh2|mh} 
                   indicates whether the central scale should be mH/2
                   [default] or mH. 

  -indepscales {0|1} 

                   If this is 1 [default], then muR and muF are varied
                   independently with the restriction 1/2 <= muR/muF
                   <= 2. If this is 0, then muR and muF are kept equal.

                   Scale variation is always performed by a factor of
                   two around the central scale.

  -Qdef  QQQ       the central Q choice in units of 100*Q/mH, usually 
                   3 digits [default: 050, i.e. 0.5*mH].

  -Qvariations XXX,YYY
                   The variations of Q0 that should be carried out,
                   in the format XXX,YYY [default: 033,075]. This means
                   that one examines Q=XXX*mH/100 and Q=YYY*mH/100.

  -R0def  R0DEF    The central value for R0 [1.0]. If it's set to 0 or "", 
                   then the string ".R0-YYY" is simply left off the filename.

  -R0variations  R1,R2
                   The variations of R0 [default: "0.5,2.0"]

Other options control the kind of answer one wants to obtain:

  -answer {efficiency|efficiencySigma|directSigma|STSigma}
                   Determines the quantity that is printed out. 

                   - "efficiency" is the jet veto efficiency

                   - "efficiencySigma" is the jet-vetoed cross
                     section, using the JVE method, i.e. the
                     hypothesis that uncertainties on the total cross
                     section are uncorrelated with those on the
                     efficiency.

                   - "directSigma" means the jet-vetoed cross section
                     with direct scale variation (i.e. for each scale
                     choice one calculates the jet-vetoed cross
                     section and then one takes the envelope of all
                     cross sections). 

                   - "STSigma" means the Stewart-Tackmann prescription
                     for the jet-vetoed cross section. I.e. the
                     uncertainties are taken to be uncorrelated
                     between the total cross section and the 1-jet
                     cross section, with the jet-vetoed cross section
                     obtained as the difference between them. This
                     only makes sense for fixed-order results ("-type
                     fixed" above).

                   All three "Sigma" options give the same central
                   values, but differ in their uncertainty bands.

  -symmetric {0|1} If it's "1" then uncertainties are symmetrised
                   (taken the maximum of lower and upper uncertainties
                   to set the range. [Default: "0"]

  -diff {0|1}      If it's "1" then the output is a differential efficiency
                   or cross section [default: "0"]. The way this works
                   is that the derivative is taken separately for each
                   input file. Then uncertainties are worked out
                   as the envelope of the differential efficiency or cross
                   section (depending on "answer" above).

Examples
--------

To obtain the HEFT jet-veto efficiency results of 1511.XXXXX at
N3LO+NNLL+LLR, run the following:

  python/gather-uncertainties.py -dir benchmarks-n3lo/ -filehead H125-LHC13-R04-heft -order N3LO+NNLL -schemes a-RF-Q-R0,b

The output indicates the values of the otions, the list of input
qfiles that was used, the total cross section and its uncertainty and
the jet veto efficiency v. pt. The units are those in the input files
(by default, nb).

To obtain the jet-veto cross section results of 1511.XXXXX at
N3LO+NNLL+LLR with full mass effects, in the JVE method (table 2):

  python/gather-uncertainties.py -dir benchmarks-n3lo/ -filehead H125-LHC13-R04-mb-475-mt-1725 -order N3LO+NNLL -schemes a-RF-Q-R0,b -answer efficiencySigma 

To instead get 

  - the fixed-order results, add "-type fixed"

  - the NNLO+NNLL results, replace "-order N3LO+NNLL" 
    with "-order NNLO+NNLL -R0def 0 -schemes a-RF-Q,b"

To reproduce the pt=25, R=0.4, 8 TeV results in the table on p.4 of
1206.4998, run

  python/gather-uncertainties.py -dir benchmarks/ -filehead H125-LHC8-R04  -order NNLO+NNLL -R0def 0 -Qvariations 025,100 -schemes a-RF-Q,b,c  -answer efficiency

(replace "efficiency" with "efficiencySigma" to get the cross
sections; this won't quite agree with the results in that table,
because 1206.4998 used a slightly different total cross section from
that in the input files here).


"""

from hfile import *
import cmdline
import numpy as np
import pprint
import sys
#from datetime import datetime
#datetime.utcnow().isoformat(' ')

#--- first decide if we need usage instructions --------
if (cmdline.present("-h")):
    print usage
    sys.exit(0)

#--- then set up default options, with a little bit of pre-processing
opt = {}
opt['dir'] = 'benchmarks-n3lo/'
opt['filehead'] = 'H125-LHC13-R04-heft'
opt['order'] = 'N3LO+NNLL'
opt['symmetric'] = False
opt['answer'] = 'efficiency' # or efficiencySigma, directSigma, STSigma
# a-RF-Q-R0, b means use scheme "a" with RF(ren,fact) scale
# variations, Q scale variations and R0 scale variations. Supplement it
# with scheme b, central value
# The dashes are optional, the commas are not
opt['schemes'] = 'a-RF-Q-R0,b'  
opt['Qdef'] = '050'
opt['Qvariations'] = '033,075'
opt['R0def'] = '1.0' # use "" or "0" if there's not an R0 string on the filename
opt['R0variations'] = '0.5,2.0'
opt['type'] = 'matched' # or 'fixed'
opt['indepscales'] = 1
opt['central'] = 'mh2' # central renormalisation and factorisation scale
opt['default_scale']='050'
opt['diff'] = False
opt['rescaling'] = False # establishes whether the NNLO and N3LO corrections are rescaled by the inclusive Born ratio factor for mass effects

orderN = 0 # the upper index of the order (3 for N3LO)

def main():

    deftype = bool

    # determine the options based on the command line
    for key in opt.iterkeys():
        isbool = isinstance(opt[key],(bool))
        opt[key] = cmdline.value("-"+key, opt[key])

    # make sure we are safe...
    cmdline.assert_all_options_used()

    # change default scale if the central scale is mh
    if (opt['central'] == 'mh'): opt['default_scale'] = '100'
    
    # write out the options
    print "# "+cmdline.cmdline()
    print "# Settings are: "
    for op in sorted(opt.iterkeys()):
        print "#  "+op+" = ", opt[op]

    # the files have sigmas up to order N3LO even when the order we're
    # dealing with is NNLO; so need to manually determine an order to use
    global orderN # needed to make sure we alter the global variable
    if   (re.match(r'^N3LO',opt['order'])): orderN = 3
    elif (re.match(r'^NNLO',opt['order'])): orderN = 2
    elif (re.match(r'^NLO', opt['order'])): orderN = 1
    elif (re.match(r'^LO',  opt['order'])): orderN = 0
    else: raise NameError("Could not interpret order option: "+opt['order'])
    
    # collect the configurations we'll be wanting;
    # start off by creating the array and defining some shorthand names
    configs = []
    xmuR='xmuR'
    xmuF='xmuF'
    scheme='scheme'
    R0='R0'

    # first put the central scales, scheme a
    #configs.append({})
        
        
    # now work through the requested schemes
    for sch in opt['schemes'].split(','):
        # scheme core name is just whatever is lowercase or contains a colon: so
        # remove everything that is not lowercase or a colon
        schname = re.sub(r'[^a-z:].*','',sch)
        #print 'schname is ', schname
        configs.append({scheme: schname})
        if (re.search('RF', sch)):
        
            ## PFM: add ren. & fact. scales with central 100-100
            ## first symmetric scales
            if opt['central'] == 'mh':
                configs.append({scheme:schname, xmuR:'200', xmuF:'200'})
                configs.append({scheme:schname, xmuR:'050', xmuF:'050'})
                # then asymmetric scales if we want them
                if (opt['indepscales']):
                    configs.append({scheme:schname, xmuR:'050', xmuF:'100'})
                    configs.append({scheme:schname, xmuR:'200', xmuF:'100'})
                    configs.append({scheme:schname, xmuR:'100', xmuF:'050'})
                    configs.append({scheme:schname, xmuR:'100', xmuF:'200'})

            # add ren. & fact. scales
            # first symmetric scales
            if opt['central'] == 'mh2':
                configs.append({scheme:schname, xmuR:'100', xmuF:'100'})
                configs.append({scheme:schname, xmuR:'025', xmuF:'025'})
                # then asymmetric scales if we want them
                if (opt['indepscales']):
                    configs.append({scheme:schname, xmuR:'025', xmuF:'050'})
                    configs.append({scheme:schname, xmuR:'100', xmuF:'050'})
                    configs.append({scheme:schname, xmuR:'050', xmuF:'025'})
                    configs.append({scheme:schname, xmuR:'050', xmuF:'100'})

        if (re.search('Q', sch)):
            # then add Q scales
            for xQ in opt['Qvariations'].split(','):
                configs.append({scheme:schname, 'xQ': xQ})
        if (re.search('R0', sch)):
            # then add R0 scales
            for R0 in opt['R0variations'].split(','):
                configs.append({scheme:schname, 'R0': R0})
            
    # # option of having all scales for scheme B, or just a single scale
    # if (opt['allScalesB']):
    #     nconfigs = len(configs)
    #     for i in range(0,nconfigs):
    #         newconfig = configs[i].copy()
    #         newconfig[scheme] = 'b'
    #         configs.append(newconfig)
    # else:
    #     configs.append({'scheme':'b'})

    # now see what we have
    pp = pprint.PrettyPrinter()
    files = []
    fileset = set() # to help avoid duplicate files (not really crucial)
    print "# Using the following files: "
    for config in configs:
        file = File(**config)
        if (file.filename not in fileset):
            files.append(file)
            fileset.add(file.filename)
            print '# ', file.filename

            
    # decide what column to use
    if (opt['type'] == 'fixed'):
        column = 3 
    elif (opt['type'] == 'matched'):
        column = 1
    elif (opt['type'] == 'resummed'):
        column = 2
    else:
        raise NameError("Did not recognize 'type' "+opt['type']);
    print "# using {0} column (col {1} in fortran numbering) to read efficiency".format(opt['type'],column+1)


    # now choose exactly what kind of answer we want
    if   (opt['answer'] == 'directSigma'):
        result = SigmaEnvelope(files, column=column);
    elif (opt['answer'] == 'direct1JetSigma'):
        result = Sigma1JetEnvelope(files, column=column);
    elif (opt['answer'] == 'efficiency'):
        result = EfficiencyEnvelope(files, column=column);
    elif (opt['answer'] == 'efficiencySigma'):
        result = EfficiencySigma(files, column=column);
    elif (opt['answer'] == 'efficiency1JetSigma'):
        result = EfficiencySigma(files, column=column, do1Jet=True);
    elif (opt['answer'] == 'STSigma'):
        result = StewartTackmann(files, column=column);
    else: raise NameError("Could not interpret 'answer' = ", opt['answer'])

    sigma_err = sigma_and_error(files)
    print '# sigma_tot and uncertainty: {0} {1} +{2}'.format(sigma_err[0],min(sigma_err[1],sigma_err[2]),max(sigma_err[1],sigma_err[2]))
    print "# columns for "+opt['answer']+": pt, central, min, max"
    print reformat(files[0].pt_values(), result[0], result[1], result[2])
        
#----------------------------------------------------------------------
def SigmaEnvelope(files, column):
    '''Return the envelope of the cross sections using correlated cross
       section and efficiency choices'''

    Sigmas = []
    for file in files:
        Sigmas.append(file.efficiencies[:,column]*file.sigma)

    return envelope(Sigmas)

#----------------------------------------------------------------------
def Sigma1JetEnvelope(files, column):
    '''Return the envelope of the 1-jet inclusive (integrated) cross sections using correlated cross
       section and efficiency choices'''

    Sigmas = []
    for file in files:
        Sigmas.append((1-file.efficiencies[:,column])*file.sigma)

    return envelope(Sigmas)
        
#----------------------------------------------------------------------
def EfficiencyEnvelope(files, column, do1Jet = False):
    '''Return the envelope of the efficiencies'''

    effs = []
    for file in files:
        if (do1Jet): effs.append(1 - file.efficiencies[:,column])
        else       : effs.append(file.efficiencies[:,column])
            
    return envelope(effs)

#----------------------------------------------------------------------
def EfficiencySigma(files, column, do1Jet=False):
    '''Return Sigma and its uncertainty based on the efficiency method

    If do1Jet is True then instead return the result for the cross
    section of having at least 1 jet above the given pt.
    '''

    effs = EfficiencyEnvelope(files, column, do1Jet)
    sigma_err = sigma_and_error(files)

    # first convert the central efficiency into a central sigma
    Sigmas = []
    Sigmas.append(effs[0] * sigma_err[0])

    # then deal with the lower and upper errors
    # we add the fractional efficiency error and fractional sigma error in quadrature,
    # and then apply the factor to the central value
    relerr_down = np.sqrt((sigma_err[1]/sigma_err[0])**2 + ((effs[1]-effs[0])/effs[0])**2)
    relerr_up   = np.sqrt((sigma_err[2]/sigma_err[0])**2 + ((effs[2]-effs[0])/effs[0])**2)
        
    Sigmas.append(Sigmas[0]*(1-relerr_down))
    Sigmas.append(Sigmas[0]*(1+relerr_up))

    return Sigmas

#----------------------------------------------------------------------
def StewartTackmann(files, column):
    'work in progress'
    SigmaBars = []
    sigma_err = sigma_and_error(files)
    for file in files:
        SigmaBars.append((1 - file.efficiencies[:,column])*file.sigma)

    # then get the envelope
    SigmaBarsEnv = envelope(SigmaBars)
    
    Sigmas = []
    # first get the central value
    Sigmas.append(sigma_err[0] - SigmaBarsEnv[0])

    # then get the lower and upper edges
    # watch out: Sigmabar upper edge adds to sigma lower edge to give Sigma lower edge
    Sigmas.append(Sigmas[0] - np.sqrt(sigma_err[1]**2+(SigmaBarsEnv[2]-SigmaBarsEnv[0])**2))
    Sigmas.append(Sigmas[0] + np.sqrt(sigma_err[2]**2+(SigmaBarsEnv[1]-SigmaBarsEnv[0])**2))

    return Sigmas

#----------------------------------------------------------------------
def envelope(arrays):
   '''Return the envelope of the arrays; if symmetric = True, then assume
   that arrays[0] is the central value and return an envelope that has
   been symmetrised around it.
   '''

   # first get the central result
   result = [ arrays[0], None, None ]

   # then add the extrema
   extrema = minmax(arrays)
   result[1] = extrema[0]
   result[2] = extrema[1]

   # then handle the situation where things are supposed to be symmetric
   if (opt['symmetric']):
       symerr = maximum([abs(result[1]-result[0]),abs(result[2]-result[0])])
       result[1] = result[0] - symerr
       result[2] = result[0] + symerr
       
   return result

#----------------------------------------------------------------------
def sigma_and_error(files):

    central = files[0].sigma
    sigmas = map(lambda f : f.sigma, files)
    error_lo = min(sigmas) - central
    error_hi = max(sigmas) - central

    if (opt['symmetric']):
        error_hi = max(abs(error_lo), abs(error_hi))
        error_lo = - error_hi

    return [central, error_lo, error_hi]

#----------------------------------------------------------------------
def filename(xmuR=None, xmuF=None, xQ=None, scheme='a', R0=None):
#def filename(xmuR = '050', xmuF='050', xQ='050', scheme='a', R0=None):
#def filename(xmuR = '100', xmuF='100', xQ='050', scheme='a', R0=None):

    match = re.search(r'([a-z]+):([a-z]+)',scheme)
    if (match):
        print "# Found composite difference scheme", match.group(0), \
              'consisting of schemes', match.group(1), match.group(2)
        return (filename(xmuR, xmuF, xQ, match.group(1), R0),
                filename(xmuR, xmuF, xQ, match.group(2), R0),
                filename(scheme=match.group(2)))

    # PM: set the central scales
    if (xmuR == None):
        xmuR=opt['default_scale']
    if (xmuF == None):
        xmuF=opt['default_scale']
    if (xQ   == None):
        xQ  =opt['Qdef']
        
    file = opt['dir']+opt['filehead']+'-xmur'+xmuR+'-xmuf'+xmuF+'-xQ'+xQ+'-'+opt['order']+scheme
    if (R0 == None): R0 = opt['R0def']
    if (R0 != '' and R0 != '0'): file += '.R0-'+R0
        
    if (opt['rescaling']):#PM: modify filename if rescaling is required
        file += '.rescaled.res'
    else:
        file += '.res'
    return file


#----------------------------------------------------------------------
class File:
    '''Class to store the information from the matched file specified
    through the **config elements (cf. filename()) and the global
    options in opt[...]
    '''

    def __init__(self, **config):
        # We allow File to be created either with an explicit filename
        # or with a more elaborate set of configuration options, which are
        # then passed to the filename(...) routine
        if ('file' in config):
            self.filename = config['file']
        else:
            self.filename = filename(**config)

        # handle special case where filename is a tuple, which for now
        # (2015-08-20) indicates that one should process it as a
        # special difference scheme. (This is for things like schemes
        # a:b, which use differences between those two schemes)
        if (isinstance(self.filename, tuple)):
            print "discovered filename is a tuple"
            # first, build subsiduary File objects.
            # We expect the triplet for b:a to be composed of:
            #   b(scale, etc.), a(scale, etc.), a(central)
            self.triplet = []
            for file in self.filename:
                self.triplet.append(File(file=file))
            # then process them to get the information we need
            # get the cross section from the "a(central)"
            self.sigmas = self.triplet[2].sigmas
            self.sigma  = self.triplet[2].sigma
            # get the efficiency from an appropriate difference
            self.efficiencies = self.triplet[0].efficiencies-self.triplet[1].efficiencies\
                                  + self.triplet[2].efficiencies
            # then exit
            return
            
        self.efficiencies = get_array(self.filename,'')
        #print opt['diff'], self.efficiencies.shape
        # get derivative to get a differential distribution
        if (opt['diff']):
            #pass
            npt  = self.efficiencies.shape[0]
            ncol = self.efficiencies.shape[1]
            deriv = np.empty((npt-1,self.efficiencies.shape[1]))
            deriv[:,0] = 0.5*(self.efficiencies[1:npt,0] + self.efficiencies[0:npt-1,0])
            for icol in range(1, ncol):
                deriv[:,icol] = (self.efficiencies[1:npt,icol] 
                                 - self.efficiencies[0:npt-1,icol])/(
                    self.efficiencies[1:npt,0] - self.efficiencies[0:npt-1,0])
            self.efficiencies = deriv

        # now reread the file to get the cross sections
        # (there must be a better way...)
        file = open(self.filename, 'r')
        self.sigma = None
        while True:
            line = file.readline();
            if (not line): break
            if (re.match(r'^# total cross sections',line)):
                newline = re.sub(r'^.*= *','',line)
                self.sigmas = np.array(map(float, newline.split()))
                self.sigma = self.sigmas[0:orderN+1].sum()
                #print self.sigma
                break
        if (self.sigma == None):
            raise hfile.Error("could not find sigmas in file "+self.filename)

    def pt_values(self):
        'Returns a numpy array of the pt values'
        return self.efficiencies[:,0]

    def Sigma(self, column):
        'Returns a numpy array for Sigma (in the requested column)'
        return self.sigma * self.efficiencies[:,column]
    
if __name__ == '__main__':
    main()
