#!/usr/bin/perl -w
# Script that runs a check that the default results for DY and Higgs 
# created by running ./jetvheto are the same as the benchmarks ones 
#
# authors: A. Banfi, P.F. Monni, G.P. Salam and G. Zanderighi
# copyright: GNU GPLv3

use strict; 
use warnings; 

system('rm -f check.diff');

for my $proc('DY','H125'){  
    my $fileold="benchmarks/$proc-LHC8-R05-xmur050-xmuf050-xQ050-NNLO+NNLLa.res";
    my $filenew="benchmarks/$proc-LHC8-R05-xmur050-xmuf050-xQ050-NNLO+NNLLa.res.tmp";
    system("./jetvheto -in fixed-order/$proc-LHC8-R05-xmur050-xmuf050-log.fxd -xQ 0.50 -out $filenew -scheme a -order 2 -ptmax 500.01"); 
    if (-e "$filenew"){
	print "Created file $filenew\n";
    } 
    else{
	print "Check FAILED: failed to create $filenew\n"; 
	exit 1;
    }
    
    # extract the non-commented parts of the reference (old) and newly
    # generated files; essentially, for each file, take all the
    # (uncommented) numbers and push them onto an array.
    open(FILEold, $fileold) || die("Could not open file $fileold!");
    my @dataold; 
    while(<FILEold>) {
	if (substr($_,0,1) ne "#"){
	    my @tmp = split(' ',$_); 
	    push(@dataold, @tmp); 
	}
    }
    close(FILEold); 
    open(FILEnew, $filenew) || die("Could not open file $filenew!");
    my @datanew; 
    while(<FILEnew>) {
	if (substr($_,0,1) ne "#"){
	    my @tmp = split(' ',$_); 
	    push(@datanew, @tmp); 
	}
    }
    close(FILEnew); 

    # make sure the arrays of numbers from the two files have the same length
    my $oldlength = scalar(@dataold); 
    my $newlength = scalar(@datanew); 
    if ($oldlength != $newlength){
	system("diff $fileold $filenew > check.diff");  
	print "File lengths does not match: $oldlength vs $newlength \n"; 
	print "Look at differences between $fileold and $filenew in check.diff\n";
	print "Check FAILED\n"; 
	exit 1; 
    }
    
    # then compare the numbers from the two files to ensure they're close enough
    for (my $i = 0; $i < $oldlength; $i++) {
	my $diff = $dataold[$i] - $datanew[$i]; 
        # be quite loose in our checks, because we have found that
        # LHAPDF gives slightly different alphas values depending on
        # the machine that we run on (at the level of 10^(-10); near
        # the pole in alphas, these differences get magnified.
	if (abs($diff) > 1E-7 
            || (abs($dataold[$i]) > 1E-10 && abs($diff) >= 1E-5*abs($dataold[$i]))) {
	    system("diff $fileold $filenew > check.diff");  
	    print "Files $fileold and $filenew do not match:\n    $dataold[$i] v. $datanew[$i]\n";
            print "For more info look at differences in check.diff\n"; 
	    print "Check FAILED\n"; 
	    exit 1; 
	}
    }
# Remove temporary files created 
    system ("rm $filenew"); 
    print "Removed temporary file $filenew\n";
}

print "Check PASSED\n"; 

