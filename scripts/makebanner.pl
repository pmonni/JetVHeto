#!/usr/bin/perl -w
#
# Script that generates automatically that extracts header and help
# information from a readme file and generates a f90 file that prints
# this information on the screen

$usage = 'makebanner.pl <README file> <f90 file>';

@ARGV || die $usage."\n";

$readme = shift @ARGV;

$outfile = shift @ARGV;

print "Automatically generating $outfile from $readme \n";

open(IN,"<$readme");
open(OUT,">$outfile");


print OUT "! This file is generated automatically by \n";
print OUT "! makebanner.pl $readme $outfile";

print OUT '
module banner 

contains 

  subroutine print_banner

';

while($line=<IN>) {
    chomp $line;
    if ($line =~ /The program is written in Fortran 95./) {last;}
    print OUT "    write(0,*) \'$line\'\n";
}

print OUT '
  end subroutine print_banner
';


print OUT '
  subroutine print_help

';

# move the file handle to the beginning of the help section of the README
while($line=<IN>) {
    if ($line =~ /Usage:/) {
	chomp $line;
	print OUT "    write (0,*) \'$line\'\n";
	last;
    }
}

# copy the remaining lines of the README file
while($line=<IN>) {
	chomp $line;
	print OUT "    write (0,*) \'$line\'\n";
}

print OUT '
  end subroutine print_help
';

print OUT '
end module banner
';
close IN;
close OUT;
