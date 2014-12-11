#!/usr/bin/perl -w
use strict;

my ($N) = @ARGV;
my $usage = "The default ALTSCAN uses a TOP 100 model. This script is to change the Top number as user wishes.
usage: $0 N
N is the number to change to.
";
die $usage if @ARGV < 1;

die("Please change directory to bin/ and run this script again!\n") unless (-e "../bin/change_top_N.pl");
system("chmod 644 ../src/zTrellis.h");
open(IN,"../src/zTrellis.h")||die("Cannot read ../zTrellis.h!\n");
my @tmp = <IN>;
close IN;
open(OUT,">../src/zTrellis.h")||die("Cannot write to ../zTrellis.h!\n");
foreach (@tmp){
    if(/^\#define top_num/){
	print OUT "\#define top_num $N\n";
    }
    else{
	print OUT $_;
    }
}

print ("Please do \:\n\tcd \.\.\n\tmake clean\n\tmake\n");

exit;
