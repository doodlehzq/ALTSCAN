#!/usr/bin/perl -w
use strict;

my ($par, $seq)=@ARGV;
my $usage="This script calls altscan and then filter the results to remove the repeat ones.
usage: perl $0 <parameter> <input.fasta>

Now we only provide parameters for human.
Pay attention that \"input.fasta\" can include only one sequences.
";
die $usage if @ARGV < 2;
#check if executable altscan exists
my $dir = "";
if($0 =~/\//){
    my @parse = split("\/",$0);
    pop @parse;
    $dir = join("\/",@parse);
    $dir .= "\/";
}
else{
    $dir=".\/";
}
my $altscan = $dir."altscan";
unless(-e $altscan){
    die("An executable altscan doesn't exist! Make sure you install it correctly!\n");
}
my $filter=$dir."del_rep.pl";
unless(-e $altscan){
    die("del_rep.pl doesn't exist!\n");
}

#generate intermediate temporary file
my $inter_file_base = "temp";
my $inter_file = "temp0";
my $inter_n = 0;
while(-e $inter_file){
    $inter_n++;
    $inter_file = $inter_file_base.$inter_n;
}

#run altscan
print("$altscan $par $seq -top >$inter_file");
system("$altscan $par $seq -top >$inter_file");

#filter
system("perl $filter $inter_file");

#remove
system("rm -rf $inter_file");

exit;
