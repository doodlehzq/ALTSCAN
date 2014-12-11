#!/usr/bin/perl -w

use strict;
use Getopt::Std;

use vars qw($opt_s $opt_d);
getopts('sd:');

my $usage = "$0 [options] <input_file> 

Options:
  -d       del_information_output_file
";


die $usage if @ARGV<1;

my ($file,$output) = @ARGV;
my $i = -1;
my @tvec;
my $name;
my $score;
open(IN,$file)||die "can't open $file:$!\n";
if(defined($opt_d)){
    open(OUTD,">$opt_d")||die "open $opt_d error:$!\n";
}
while(<IN>){
    chomp;
    if(/\# Target Sequence: \>/){
	$name = $_;
    }
    if(/\# Score: ([+-\.\w]+)/){
	$score = $1;
    }
    if(/\#Transcript/){
	$i++;
	$tvec[$i] = init_tran(); 
    }
    if(/^[0-9a-z]+/){
	my @temp = split("\t",$_);
	$temp[5] = $score;
	$temp[8]=~/(gene_id \"[\w\.]+\"); (transcript_id \"[\w\.]+\")/;
	put_tran($tvec[$i],new_sf($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$1,$2));
    }
}
close IN;

#check identical

my $j;
my $k;
for($j=0;$j<=$i;$j++){
    if($tvec[$j] != 0){
	for($k=$j+1;$k<=$i;$k++){
	    if($tvec[$k] != 0){
		if(check_tr($tvec[$j], $tvec[$k])){
		    $tvec[$k] = 0;
		}
	    }
	}
    }
}
if(defined($opt_d)){
print OUTD "$name\n";
}

my $num = Output(\@tvec, $name, $output);
if(defined($opt_d)){
print OUTD "\ndelete number: $num\n";
}


 


sub new_sf{             
    my ($chr,$src,$type,$start,$end,$score,$sym,$phase,$gid,$tid)=@_;
    my %hash = (
	"chr" => $chr,
	"source" => $src,
	"type" => $type,
	"score" => $score,
	"gid" => $gid,
	"tid" => $tid,
	"start" => $start,
	"end" => $end,
	"sym" => $sym,
	"phase" => $phase,
    );
    return \%hash;
}

sub init_tran{
    my @a;
    return \@a;
}
sub put_tran{       #tran array address,sf         
    my($tran,$sf)=@_;
    my $l = @{$tran};
    ${$tran}[$l] = $sf;    
}

sub check_tr{          #transcript 1, transcript 2;
                              # identical flag=1
    my ($tr1, $tr2)=@_;
    my $flag = 1;
    if(scalar(@{$tr1})==scalar(@{$tr2})){
	my $i;
	for($i=0;$i<@{$tr1};$i++){
	    if(!check_sf(${$tr1}[$i],${$tr2}[$i])){
		$flag = 0;
		last;
	    }
	}
    }
    else{
	$flag = 0;
    }
    return $flag;
}


sub check_sf{                     #identical flag=1;
    my ($sf1, $sf2)=@_;
    my $flag = 1;
    if(${$sf1}{start} != ${$sf2}{start} 
       || ${$sf1}{end} != ${$sf2}{end} ){
	$flag =0;
    }
    return $flag;
}


sub Output{              #\@tvec,$name ,output file
    my ($ad, $title, $outfile) = @_;
    print "$name\n";
    my $num = 0;
    my $d1;
    for($d1=1;$d1<=@{$ad};$d1++){ 
	if(${$ad}[$d1-1] != 0){
	    $num++;
	    if(defined($opt_d)){
		print OUTD "$d1\t";
	    }
	    print "\#Transcript $num: \n";
	    my $d2;
	    for($d2=0;$d2<@{${$ad}[$d1-1]};$d2++){

		print  "${${${$ad}[$d1-1]}[$d2]}{chr}\t";
		print  "${${${$ad}[$d1-1]}[$d2]}{source}\t";
		print  "${${${$ad}[$d1-1]}[$d2]}{type}\t";
		print  "${${${$ad}[$d1-1]}[$d2]}{start}\t";
		print  "${${${$ad}[$d1-1]}[$d2]}{end}\t";
		print  "${${${$ad}[$d1-1]}[$d2]}{score}\t";
		print  "${${${$ad}[$d1-1]}[$d2]}{sym}\t";
		print  "${${${$ad}[$d1-1]}[$d2]}{phase}\t";
		print  "${${${$ad}[$d1-1]}[$d2]}{gid}; ";
		print  "${${${$ad}[$d1-1]}[$d2]}{tid};\n";
	    }
	}
    }
    my $del = @{$ad} - $num;
    return $del;
}
