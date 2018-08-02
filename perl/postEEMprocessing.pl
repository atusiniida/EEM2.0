#!/usr/local/bin/perl
use strict;
use warnings;

$0 =~ /.*(\/){0,1}perl/ and unshift(@INC, $&);
require  EEM;
set_environment();


my $home = get_home_dir();
my $usage = "usage: $0  [-p pvalueCutoff -c correlationCutoff -m minCulusterSize] eemFile  expFile\n";

my $p = 6;
my $c = 0.70;
my $m = 0;

my $tmp = join(" ",@ARGV);
if($tmp =~ s/-p (\S+)//){
    $p = $1;
}
if($tmp =~ s/-c (\S+)//){
    $c = $1;
}
if($tmp =~ s/-m (\S+)//){
    $m = $1;
}
my @tmp = split(" ", $tmp);
@tmp < 2 and die $usage;
my $eemFile = shift(@tmp);
my $expFile = shift(@tmp);

open(OUT, ">tmp${$}.R");
print OUT "eemFile <- '$eemFile'\n";
print OUT "expFile <- '$expFile'\n";
print OUT "Pcut <- $p\n";
print OUT "treeCut <- $c\n";
print OUT "minClusterSize  <- $m\n";
print OUT "source('$home/R/processEEMresult.R')\n";
close(OUT);
my $status = system("R --no-save < tmp${$}.R >& /dev/null");
if($status==0){
    `rm tmp${$}.R`;
}
exit($status);
