#!/usr/local/bin/perl 

use strict;
use warnings;

my $usage = "usage: $0  A_B.eem C.gmt D.gmt....\n";

my $eem = shift or die $usage;
my $gmt = join(" ", @ARGV) or die $usage;
chomp(my @tmp = `cut -f 1 $gmt`);
my %seen;
map { $seen{$_}=1 } @tmp;
open(IN, $eem);
while(<IN>){
    my @tmp = split("\t");
    if($seen{$tmp[0]}){
	print;
    }
}
