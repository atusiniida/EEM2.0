#!/usr/local/bin/perl 

use strict;
use warnings;

my %out;
my %p;

while(<>){
    my @tmp = split("\t");
    $out{$tmp[0]} = $_;
    $p{$tmp[0]} = $tmp[1];
}

foreach(sort {$p{$b} <=> $p{$a} } keys %out){
    print $out{$_};
}
