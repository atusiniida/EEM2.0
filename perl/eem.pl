#!/usr/bin/perl

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin";
use Cwd;

my $base  = $FindBin::Bin;
$base =~ s/\/[^\/]*$//;

my $usage = "usage: $0 [-n numberOfChildren -m memorySize(Gbite)] expFile(*.tab) gsFile(*.gmt)\n" ;

my $mem = 12;  #memory (G bite)
my $n = 20;  # number of children

my $argv = join(" ", @ARGV);
if($argv =~ s/-n\s+([\d]+)//){
    $n = $1;
}
if($argv =~ s/-m\s+([\d]+)//){
    $mem = $1;
}
if($argv =~ /-\w\s/){
    die $usage;
}

@ARGV = split(" ", $argv);

my @exp = map {Cwd::abs_path($_)} grep {/\.tab$/ and -s} @ARGV or die $usage;
my @geneset = map {Cwd::abs_path($_)} grep {/\.gmt$/ and -s} @ARGV or die $usage;
my $eemBin = "$base/eemParallel_1.0.1/src/eem/eemParallel";
#my $mpirunBin = "/usr/local/package/openmpi/current_gcc_static/bin/mpirun";
#my $ldlbpath = "/usr/local/package/openmpi/current_gcc_static/lib:".$ENV{LD_LIBRARY_PATH};
my $mpirunBin = "/usr/local/package/openmpi/1.8.4-static-gcc/bin/mpirun";
my $ldlbpath = "/usr/local/package/openmpi/1.8.4-static-gcc/lib:".$ENV{LD_LIBRARY_PATH};
$ldlbpath = "/usr/local/package/java/jdk1.7.0_72_32/jre/lib:".$ENV{LD_LIBRARY_PATH};

my @script;

$n = (int($n/3)+1)*3;

#my  $qsuboption = "-pe mpi-fillup $n";
my  $qsuboption = "-pe mpi $n";

foreach my $expFile  (@exp){
    foreach my $gsFile (@geneset){
	$gsFile =~ /(.*?\/)*([^.]+)/;
        my $geneset = $2;
        $expFile =~ /(.*?\/)*([^.]+)/;
        my $exp = $2;
        my $id = $exp."_".$geneset;
	my $outFile = "${id}.eem";
	-s $outFile and next;
	my $scriptFile = "tmp${$}.${id}.sh";
	push(@script, $scriptFile);

open(OUT, ">$scriptFile");
print OUT<<"EOF";
#!/bin/tcsh
#\$ -S /bin/tcsh
#\$ -cwd $qsuboption
#\$ -v LD_LIBRARY_PATH=$ldlbpath

$mpirunBin   -np \$NSLOTS   -machinefile \$TMPDIR/machines   $eemBin  -R 0.05,0.1,0.15 -o $outFile  $expFile  $gsFile

EOF

     }
}

foreach(@script){
   #while(system("qsub  $_")){
    while(system("qsub -l s_vmem=${mem}G,mem_req=${mem}G,os6  $_")){   
     sleep(10)
     }
}

wait_for_SGE_finishing("tmp${$}.");
`rm tmp${$}.*`;


sub  wait_for_SGE_finishing{
    my $script = shift;
    my $cutoff;
    if(@_){
	$cutoff = shift;
    }else{
	$cutoff = 1;
    }
    $script = substr($script,0,10);
    while(1){
	while(system("qstat > /dev/null") != 0){
	    sleep(10);
	}
	my $out = `qstat| grep $script | wc`;
	$out =~ /\d+/;
	if($& < $cutoff ){
	    return;
	}else{
	    sleep(10);
	}
    }
}
