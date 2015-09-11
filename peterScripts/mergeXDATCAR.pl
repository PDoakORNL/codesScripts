#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $firstX = "1";
my $lastX = "1";
GetOptions("firstX:i", => \$firstX,
	   "lastX:i", => \$lastX);

my $configCount = 1;
open XDM, ">XDATCAR_merged";

for(my $i = $firstX; $i<=$lastX; $i++){
    my $sX = sprintf("%02i", $i);
    open XDAT, "<XDATCAR_${sX}" or die "failed to open XDATCAR_${sX} $| \n"; 
    if($i != $firstX) { # eat header for later XDATCAR
	for(my $j = 0; $j < 7; $j++){
	    my $line = <XDAT>;
	    print $line;
	}
    }
    while(<XDAT>) {
	if($_=~/(Direct configuration=\s+)(\d+)/) {
	    print XDM "${1}${configCount}\n";
	    $configCount++;
	    next;
	}
	print XDM $_;
    }
}
