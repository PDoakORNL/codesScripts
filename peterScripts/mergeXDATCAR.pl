use strict;
use warnings;
use Getopt::Long;
use PDL;

our $firstX = 0;
our $lastX = 0;
GetOptions("firstX:i", => \$firstX,
	   "lastX:i", => \$lastX);

our $configCount = 1;
our $header = "";
open my $xdmf, '>', "XDATCAR_merged";

print "Merging XDATCAR: ";

if($firstX != 0) {
    for(my $i = $firstX; $i<=$lastX; $i++){
	my $sX = sprintf("%02i", $i);
	print "XDATCAR_${sX} ";
	read_write($xdmf, "XDATCAR_${sX}");
    }
}


my @slist = sort {($a =~ /([\d\.]+)\/XDATCAR/)[0] <=> ($b =~ /([\d\.]+)\/XDATCAR/)[0]} @ARGV;

foreach my $xdatFile (@slist) {
    print "$xdatFile ";
    read_write($xdmf, $xdatFile);
}
print "\n";

close $xdmf;

# just to avoid code duplication
sub read_write {
    my $xdmFile = shift(@_);
    my $xdatFile = shift(@_);
    open my $xdat, '<', $xdatFile or die "failed to open $xdatFile $| \n"; 
    my $tempHeader = $header;
    for(my $j = 0; $j < 7; $j++){
	$_ = <$xdat>;
	if($header eq "" ) { # eat header for later $xdatCAR
	    print $xdmFile $_;
	    $tempHeader .= $_;
	}
    }
    $header = $tempHeader;

    while(<$xdat>) {
	if($_=~/(Direct|Cartesian)(.*\s*=\s+)(\d+)/) {
	    print $xdmFile "${1}${2}${configCount}\n";
	    $configCount++;
	    next;
	}
	print $xdmFile $_;
    }
    close $xdatFile;
}

