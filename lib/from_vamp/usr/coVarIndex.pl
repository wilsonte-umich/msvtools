#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames %bowtieFields));

my $dir = "/home/wilsonte/wga/covariance";
my $scriptFile = "$dir/covariance.R";
my @increments = (1);
#-------------------------------------------------------------
#my $chromFilter = "CHROMOSOME <= 19 AND CHROMOSOME != 2";

#my @samples = ("F112R", "B3R", "B5R");
#my $arraySuffix = "_AFFY_MM9";


my $chromFilter = "CHROMOSOME <= 22";

#my @samples = ("sperm09");
my @samples = ("diploid1","diploid2","diploid3","diploid4","diploid5","diploid6",
               "sperm01","sperm02","sperm03","sperm04","sperm05","sperm06","sperm07","sperm08","sperm09");
my $arraySuffix = "_NIM_HG18";

#my @samples = ("7974", "7975");
#my $arraySuffix = "_ILL_HG18";

#my @samples = ("A3A2", "A1A1");
#my $arraySuffix = "_NIM_HG18";

#-------------------------------------------------------------

sub coVarIndex{
    my $file = "$dir/coVarIndex.csv";
    open my $fileH, ">", $file or die "could not open $file\n";
    print $fileH "SAMPLE,INCREMENT,SLOPE,INTERCEPT,R2\n";
    foreach my $sample(@samples){
        my $table = "ARRAY_$sample$arraySuffix";    
	foreach my $increment(@increments){
	    my $sql = "SELECT round(REGR_SLOPE(tar,idx),3) slope, round(REGR_INTERCEPT(tar,idx),3) intercept, round(REGR_R2(tar,idx),3) r2
			FROM (  SELECT val idx, lead(val,$increment,0) OVER (ORDER BY CHROMOSOME, POSITION) tar
				FROM (  SELECT CHROMOSOME, POSITION, log(2,RATIO) val
					FROM $table 
					WHERE RATIO > 0 
					  AND $chromFilter) 
			        ORDER BY CHROMOSOME, POSITION)";
	    runSQL($sql,\my($slope,$intercept,$r2));
	    fetchRow();
	    print $fileH join(",", ($sample, $increment, $slope, $intercept, $r2))."\n";
	}
    }
    close $fileH;
}


1;



