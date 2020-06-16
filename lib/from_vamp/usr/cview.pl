#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames %bowtieFields));

my $dir = "/home/wilsonte/wga/cview";

#-------------------------------------------------------------
my @samples = ("diploid1","diploid2","diploid3","diploid4","diploid5","diploid6",
               "sperm01","sperm02","sperm03","sperm04","sperm05","sperm06","sperm07","sperm08","sperm09");
my $arraySuffix = "_NIM_HG18";
#-------------------------------------------------------------

sub cview{
    foreach my $sample(@samples){
        my $table = "ARRAY_$sample$arraySuffix";  
	my $mnTable = "ARRAY_$sample"."mn$arraySuffix";  
        my $sql = "SELECT arr.CHROMOSOME, arr.POSITION, log(2,nullif(arr.RATIO,0)) LRR, log(2,nullif(mn.RATIO,0)) LRR_MN
			FROM $table arr, $mnTable mn
			WHERE arr.CHROMOSOME = mn.CHROMOSOME
  			  AND arr.POSITION = mn.POSITION
			ORDER BY arr.CHROMOSOME, arr.POSITION";
        runSQL($sql, \my($chrom,$pos,$lrr,$mnLrr)); 
        my $file = "$dir/$table.cview.csv";      
        open my $fileH, ">", $file;
        print $fileH "CHROMOSOME,POSITION,LRR,LRR_MN\n";
        while (fetchRow()){ print $fileH "$chrom,$pos,$lrr,$mnLrr\n" }
        close $fileH;
    }
}

1;



