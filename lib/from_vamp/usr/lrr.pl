#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames %bowtieFields));

my $dir = "/home/wilsonte/wga/lrr";

#-------------------------------------------------------------
#my @samples = ("F112R", "B3R", "B5R");
#my $arraySuffix = "_AFFY_MM9";

#my @samples = ("diploid1","diploid2","diploid3","diploid4","diploid5","diploid6",
#               "sperm01","sperm02","sperm03","sperm04","sperm05","sperm06","sperm07","sperm08","sperm09");
#my $arraySuffix = "_NIM_HG18";

#my @samples = ("7974", "7975");
#my $arraySuffix = "_ILL_HG18";

#my @samples = ("A3A2", "A1A1");
#my $arraySuffix = "_NIM_HG18";

#my @samples = ("diploid1","diploid2","diploid3","diploid4","diploid5","diploid6");
my @samples = ("sperm01","sperm02","sperm03","sperm04","sperm05","sperm06","sperm07","sperm08","sperm09");
my $arraySuffix = "_NIM_HG18";
#-------------------------------------------------------------

sub lrrPlot{
    foreach my $sample(@samples){
        my $binSize = 0.25;    
        my $table = "ARRAY_$sample$arraySuffix";  
          
        #my $tableSql = $table;
        my $tableSql ="SELECT arr.RATIO / decode(mn.MEDIAN_RATIO,0,0.001,mn.MEDIAN_RATIO) RATIO, arr.CHROMOSOME, arr.POSITION
                       FROM $table arr, SPERM_MEDIAN_RATIO mn
                       WHERE arr.CHROMOSOME = mn.CHROMOSOME
                         AND arr.POSITION = mn.POSITION";
        
        my $val = " round(log(2,decode(arr.RATIO,0,0.001,arr.RATIO))/$binSize)*$binSize ";
        my $sqlCommon = " FROM ($tableSql) arr WHERE arr.CHROMOSOME <= 19 AND arr.CHROMOSOME != 2 ";
        #my $sqlCommon = " FROM ($tableSql) arr WHERE arr.CHROMOSOME <= 22 ";
        my $valSql = "SELECT $val val $sqlCommon ";       
        my $countSql = "SELECT count(*) $sqlCommon ";
        my $sql = "SELECT val, count(*) / ($countSql) freq FROM ($valSql) GROUP BY val ORDER BY val";
        runSQL($sql, \my($value, $freq)); 
                 
        #my $file = "$dir/$table.lrr.csv";
        my $file = "$dir/$table.lrr.mn.csv";
        
        open my $fileH, ">", $file;
        print $fileH "LRR,FREQUENCY\n";
        while (fetchRow()){ print $fileH "$value,$freq\n" }
        close $fileH;
    }
}

1;



