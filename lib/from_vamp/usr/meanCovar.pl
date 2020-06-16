#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames %bowtieFields));

my $dir = "/home/wilsonte/wga/covariance";
#my @increments = (1, 2, 5, 10, 50, 100);
my @increments = (1,10,100,1000); 
my $binSize = 0.05;
my $bin = " round(lrr/$binSize)*$binSize ";

my @samples = (
    "F112R_AFFY_MM9",
    "B3R_AFFY_MM9",
    "B5R_AFFY_MM9",
    "diploid1_NIM_HG18",
    "diploid2_NIM_HG18",
    "diploid3_NIM_HG18",
    "diploid4_NIM_HG18",
    "diploid5_NIM_HG18",
    "diploid6_NIM_HG18",
    "sperm01_NIM_HG18",
    "sperm02_NIM_HG18",
    "sperm03_NIM_HG18",
    "sperm04_NIM_HG18",
    "sperm05_NIM_HG18",
    "sperm06_NIM_HG18",
    "sperm07_NIM_HG18",
    "sperm08_NIM_HG18",
    "sperm09_NIM_HG18",
    "7974_ILL_HG18",
    "7975_ILL_HG18",
    "A3A2_NIM_HG18", 
    "A1A1_NIM_HG18"
);

sub meanCoVar {
    my $file = "$dir/meancovar.csv";
    open my $fileH, ">", $file;
    print $fileH "SAMPLE";
    foreach my $increment(@increments){ print $fileH ",$increment" }
    print $fileH "\n";
    
    foreach my $sample(@samples){
        print $fileH "$sample";
        my $table = "ARRAY_$sample";  
        my $lrrTable = "$table\_mcv";
        runSQL("CREATE TABLE $lrrTable AS 
                        SELECT CHROMOSOME, POSITION, log(2,RATIO) lrr, 0 bin
                        FROM $table 
                        WHERE RATIO > 0
			  AND CHROMOSOME <= 22");  
        runSQL("UPDATE $lrrTable SET bin = $bin");   
        my ($mean, $stDev) = getStatsMCV($lrrTable);
        runSQL("UPDATE $lrrTable SET lrr = (lrr - $mean) / $stDev");
        runSQL("UPDATE $lrrTable SET bin = $bin");
        runSQL("SELECT count(*) FROM $lrrTable", \my($n));         
        fetchRow();        
        my $minN = $n / 1000;
        foreach my $increment(@increments){   
            my $sql = " SELECT sum(bin*adjMean)/count(*) cv
                        FROM (
                            SELECT bin, avg(adj) adjMean
                            FROM(
                                SELECT bin, lead(lrr,$increment,0) OVER (ORDER BY CHROMOSOME, POSITION) adj
                                FROM $lrrTable
                                ORDER BY CHROMOSOME, POSITION
                            )
                            GROUP BY bin
                            HAVING count(*) >= $minN
                        )";         
            runSQL($sql, \my($cv));
            fetchRow();  
            print $fileH ",$cv";      
        } 
        print $fileH "\n";
        dropTable($lrrTable);
    }  
    close $fileH; 
}

sub popCoVar {
    my $file = "$dir/popcovar.csv";
    open my $fileH, ">", $file;
    print $fileH "SAMPLE";
    foreach my $increment(@increments){ print $fileH ",$increment" }
    print $fileH "\n";

    foreach my $sample(@samples){
        print $fileH "$sample";
        my $table = "ARRAY_$sample";  
        my $lrrTable = "$table\_pcv";
        runSQL("CREATE TABLE $lrrTable AS 
                        SELECT CHROMOSOME, POSITION, log(2,RATIO) lrr, 0 bin
                        FROM $table 
                        WHERE RATIO > 0
                          AND CHROMOSOME <= 22"); 
        runSQL("UPDATE $lrrTable SET bin = $bin");  
        my ($mean, $stDev) = getStatsMCV($lrrTable);
        runSQL("UPDATE $lrrTable SET lrr = (lrr - $mean) / $stDev");
        foreach my $increment(@increments){   
            my $sql = " SELECT sum(idx*adj)/count(*) cv
                        FROM (
                            SELECT lrr idx, lead(lrr,$increment,0) OVER (ORDER BY CHROMOSOME, POSITION) adj
                            FROM $lrrTable
                            ORDER BY CHROMOSOME, POSITION
                        )";   
            runSQL($sql, \my($cv));
            fetchRow();  
            print $fileH ",$cv";      
        } 
        print $fileH "\n";
        dropTable($lrrTable);
    }  
}

sub getStatsMCV {
    my ($lrrTable) = @_;
    my $histSQL = " SELECT bin X, Count(bin) Y FROM ($lrrTable) WHERE bin != 0 GROUP BY bin ";
    runSQL("SELECT Max(Y) FROM ($histSQL)", \my($maxY));
    fetchRow();  
    defined $maxY or return undef;
    my $minY = 0.05 * $maxY;
    $histSQL = " SELECT X, Y FROM ($histSQL) WHERE Y >= $minY ";   
    runSQL("SELECT Min(X), Max(X) FROM ($histSQL)", \my($minX, $maxX));
    fetchRow();
    (defined $minX and defined $maxX) or return undef;
    my ($mean, $stDev);
    my $sGuess = ($maxX - $minX) / (2 * 2.45); #expectations of stdDev for perfect Gaussian given minY at 5% of maxY, mean = 0    
    my $sGuess3 = $sGuess * 3;
    runSQL("SELECT Avg(lrr), StdDev(lrr) 
            FROM ($lrrTable)
            WHERE bin <= $sGuess3 AND bin >= -$sGuess3", \($mean, $stDev));
    fetchRow();      
    return ($mean, $stDev);
}

1;
