#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames %bowtieFields));

my $dir = "/home/wilsonte/wga/covariance";
my $scriptFile = "$dir/covariance.R";
#my @increments = (1, 2, 5, 10, 50, 100);
my @increments = (1,10); 
#-------------------------------------------------------------
#my @samples = ("F112R", "B3R", "B5R");
#my $arraySuffix = "_AFFY_MM9";

my @samples = ("diploid1","sperm09");
my $arraySuffix = "_NIM_HG18";

#my @samples = ("diploid1","diploid2","diploid3","diploid4","diploid5","diploid6",
#               "sperm01","sperm02","sperm03","sperm04","sperm05","sperm06","sperm07","sperm08","sperm09");
#my $arraySuffix = "_NIM_HG18";

#my @samples = ("7974", "7975");
#my $arraySuffix = "_ILL_HG18";

#my @samples = ("A3A2", "A1A1");
#my $arraySuffix = "_NIM_HG18";

#-------------------------------------------------------------
#my $binSize = 0.1;
#my @value = ("BAF", "Round(arr.BFREQUENCY,1)", 0.1, 0.9);
#my $infTable = "MASK_CASTB6";
my $binSize = 0.25;
my @value = ("LRR", "round(log(2,RATIO)/$binSize)*$binSize", -1.5, 1.5);
my $infTable;
#-------------------------------------------------------------

sub covariance {
    foreach my $sample(@samples){
        my $table = "ARRAY_$sample$arraySuffix";    
        my ($values, $results, $counts, $means, $csvFiles, $csvFileHs);
        $values = getBinnedValuesCV($table);
        ($results, $counts, $means) = getCovarianceCV($values);
        printMeansCV ($table, $counts, $means);
        # ($csvFiles, $csvFileHs) = createCsvFilesCV($table);
        # printCsvFilesCV($csvFileHs, $results, $counts);
        # closeCsvFilesCV($csvFileHs);
    }
}

sub getBinnedValuesCV{
    my ($table) = @_;
    my $sql;
    ##mouse
    # if ($infTable){
        # $sql = "SELECT $value[1] val
                # FROM $table arr, $infTable inf
                # WHERE RATIO > 0
                # AND arr.CHROMOSOME = inf.CHROMOSOME
                # AND arr.POSITION = inf.POSITION
                # AND arr.CHROMOSOME <= 19 AND arr.CHROMOSOME != 2
                # ORDER BY arr.CHROMOSOME, arr.POSITION";
    # } else {
        # $sql = "SELECT $value[1] val
                # FROM $table arr
                # WHERE RATIO > 0
                # AND arr.CHROMOSOME <= 19 AND arr.CHROMOSOME != 2
                # ORDER BY arr.CHROMOSOME, arr.POSITION";
    # }  
    ##human
    if ($infTable){
        $sql = "SELECT $value[1] val
                FROM $table arr, $infTable inf
                WHERE RATIO > 0
                  AND arr.CHROMOSOME = inf.CHROMOSOME
                  AND arr.POSITION = inf.POSITION
                  AND arr.CHROMOSOME <= 22 
                ORDER BY arr.CHROMOSOME, arr.POSITION";
    } else {
        $sql = "SELECT $value[1] val
                FROM $table arr
                WHERE RATIO > 0
                  AND arr.CHROMOSOME <= 22
                ORDER BY arr.CHROMOSOME, arr.POSITION";
    }          
    runSQL($sql, \my$value);    
    my @values;  
    while (fetchRow()){push @values, $value}
    return \@values;
}

sub getCovarianceCV{
    my ($values) = @_;
    my (%results, %counts, %means);
    foreach my $index(0..(scalar(@$values)-101)){
        my $indexBin = getBinCV($$values[$index]);
        foreach my $increment(@increments){
            my $targetValue = $$values[$index + $increment];
            my $targetBin = getBinCV($targetValue);
            $results{$increment}{$indexBin}{$targetBin}++;
            $counts{$increment}{$indexBin}++;
            $means{$increment}{$indexBin} += $targetValue;
        }   
    }
    return (\%results, \%counts, \%means);
}

sub getBinCV{
    my ($bin) = @_;
    return $bin + 0; #forces bin to be handled as a number
}

sub printMeansCV {
    my ($table, $counts, $means) = @_;
    my $d = getOutDirCV($table);
    my $file = "$d/$table.means.csv";
    open my $fileH, ">", $file;
    print $fileH "INCREMENT";
    for (my $indexBin = $value[2]; $indexBin <= $value[3]; $indexBin += $binSize){ print $fileH ",$indexBin" }
    print $fileH "\n";
    foreach my $increment(sort {$a <=> $b} keys %$counts) {
        print $fileH "$increment";
        for (my $indexBin = $value[2]; $indexBin <= $value[3]; $indexBin += $binSize){ 
            my $mean = sprintf("%.3f", $$means{$increment}{$indexBin} / $$counts{$increment}{$indexBin});
            print $fileH ",$mean";
        }
        print $fileH "\n";
    } 
    close $fileH;
}

sub createCsvFilesCV{
    my ($table) = @_;
    my (@csvFiles, %csvFileHs);
    foreach my $increment(@increments){
        my $d = getOutDirCV($table);
        my $file = "$d/$table.$increment.csv";
        push @csvFiles, $file;
        open my $fileH, ">", $file;
        $csvFileHs{$increment} = $fileH;
    }
    return (\@csvFiles, \%csvFileHs);
}

sub getOutDirCV {
    my ($table) = @_;
    my $d = "$dir/$value[0]/$table";
    -d $d or mkdir $d;   
    return $d;
}

sub printCsvFilesCV{
    my ($fileHs, $results, $counts) = @_;
    foreach my $increment(@increments){
        my $fileH = $$fileHs{$increment};
        for (my $indexBin = $value[2]; $indexBin <= $value[3]; $indexBin += $binSize){
            print $fileH "$indexBin";
            for (my $targetBin = $value[2]; $targetBin <= $value[3]; $targetBin += $binSize){
                my $freq = 0;
                $$results{$increment}{$indexBin}{$targetBin} 
                    and $freq = $$results{$increment}{$indexBin}{$targetBin} / $$counts{$increment}{$indexBin};
                print $fileH ",$freq";  
            }
            print $fileH "\n";
        }        
    }
}

sub closeCsvFilesCV{
    my ($csvFileHs) = @_;;
    foreach my $increment(@increments){
        my $fileH = $$csvFileHs{$increment};
        close $fileH;
    }
}

# sub popcovar {
    # 
    # my @samples = (
# "F112R_AFFY_MM9",
# "B3R_AFFY_MM9",
# "B5R_AFFY_MM9",
# "diploid1_NIM_HG18",
# "diploid2_NIM_HG18",
# # "diploid3_NIM_HG18",
# # "diploid4_NIM_HG18",
# # "diploid5_NIM_HG18",
# # "diploid6_NIM_HG18",
# # "sperm01_NIM_HG18",
# # "sperm02_NIM_HG18",
# "sperm03_NIM_HG18",
# "sperm04_NIM_HG18",
# # "sperm05_NIM_HG18",
# # "sperm06_NIM_HG18",
# "sperm07_NIM_HG18",
# "sperm08_NIM_HG18",
# "sperm09_NIM_HG18",
# "7974_ILL_HG18",
# # "7975_ILL_HG18",
# # "A3A2_NIM_HG18", 
# "A1A1_NIM_HG18"
    # );
# 
    # my $file = "$dir/pearsons.csv";
    # open my $fileH, ">", $file;
    # print $fileH "SAMPLE";
    # foreach my $increment(@increments){ print $fileH ",$increment" }
    # print $fileH "\n";
# 
    # foreach my $sample(@samples){
        # print $fileH "$sample";
        # my $table = "ARRAY_$sample";  
        # my $lrrTable = "$table\_lrr";
        # runSQL("CREATE TABLE $lrrTable AS 
                        # SELECT CHROMOSOME, POSITION, log(2,RATIO) lrr
                        # FROM $table 
                        # WHERE RATIO > 0");     
        # runSQL("SELECT Avg(lrr), StdDev(lrr), count(*) FROM $lrrTable", \my($mean,$sd,$n));
        # fetchRow();
        # foreach my $increment(@increments){   
            # #my $sql = " SELECT CORR(idx,tar)
            # my $sql = " SELECT sum(((idx - $mean)/$sd) *  ((tar - $mean)/$sd))/$n cv
                        # FROM (
                            # SELECT lrr idx, lead(lrr,$increment,0) OVER (ORDER BY CHROMOSOME, POSITION) tar
                            # FROM $lrrTable
                            # ORDER BY CHROMOSOME, POSITION 
                        # )";
            # runSQL($sql, \my($cv));
            # fetchRow();  
            # print $fileH ",$cv";      
        # }  
        # print $fileH "\n";
        # dropTable($lrrTable);
    # }  
    # close $fileH; 
# }

1;
