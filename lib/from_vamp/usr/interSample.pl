#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames %bowtieFields));

my $dir = "/home/wilsonte/interSample";
#-------------------------------------------------------------
#my @samples = (["B5R","B6R"], ["B5R","B7R"], ["B6R","B7R"]);
#my $arraySuffix = "_AFFY_MM9";

#my @samples = (["diploid4","diploid5"],["sperm01","sperm02"],["sperm03","sperm04"],["sperm05","sperm06"],["sperm07","sperm08"]);
#my $arraySuffix = "_NIM_HG18";

#my @samples = (["A3A2","A1A1"]);
#my $arraySuffix = "_NIM_HG18";

my @samples = (["7974","7975"]);
my $arraySuffix = "_ILL_HG18";
#-------------------------------------------------------------
#my $binSize = 0.1;
#my @value = ("BAF", "Round(arr1.BFREQUENCY,1)", "Round(arr2.BFREQUENCY,1)", 0.1, 0.9);
#my $infTable = "MASK_CASTB6";
my $binSize = 0.25;
my @value = ("LRR", "round(log(2,decode(arr1.RATIO,0,0.001,arr1.RATIO))/$binSize)*$binSize", "round(log(2,decode(arr2.RATIO,0,0.001,arr2.RATIO))/$binSize)*$binSize", -1.5, 1.5);
my $infTable;
#-------------------------------------------------------------

sub interSample{
    foreach my $samples(@samples){
        my ($sample1, $sample2) = @$samples;
        my $table1 = "ARRAY_$sample1$arraySuffix";    
        my $table2 = "ARRAY_$sample2$arraySuffix"; 
        my ($values, $results, $counts, $csvFiles, $csvFileHs);
        $values = getBinnedValuesIS($table1, $table2);
        ($results, $counts) = getCovarianceIS($values);
        ($csvFiles, $csvFileHs) = createCsvFilesIS("$sample1\_$sample2");
        printCsvFilesIS($csvFileHs, $results, $counts);
        closeCsvFilesIS($csvFileHs);
    }
}

sub getBinnedValuesIS{
    my ($table1, $table2) = @_;
    my $sql;
#    if ($infTable){
#        $sql = "SELECT $value[1] val, $value[2] val
#                FROM $table1 arr1, $table2 arr2, $infTable inf
#                WHERE arr1.CHROMOSOME = inf.CHROMOSOME
#                  AND arr1.POSITION = inf.POSITION
#                  AND arr1.CHROMOSOME <= 19 AND arr1.CHROMOSOME != 2
#                  AND arr1.CHROMOSOME = arr2.CHROMOSOME
#                  AND arr1.POSITION = arr2.POSITION
#                ORDER BY arr1.CHROMOSOME, arr1.POSITION";
#    } else {
#        $sql = "SELECT $value[1] val, $value[2] val
#                FROM $table1 arr1, $table2 arr2
#                WHERE arr1.CHROMOSOME <= 19 AND arr1.CHROMOSOME != 2
#                  AND arr1.CHROMOSOME = arr2.CHROMOSOME
#                  AND arr1.POSITION = arr2.POSITION
#                ORDER BY arr1.CHROMOSOME, arr1.POSITION";
#    }  
    if ($infTable){
        $sql = "SELECT $value[1] val, $value[2] val
                FROM $table1 arr1, $table2 arr2, $infTable inf
                WHERE arr1.CHROMOSOME = inf.CHROMOSOME
                  AND arr1.POSITION = inf.POSITION
                  AND arr1.CHROMOSOME <= 22 
                  AND arr1.CHROMOSOME = arr2.CHROMOSOME
                  AND arr1.POSITION = arr2.POSITION
                ORDER BY arr1.CHROMOSOME, arr1.POSITION";
    } else {
        $sql = "SELECT $value[1] val, $value[2] val
                FROM $table1 arr1, $table2 arr2
                WHERE arr1.CHROMOSOME <= 22
                  AND arr1.CHROMOSOME = arr2.CHROMOSOME
                  AND arr1.POSITION = arr2.POSITION
                ORDER BY arr1.CHROMOSOME, arr1.POSITION";
    }          
    runSQL($sql, \my($value1,$value2));    
    my @values;  
    while (fetchRow()){push @values, [$value1,$value2]}
    return \@values;
}

sub getCovarianceIS{
    my ($values) = @_;
    my (%results, %counts);
    foreach my $i(0..(scalar(@$values)-1)){
        my $indexBin = getBinIS($$values[$i][0]);
        my $targetBin = getBinIS($$values[$i][1]);
        $results{$indexBin}{$targetBin}++;
        $counts{$indexBin}++;
    }
    return (\%results, \%counts);
}

sub getBinIS{
    my ($bin) = @_;
    return $bin + 0;
}

sub createCsvFilesIS{
    my ($table) = @_;
    my (@csvFiles, %csvFileHs);
    my $d = "$dir/$value[0]/$table";
    -d $d or mkdir $d;
    my $file = "$d/$table.csv";
    push @csvFiles, $file;
    open my $fileH, ">", $file;
    $csvFileHs{interSample} = $fileH;
    return (\@csvFiles, \%csvFileHs);
}

sub printCsvFilesIS{
    my ($fileHs, $results, $counts) = @_;
    my $fileH = $$fileHs{interSample};
    for (my $indexBin = $value[3]; $indexBin <= $value[4]; $indexBin += $binSize){
        print $fileH "$indexBin";
        for (my $targetBin = $value[3]; $targetBin <= $value[4]; $targetBin += $binSize){
            my $freq = 0;
            $$results{$indexBin}{$targetBin} and $freq = $$results{$indexBin}{$targetBin} / $$counts{$indexBin};
            print $fileH ",$freq";  
        }
        print $fileH "\n";
    }        
}

sub closeCsvFilesIS{
    my ($csvFileHs) = @_;;
    foreach my $increment(keys %$csvFileHs){
        my $fileH = $$csvFileHs{$increment};
        close $fileH;
    }
}

1;



