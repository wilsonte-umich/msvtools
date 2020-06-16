#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{xx} or $command{xx} = ['multiThread', '24:00:00', 5000, 0];
defined $command{compareIntrons} or $command{compareIntrons} = ['singleThread', '4:00:00', 2000, 0];

sub compareIntrons {
    my ($sample1, $sample2) = @_;
    my ($crosstabTable, $crosstabFile) = crosstabIntrons($sample1, $sample2);
}

sub crosstabIntrons {
    my ($sample1, $sample2) = @_;
    status("comparing hit densities per intron for samples $sample1 and $sample2...\n");
    my $iMapTable1 = getTableName('IMap', "$sample1");
    my $iMapTable2 = getTableName('IMap', "$sample2");
    my $file = "$param{inputPath}/compareIntrons_$sample1\_$sample2.csv";
    open my $fileH, ">", $file or die "could not open $file\n";   
    my $rankSQL = getTwoSampleCrosstab($sample1, $sample2, $iMapTable1, $iMapTable2,
                                       "NAME2, CHROMOSOME, START_, CORREND_", "NORMALIZEDDENSITY");
    my $crosstabTable = "imap_$sample1\_$sample2\_ct";
    dropTable($crosstabTable);
    runSQL("CREATE TABLE $crosstabTable AS $rankSQL");
    print $fileH "RANK,NAME2,CHROMOSOME,START_,CORREND_,$sample1,$sample2,$sample1\/$sample2,$sample2\/$sample1\n";
    runSQL($rankSQL);
    my $genes = fetchAllHashRef();
    foreach my $gene(@$genes){
        print $fileH "$$gene{RANK},$$gene{NAME2},$$gene{CHROMOSOME},$$gene{START_},$$gene{CORREND_},";
        my ($s1, $s2, $s12, $s21) = 
        ($$gene{"\U$sample1"},           $$gene{"\U$sample2"},
         $$gene{"\U$sample1\_$sample2"}, $$gene{"\U$sample2\_$sample1"});
        defined $s1 or $s1 = "";
        defined $s2 or $s2 = "";
        defined $s12 or $s12 = "";
        defined $s21 or $s21 = "";
        print $fileH "$s1,$s2,$s12,$s21\n";
    }
    close $fileH;
    return ($crosstabTable, $file);
}


1;


