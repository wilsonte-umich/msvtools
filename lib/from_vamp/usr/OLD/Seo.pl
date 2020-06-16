#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %fieldNames %refSeqs));

sub loadSeo{
    my $inputFile = "$param{inputPath}/Seo.csv";
    open my $inputH, "<", $inputFile;
    my $outTable = newTable('CNVs', 'Seo_hg18');
    my $outFile = "$outTable.csv";
    open my $outputH, ">", $outFile;
    while (<$inputH>){
        chomp $_;
        my ($index, $sample, $chr, $length, $start, $end, $log2R) = split(",", $_);
        $end or next;
        my $cnvType = 1;
        $log2R and $log2R > 0 and $cnvType = 8;
        print $outputH "$index\::$sample,$cnvType,$chr,$start,$end\n";
    }
    close $inputH;
    close $outputH;
    loadData($outFile, $outTable, ",", $fieldNames{CNVs});
}

1;


