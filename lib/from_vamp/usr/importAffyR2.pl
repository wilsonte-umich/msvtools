#!/usr/local/bin/perl
use strict;
use warnings;
use Math::Trig;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

my $delimiter = ",";

sub importAffyR2{
    my (@inSamples) = @_;
    foreach my $inSample(@inSamples){
        status("======================================\n$inSample\n");
        loadAffySampleR2($inSample);
    }    
}

sub loadAffySampleR2 {
    my ($inSample) = @_;
    my ($fileH, $header) = getAffyHeader("$param{inputPath}/cel/$inSample.CEL.ABRT.csv", $delimiter);
    my $arrayTable = newTable('Array', arrayTableName($inSample."R", 'Affymetrix', $param{refSeqBase}));
    my $arrayFile = "$arrayTable.csv";
    open my $outH, ">", $arrayFile or die "could not open $arrayFile"; 
    while (my $line = <$fileH>) {
        chomp $line;
        $line =~ s/\"//g;
        my @line = split($delimiter, $line);
        my $ratio = $line[$$header{RR}];
        $ratio <= 0 and $ratio = 0.001;
        my $bFreq = $line[$$header{BAF}];
        $bFreq <= 0 and $bFreq = 0.001;
        my $zygosity = abs(0.5 - $bFreq) + 0.5;      
        print $outH join(",", $line[$$header{CHROMOSOME}], $line[$$header{POSITION}], 
                              $ratio, $ratio, $bFreq, $zygosity, $zygosity, 
                              -1, 1, $line[$$header{A}], $line[$$header{B}], $line[$$header{T}], $ratio, $bFreq)."\n";
##                             CHROMOSOME, POSITION, 
##                             RATIO, NORMALIZEDRATIO, BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, 
##                             CNVVALUE, GCSCORE, X, Y, THETA, R, BF
    }
    close $outH;  
    loadData($arrayFile, $arrayTable, ',', $fieldNames{Array});
}

1;
