#!/usr/local/bin/perl
use strict;
use warnings;
use Math::Trig;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

my ($offset, $delimiter)= (150, ",");
my %probes;

sub importAffyR{
    my (@inSamples) = @_;
    loadAffyProbes();
    foreach my $inSample(@inSamples){
        status("======================================\n$inSample\n");
        loadAffySampleR($inSample);
    }    
}

sub loadAffyProbes {
    status("loading Affy probes...\n");
    my ($fileH, $header) = getAffyHeader("$param{inputPath}/def/mdgProbes.csv", $delimiter);
    while (my $line = <$fileH>) {
        chomp $line;
        my @line = split($delimiter, $line);  
        my $id = $line[$$header{snpId}];
        $probes{$id}{chrom} = $line[$$header{CHROMOSOME}]; 
        $probes{$id}{pos} = $line[$$header{POSITION}]; 
        $probes{$id}{R} = $line[$$header{R}]; 
        $probes{$id}{T} = $line[$$header{T}]; 
    }    
    close $fileH;
}

sub loadAffySampleR {
    my ($inSample) = @_;
    my ($fileH, $header) = getAffyHeader("$param{inputPath}/cel/$inSample.CEL.AB.csv", $delimiter);
    my $arrayTable = newTable('Array', arrayTableName($inSample."R", 'Affymetrix', $param{refSeqBase}));
    my $arrayFile = "$arrayTable.csv";
    open my $outH, ">", $arrayFile or die "could not open $arrayFile"; 
    while (my $line = <$fileH>) {
        chomp $line;
        $line =~ s/\"//g;
        my ($id, $A, $B) = split($delimiter, $line); 
        $A = 2 ** $A - $offset;
        $B = 2 ** $B - $offset;
        my $R = sqrt($A**2 + $B **2);
        my $T = atan2($B, $A);
        my $ratio = $R / $probes{$id}{R};
        my $bFreq = calculateBAF($T, $probes{$id}{T});    
        my $zygosity = abs(0.5 - $bFreq) + 0.5; 
        print $outH join(",", $probes{$id}{chrom}, $probes{$id}{pos}, 
                              $ratio, $ratio, $bFreq, $zygosity, $zygosity, 
                              -1, 1, $A, $B, $T, $ratio, $bFreq)."\n";
#                             CHROMOSOME, POSITION, 
#                             RATIO, NORMALIZEDRATIO, BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, 
#                             CNVVALUE, GCSCORE, X, Y, THETA, R, BF
    }
    close $outH;  
    loadData($arrayFile, $arrayTable, ',', $fieldNames{Array});
}

1;

