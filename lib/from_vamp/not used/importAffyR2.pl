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
    
#    print join(" ", keys %$header);
#    print "\n";

foreach my $key(keys %$header){print"$key = $$header{$key}\n"}

    
    my $arrayTable = newTable('Array', arrayTableName($inSample."R", 'Affymetrix', $param{refSeqBase}));
    my $arrayFile = "$arrayTable.csv";
    open my $outH, ">", $arrayFile or die "could not open $arrayFile"; 
    while (my $line = <$fileH>) {
        
        print "\nline = $line\n";
        
        chomp $line;
        $line =~ s/\"//g;
        my @line = split($delimiter, $line);
        
foreach my $key(keys %$header){print"$key = $line[$$header{$key}]\n"}
        
#    print join("\n value =", @line);
#    print "\n";

    
    
        my $ratio = $line[$$header{RR}];
        my $bFreq = $line[$$header{BAF}];
        my $zygosity = abs(0.5 - $bFreq) + 0.5; 
        
    exit;
              
        print $outH join(",", $line[$$header{CHROMOSOME}], $line[$$header{POSITION}], 
                              $ratio, $ratio, $bFreq, $zygosity, $zygosity, 
                              -1, 1, $line[$$header{A}], $line[$$header{B}], $line[$$header{T}], $ratio, $bFreq)."\n";
##                             CHROMOSOME, POSITION, 
##                             RATIO, NORMALIZEDRATIO, BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, 
##                             CNVVALUE, GCSCORE, X, Y, THETA, R, BF


exit;

    }
    close $outH;  
    loadData($arrayFile, $arrayTable, ',', $fieldNames{Array});
}

1;
