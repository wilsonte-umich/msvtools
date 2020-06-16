#!/usr/local/bin/perl
use strict;
use warnings;
use Math::Trig;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

my $theta100 = atan2(1, 0); #theta value for 100% B allele frequency, i.e. A = 0;
my $theta50 = atan2(1, 1); #theta value for 50% B allele frequency, i.e. A = B
my (%refs, %pos);

#$param{inputPath} = "/home/shared/data";
$param{affyRef} = "B3_SNP10_198";
$param{affyDefPath} = "/home/wilsonte/mdg/def";
$param{affyArray} = "MOUSEDIVm520650.na31";
#$param{refSeqBase} = "mm9";
#importAffy("B3_SNP10_198", "B4_SNP10_199", "B5_SNP10_200", "B6_SNP10_201", "B7_SNP10_202", "F1cross12");

sub importAffy{
    my (@inSamples) = @_;
    loadAffyRefs();
    loadAffyPositions();
    foreach my $inSample(@inSamples){
        print "======================================\n$inSample\n";
        loadAffySample($inSample);
    }    
}

sub loadAffyRefs{
    print "loading Affy reference values...\n";
    my $delimiter = "\t";
    my ($fileH, $header) = getAffyHeader("$param{inputPath}/$param{affyRef}/R_theta_refs.txt", $delimiter);
    while (my $line = <$fileH>) {
        chomp $line;
        my @line = split($delimiter, $line);  
        my $id = $line[$$header{probeset_id}];
        foreach my $ref (keys %$header){
            $ref eq "probeset_id" and next;
            $refs{$id}{$ref} = $line[$$header{$ref}]; 
        }  
    }    
    close $fileH;
}

sub loadAffyPositions{
    print "loading Affy chromosome positions...\n";
    my $delimiter = ",";
    my ($fileH, $header) = getAffyHeader("$param{affyDefPath}/$param{affyArray}.annot.csv", $delimiter);
    my $idHeader = "\"Probe Set ID\"";
    my $chromHeader = "\"Chromosome\"";
    my $posHeader = "\"Physical Position\"";  
    while (my $line = <$fileH>) {
        chomp $line;
        my @line = split($delimiter, $line);  
        my $id = purgeQuotes($line[$$header{$idHeader}]);
        $pos{$id}{chrom} = purgeQuotes($line[$$header{$chromHeader}]);
        $pos{$id}{pos} = purgeQuotes($line[$$header{$posHeader}]); 
    }    
    close $fileH;
}

sub purgeQuotes{
    my ($value) = @_;
    $value =~ m/\"(.*)\"/ or return $value;
    return $1;
}

sub loadAffySample{
    my ($inSample) = @_;
    print "  normalizing and printing sample values\n";
    my ($fileH, $header) = getAffyHeader("$param{inputPath}/$inSample/quant-norm.pm-only.brlmm-p.summary.txt", "\t");
    my $arrayTable = newTable('Array', arrayTableName($inSample, 'Affymetrix', $param{refSeqBase}));
    my $arrayFile = "$arrayTable.csv";
    open my $outH, ">", $arrayFile or die "could not open $arrayFile"; 
    while (my $line = <$fileH>) {
        my ($id, $alleleA, $aVals) = parseSummaryLine($line, $header);
        $alleleA eq "A" or die "allele A not where expected";
        $line = <$fileH>;
        my ($idB, $alleleB, $bVals) = parseSummaryLine($line, $header);
        $id eq $idB or die "failed id pairing in summary file";
        foreach my $sample (keys %{$header}){
        
            ($sample eq $inSample or $sample eq "$inSample.CEL") or next;
            
            my $pos = $pos{$id}{pos};                       
            my $chr = $pos{$id}{chrom};      
            my $chrom = $refSeqs{$param{refSeqBase}}{"chr$chr"};            
            my $a = $$aVals[$$header{$sample}];
            my $b = $$bVals[$$header{$sample}];
            my ($R, $theta) = transform($a, $b);
            my $ratio = $R / $refs{$id}{R_ref};
            my $bFreq = calculateBAF($theta, $refs{$id}{theta_ref});
            my $zygosity = abs(0.5 - $bFreq) + 0.5; 
            print $outH join(",", $chrom, $pos, 
                                  $ratio, $ratio, $bFreq, $zygosity, $zygosity, 
                                  -1, 1, $a, $b, $theta, $ratio, $bFreq)."\n";
    #                             CHROMOSOME, POSITION, 
    #                             RATIO, NORMALIZEDRATIO, BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, 
    #                             CNVVALUE, GCSCORE, X, Y, THETA, R, BF
        }
    }
    close $outH;  
    loadData($arrayFile, $arrayTable, ',', $fieldNames{Array});
}

sub calculateBAF{
    my ($theta, $theta_ref) = @_;
    if ($theta > $theta_ref) {
        return (($theta - $theta_ref) / ($theta100 - $theta_ref) * 0.5) + 0.5;
    } elsif ($theta < $theta_ref) {
        return $theta / $theta_ref * 0.5;
    } else {
        return 0.5;
    }
}

1;

