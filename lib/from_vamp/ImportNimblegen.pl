#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#microarray analysis, import Nimblegen aCGH into VAMP
#########################################################################

use vars(qw(%param %types %fields %fieldNames %refSeqs));  #common parameters used by most analysis scripts
my ($sample, $sampleDir, %dataFields);

sub importNimblegen{
    ($sample) = @_;
    status("loading Nimblegen array sample $sample...\n");
    $sampleDir = "$param{inputPath}/Array/$sample";
    my $arrayTable = newTable('Array', arrayTableName($sample, 'Nimblegen', $param{refSeqBase}));    
    my $outFile = "$arrayTable.csv"; 
    open my $outFileH, ">", $outFile;    
    foreach my $inFile(<$sampleDir/*unavg_segMNT.txt>){ 
        open my $inFileH, "<", $inFile;
        getArrayFields($inFileH, \%dataFields);      
        writeNimblegenData($inFileH, $outFileH);
        close $inFileH;
    }  
    close $outFileH;     
    loadData($outFile, $arrayTable, ',', $fieldNames{Array});      
}

sub writeNimblegenData{
    my ($inFileH, $outFileH) = @_;
    my $noValue = -1;
#    my $isSacCer2 = ($param{refSeqBase} eq 'sacCer2');
#    my %fixSacCer2 = (chr1=>'chrI',chr2=>'chrII',chr3=>'chrIII',chr4=>'chrIV',chr5=>'chrV',
#                      chr6=>'chrVI',chr7=>'chrVII',chr8=>'chrVIII',chr9=>'chrIX',chr10=>'chrX',
#                      chr11=>'chrXI',chr12=>'chrXII',chr13=>'chrXIII',chr14=>'chrXIV',chr15=>'chrXV',
#                      chr16=>'chrXVI',chr17=>'chrXVII',chr18=>'chrXVIII',chr19=>'chrXIX',chr20=>'chrXX');
    while (my $line = <$inFileH>){
        chomp $line; 
        $line =~ s/\r//g;
        my @fields = split("\t", $line);
        my $pos = $fields[$dataFields{CHR_POSITION}];                  
        $pos =~ m/^\d+$/ or next;   
        my $chr = $fields[$dataFields{CHROMOSOME}];
#        $isSacCer2 and $chr = $fixSacCer2{$chr};
        my $chrom = $refSeqs{$param{refSeqBase}}{$chr};  
        $chrom or next;
        my $ratio = 2 ** $fields[$dataFields{'RATIO_CORRECTED'}]; 
        print $outFileH join(",", $chrom, $pos, 
                              $ratio, $ratio, $noValue, $noValue, $noValue, 
                              $noValue, 1, $noValue, $noValue, $noValue, 
                              $ratio, $noValue)."\n";
#                             CHROMOSOME, POSITION, 
#                             RATIO, NORMALIZEDRATIO, BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, 
#                             CNVVALUE, GCSCORE, X, Y, THETA, 
#                             R, BF
    }
}

#Nimblegen unavg_segMNT.txt fields
#============================================================================
#CHROMOSOME
#CHR_POSITION
#RAW_DATAPOINTS
#POSITION_COUNT
#RATIO_CORRECTED
#WINDOW_SIZE
#============================================================================

1;

