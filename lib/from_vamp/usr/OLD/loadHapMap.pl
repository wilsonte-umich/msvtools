#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames %bowtieFields));

my ($runDir, @samples, $fileH, %dataFields);

sub loadHapMap{

    my %probes;
    runSQL("SELECT SNPNAME, CHROMOSOME, POSITION FROM HUMANOMNI1QUAD_HG18", \my($n,$c,$p));
    while ( fetchRow() ) { $probes{$n} = [$c,$p]; }

    my $file = "/home/wilsonte/vamp/data/Array/HapMap_CEU/GSE17197_HumanOmni1-Quad_v1_88CEU_FinalReport.csv";
    open $fileH, "<", $file;
    scrollToIlluminaHeaderYY();
    getArrayFieldsYY($fileH, \%dataFields);
    my %samples;
    while (my $line = <$fileH>){
        chomp $line; 
        $line =~ s/\r//g;
        my @fields = split(",", $line);  
        my $sample = $fields[$dataFields{"Sample ID"}];   
        $sample or next;
        $sample =~ s/ /_/g;
        $samples{$sample} or open $samples{$sample}, ">", "$sample.csv";
        my ($chrom, $pos) = @{$probes{$fields[$dataFields{'SNP Name'}]}};
        my $ratio = 2 ** $fields[$dataFields{'Log R Ratio'}];        
        my $bFreq = $fields[$dataFields{'B Allele Freq'}];                 
        my $zygosity = abs(0.5 - $bFreq) + 0.5;  
        my $cnvValue = -1;
        my $fileH = $samples{$sample};
        print $fileH join(",", 
                              $chrom, $pos,
                              $ratio, $ratio, $bFreq, $zygosity, $zygosity, 
                              $cnvValue, $fields[$dataFields{'GC Score'}], 
                              $fields[$dataFields{'X'}], $fields[$dataFields{'Y'}], -1, 
                              -1, $bFreq)."\n";
#                             CHROMOSOME, POSITION, 
#                             RATIO, NORMALIZEDRATIO, BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, 
#                             CNVVALUE, GCSCORE, X, Y, THETA, 
#                             R, BF
    }
    close $fileH;
    foreach my $sample (keys %samples) {
        my $fileH = $samples{$sample};
        close $fileH;
        my $arrayTable = newTable('Array', arrayTableName("HM_$sample", 'Illumina', "hg18"));
        loadData("$sample.csv", $arrayTable, ",", $fieldNames{Array});
    } 
}

sub scrollToIlluminaHeaderYY{
    my $identifier;
    while(!$identifier or !($identifier =~ m/\[Data\]/)){
        my $line = <$fileH>; 
        ($identifier) = split(",", $line);
    }
}

sub getArrayFieldsYY{
    my ($fileH, $fieldsRef) = @_;
    my $line = <$fileH>;
    chomp $line;
    $line =~ s/\r//g;
    my @labels = split(",", $line);
    for my $i(0..((scalar @labels) - 1)){$$fieldsRef{$labels[$i]} = $i}
}



1;



