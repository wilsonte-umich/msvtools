#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %bowtieFields));

my ($runDir, @samples, $fileH, %dataFields);

sub loadOmni{
    my $file = "/home/wilsonte/vamp/data/Array/Prjt_65/Prjt_65_Glaser_GLOVER_FinalReport_7974.txt";
    open $fileH, "<", $file;
    scrollToIlluminaHeaderXX();
    getArrayFieldsXX($fileH, \%dataFields);
    my $arrayTable = "HUMANOMNI1QUAD_hg18";
    my $arrayFile = "$arrayTable.csv";
    open my $arrayFileH, ">", $arrayFile;
    while (my $line = <$fileH>){
        chomp $line; 
        $line =~ s/\r//g;
        my @fields = split("\t", $line);
        my $pos = $fields[$dataFields{Position}];                  
        $pos =~ m/^\d+$/ or next;         
        my $chr = $fields[$dataFields{Chr}];      
        $chr eq 'XY' and $chr = 'X'; #????
        $chr eq 'MT' and $chr = 'M'; 
        my $chrom = $refSeqs{hg18}{"chr$chr"};         
        print $arrayFileH join(",",   $fields[$dataFields{"SNP Name"}],
                                      $chrom, 
                                      $pos, 
                                      $fields[$dataFields{"Allele1 - Top"}],
                                      $fields[$dataFields{"Allele2 - Top"}] )."\n";
    }
    close $fileH;
    close $arrayFileH;
    loadData($arrayFile, $arrayTable, ",", "SNPNAME, CHROMOSOME, POSITION, ALLELE1, ALLELE2");
}

sub scrollToIlluminaHeaderXX{
    my $identifier;
    while(!$identifier or !($identifier =~ m/\[Data\]/)){
        my $line = <$fileH>; 
        ($identifier) = split("\t", $line);
    }
}

sub getArrayFieldsXX{
    my ($fileH, $fieldsRef) = @_;
    my $line = <$fileH>;
    chomp $line;
    $line =~ s/\r//g;
    my @labels = split("\t", $line);
    for my $i(0..((scalar @labels) - 1)){$$fieldsRef{$labels[$i]} = $i}
}



1;



