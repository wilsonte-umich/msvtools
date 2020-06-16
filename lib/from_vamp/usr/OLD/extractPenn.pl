#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %bowtieFields));

##########################
##need vamp.pl to assemble the following job threads
#mapFixed, one thread for each read
#prepareUnfixed, one thread total
#mapUnfixed, one thread for each read
#mergeFixedMaps (can just cascade from mapUnfixed)
##########################



sub extractPenn{
#    my ($sample, $read) = @_;

    my $inFile = "/home/wilsonte/vamp/penncnv/7974.rawcnv";
    open my $inFileH, "<", $inFile;
    my $outFile = "outFile.csv";
    open my $outFileH, ">", $outFile;
    

#chr1:17548878-17549783        numsnp=15     length=906         state1,cn=0 /home/wilsonte/vamp/data/Array/Prjt_65/pennCNV/7974.txt startsnp=cnvi0055644 endsnp=cnvi0055655

    while (my $line = <$inFileH>){
        chomp $line; 
        $line =~ s/\r//g;
        while($line =~ s/  / /g){}
        my @fields = split(" ", $line);        
        my $chrPos = $fields[0];
        $chrPos =~ m/(.*):(.*)-(.*)/ or next;
        my $chr = $1;
        my $start = $2;
        my $end = $3;
        $chr =~ m/chr(.*)/ or next;
        my $chrom = $1;
        my $stateCN = $fields[3];
        $stateCN =~ m/.*,cn=(.*)/ or next;
        my $copyNumber = $1;
        print $outFileH join(",", ($chrom, $start, $end, $copyNumber))."\n";  
    }
    close $inFile;
    close $outFileH;

    
    loadData($outFile, "CNVS_7974_PENNCNV_HG18", ",", "CHROMOSOME, START_, END_, COPYNUMBER"); 

    
}


#this converted Illumina raw data to penncnv input

#sub extractPenn{
##    my ($sample, $read) = @_;

#    my $inFile = "/home/wilsonte/vamp/data/Array/Prjt_65/Prjt_65_Glaser_GLOVER_FinalReport_7974.txt";
#    open my $inFileH, "<", $inFile;
#    my $outFile = "/home/wilsonte/vamp/data/Array/Prjt_65/pennCNV/7974.txt";
#    open my $outFileH, ">", $outFile;
#    my $header = <$inFileH>;
#    
#    
#    my $identifier;
#    while(!$identifier or !($identifier =~ m/\[Data\]/)){
#        my $line = <$inFileH>; 
#        ($identifier) = split("\t", $line);
#    }
#    my %dataFields;
#    my $line = <$inFileH>;
#    chomp $line;
#    $line =~ s/\r//g;
#    my @labels = split("\t", $line);
#    for my $i(0..((scalar @labels) - 1)){$dataFields{$labels[$i]} = $i}
#    
#    print $outFileH join("\t", ("Name", "7974.Log R Ratio", "7974.B Allele Freq"))."\n";

#    while (my $line = <$inFileH>){
#        chomp $line; 
#        $line =~ s/\r//g;
#        my @fields = split("\t", $line);
#        my $pos = $fields[$dataFields{Position}];                  
#        $pos =~ m/^\d+$/ or next;         
#        my $name = $fields[$dataFields{'SNP Name'}]; 
#        my $lrr = $fields[$dataFields{'Log R Ratio'}];        
#        my $bFreq = $fields[$dataFields{'B Allele Freq'}];       
#        print $outFileH join("\t", ($name, $lrr, $bFreq))."\n";
#    }
#    
#    
#    close $inFile;
#    close $outFileH;

#    
#}



1;



