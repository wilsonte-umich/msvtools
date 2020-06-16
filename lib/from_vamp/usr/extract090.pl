#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames %bowtieFields));

sub extract090{

    my $dir = "/home/wilsonte/vamp/data/Array/Prjt_65";

    my @files = (
    "Prjt_65_Glaser_GLOVER_FinalReport_7974.txt",
    "Prjt_65_Glaser_GLOVER_FinalReport_7975.txt"
    );
    
    foreach my $file (@files) {
        open my $inFileH, "<", "$dir/$file" or die "crap";
        open my $outFileH, ">", "$dir/$file.geo" or die "crap";
        scrollToIlluminaHeaderYY($inFileH);
        getArrayFieldsYY($inFileH, \my%dataFields);
        print $outFileH "ID_REF\tGC_SCORE\tTheta\tR\tB_Allele_Freq\tLog_R_Ratio\n";
        while (my $line = <$inFileH>){
            chomp $line; 
            $line =~ s/\r//g;
            my @fields = split("\t", $line);
            print $outFileH
                join("\t",
                ($fields[$dataFields{"SNP Name"}],
                $fields[$dataFields{"GC Score"}],
                $fields[$dataFields{"Theta"}],
                $fields[$dataFields{"R"}],
                $fields[$dataFields{"B Allele Freq"}],
                $fields[$dataFields{"Log R Ratio"}]))."\n";
        }
        close $inFileH;
        close $outFileH; 
    }
}


sub scrollToIlluminaHeaderYY{
    my ($fileH) = @_;
    my $identifier;
    while(!$identifier or !($identifier =~ m/\[Data\]/)){
        my $line = <$fileH>; 
        ($identifier) = split("\t", $line);
    }
}

sub getArrayFieldsYY{
    my ($fileH, $fieldsRef) = @_;
    my $line = <$fileH>;
    chomp $line;
    $line =~ s/\r//g;
    my @labels = split("\t", $line);
    for my $i(0..((scalar @labels) - 1)){$$fieldsRef{$labels[$i]} = $i}
}

1;



