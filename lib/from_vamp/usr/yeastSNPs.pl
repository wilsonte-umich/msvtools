#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

#--AC-CACACCCACACC-CAC chr10 19
#ACACCCACACACACACCACAC chr10 237 YS4 C>A 31 
#@@@@@@@@@@@M@@@@@@@@@

#CACAC--CA-ACTC-TCT--C chr10 61
#CACACAT--TTCACATTTCA- chr10 289 YS4 A>T 31 
#M@@@@@@GG@@@@@@@@@@@H

#CAC--CA-ACTC-TCT--CTC chr10 63
#CACAT--TTCACATTTCA--C chr10 291 YS4 T>A 31 
#@@@@@GG@@@@@@@@@@@HI@

my %chroms = (chr01=>1, chr02=>2, chr03=>3, chr04=>4, chr05=>5, 
              chr06=>6, chr07=>7, chr08=>8, chr09=>9, chr10=>10, 
              chr11=>11, chr12=>12, chr13=>13, chr14=>14, chr15=>15, 
              chr16=>16);

sub yeastSNPs{
    my $snpsTable = newTable("SNPs", "sacCerSNP");
    my $snpsFile = "$snpsTable.csv";
    open my $snpFileH, ">", $snpsFile;
    my $folder = "/home/wilsonte/vamp/refSeqs/sacCerSNP/";
#    my $counter = 0;
    foreach my $inFile (<$folder/*_sequencedSNPs.txt>){
        print "$inFile\n";
        open my $inFileH, "<", $inFile;
        while (my $line = getNextLine($inFileH)){
            $line =~ m/.{22}(.*)/;
            $1 or die "BAD LINE\n$line";
            my @fields = split(" ", $1);
            my $chrom = $chroms{$fields[0]};
            $line = <$inFileH>;
            if ($chrom){
                my $pos = $fields[1] + 1;
                $line =~ m/.{22}(.*)/;
                @fields = split(" ", $1);
                my $strain = $fields[2];
                my $snp = $fields[3];
                $snp =~ m/.>(.)/;
                if ($1){
                    my $disc = $types{Discs}{$1};
#                    print "$strain,$chrom,$pos,$disc\n";
                    print $snpFileH "$strain,$chrom,$pos,$disc\n";   
                } else {
                    die "BAD LINE:\n$line\n";
                }
            }
            $line = <$inFileH>;
#            $counter > 100 and exit;
#            $counter++;
        }
        close $inFileH;
    }
    close $snpFileH;
    loadData($snpsFile, $snpsTable, ",", $fieldNames{SNPs});
}

sub getNextLine{
    my ($fileH) = @_;
    my $line = undef;
    while (!$line) { 
        $line = <$fileH>;
        $line or return undef;
        $line = parseLine($line);
    }
    return $line;
}

sub parseLine{
    my ($line) = @_;
    $line or return undef;
    chomp $line;
    $line =~ s/\r//;
    $line or return undef;
    length($line) > 0 or return undef;
    return $line;
}

1;

