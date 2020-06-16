#!/usr/bin/perl
use strict;
use warnings;

my $file = "/home/GROUPS/wilsontelab/data/Ljungman/Run_198/nf0h1_mRna/1/nf0h1_mRna_1.bowtie";
my $chr = "chr1";
my $min = 160723000;
my $max = 160725000;

open my $fileH, "<", $file or die "could not open $file: $!";
my (%reads, %pairIDs, %maps);
while (my $line = <$fileH>){
#2019	-	chr15	39885857	ACACCCCTGGCCAGGTGCGCACCCTGAGGCATGACCCTCG		0	13:T>A
    my @line = split("\t", $line);
    $line[2] eq $chr or next;
    $line[3] > $min or next;
    $line[3] < $max or next;
    push @{$reads{$line[4]}}, $line;
    $pairIDs{$line[0]}++;
}
close $fileH;
foreach my $read(keys %reads){
    print $read, "\n\t", join("\n\t", @{$reads{$read}}), "\n";
}

$file = "/home/GROUPS/wilsontelab/data/Ljungman/Run_198/nf0h1_mRna/1/nf0h1_mRna_1.bowtie.mRna";
open $fileH, "<", $file or die "could not open $file: $!";
while (my $line = <$fileH>){
    my @line = split("\t", $line);
    defined $pairIDs{$line[0]} or next;
    push @{$maps{$line[0]}}, $line;
}
close $fileH;
foreach my $pairID(keys %maps){
    print $pairID, "\n\t", join("\n\t", @{$maps{$pairID}}), "\n";
}


#    AAGGCATGAGAATCGCTTGAACCTGGGAGGTGGAGGCTGC	160724144 160724144
#    CCCAGCATGGTGGCTCACACCTGTAATCCCAGCACTTTGC	160723963
#    CCGGCATGGTGGCTCACACCTGTAATCCCAGCACTTTGGC	160723964
#    TAGGCGGCCATGCATGGTGGCTCACACCTGTAATCCCAGC	160723956
#    TCTGAGTGTGGTGGTGTGCACCTGTAATCCCAGCTACTCG	160724097
#    CGAGGCCAAGGCACGAGAATCGCTTGAACCTGGGAGGTGG	160724137
#    GCTGAGTGTGGTGGTGTGCACCTGTAATCCCAGCTAGTCG	160724097
#    GGGTCATCTGAGGTCAGGAGTTCAAGACCAGCCTGCCCAC	160724018
#    CCAGCCTTGCCAACATGGTGAAACCCCATGTCTACTAAAG	160724045
#    CCGGCATGCTGGCTCACACCTGTAATCCCAGCACTTTTGG	160723964
#    GTTGTGCACCTGTATTCCCAGCTACTCGA	160724109
#    GCAGCTACTCGGGAGGCCAAGGCATGAGAATCGCTTGAAC	160724126
#    CCCCGCTTGGTGGCCCACACCTGTAATCCCAGCACTTTGG	160723963
#    TGAGAATCGCTTGAACCTGGGAGGTGGAGGTTGCAGTGAC	160724150
#    GTGTGCACCTGTAATCCCAGCTATTCGGGAGCCCAAGGCA	160724110
#    GGGAGGCCAAGGCATGAGAATCGCTTGAACCTGGGAGGTC	160724136
#    TAGGCATGAGAATCGCTTGAACCTGGGAGGAGGAGGTTGC	160724144
#    TGAGAATCGCTTGAGCCTGGGAGGTGGAGGTTGCCGTGAC	160724150
#    GCCCGGCATGGTGGCTCACACCTGTAATCCCAGCACTCTC	160723962
#    GTGTAGTCCCAGCACTTTGGGAGGCCAAGGTGGATGGGTC	160723983







1;
