#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));

my @pairIDs = (
6412633,
15770992,
297722,
16853003,
1863060,
297362,
8463676,
3777121,
19708671,
7399492,
5556034,
13650535,
20037654,
814656,
12947213,
20831711,
876885,
2889103,
14635340,
1794363,
6769919,
732334,
);




my %pairIDs;

sub pullRooted2{
    foreach my $pairID(@pairIDs){$pairIDs{$pairID} = []}
    my $root = "L1";    
    foreach my $sample(qw(D2a D2b 090a 090b)){
        my $directory = "$param{inputPath}/$sample\_$root";
        my $pairedReadsFile = "$directory/$sample\_$root\_$root\_pairedReads.csv";    
        open my $inputFileH, "<", $pairedReadsFile;  
        while(my $line = <$inputFileH>){
            chomp $line;
            my ($inputPairID, $root, $flanking) = split(",", $line);#pairID,alignedSequence,partnerSequence 
            if($pairIDs{$inputPairID}){  
                push @{$pairIDs{$inputPairID}}, "$inputPairID\t$sample\t$flanking\t$root\n";    
            }
        }
        close $inputFileH;
    }
    foreach my $pairID(@pairIDs){ print join('', @{$pairIDs{$pairID}}) }
}




1;



