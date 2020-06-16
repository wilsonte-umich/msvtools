#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames %bowtieFields));

my $dir = "/home/wilsonte";
my @increments = (1,2,5,10,50,100);

sub decile{
    foreach my $sample ("F112R","B3R","B5R"){
        my %deciles;
        my %counts;     
        my $table = "ARRAY_$sample\_AFFY_MM9";
        my %fileHs;
        my @deciles;
        foreach my $increment(@increments){
            my $file = "$dir/$table.$increment.csv";
            open my $fileH, ">", $file;
            $fileHs{$increment} = $fileH;
        }
#        runSQL("SELECT Round(arr.BFREQUENCY,1) 
#                FROM (select * from $table) arr, 
#                     (select * from mask_castb6) inf
#                WHERE arr.chromosome = inf.chromosome
#                  AND arr.position = inf.position
#                  AND arr.CHROMOSOME <= 19 AND arr.CHROMOSOME != 2
#                ORDER BY arr.CHROMOSOME, arr.POSITION", \my$decile);   
        runSQL("SELECT Round(log(2,RATIO),1) 
                FROM $table
                WHERE CHROMOSOME <= 19 AND CHROMOSOME != 2
                ORDER BY CHROMOSOME, POSITION", \my$decile);            
                          
        while (fetchRow()){push @deciles, $decile}
        foreach my $index(0..(scalar(@deciles)-101)){
            foreach my $increment(@increments){
                my $indexDecile = $deciles[$index];
                $deciles{$increment}{$indexDecile}{$deciles[$index + $increment]}++;
                $counts{$increment}{$indexDecile}++;
            }   
        }

        foreach my $increment(@increments){
            my $fileH = $fileHs{$increment};
            foreach my $indexDecile(sort {$a <=> $b} keys %{$deciles{$increment}}){ 
                print $fileH "$indexDecile,";
                foreach my $decile(sort {$a <=> $b} keys %{$deciles{$increment}}){ 
                    my $freq = 0;
                    $deciles{$increment}{$indexDecile}{$decile} 
                        and $freq = $deciles{$increment}{$indexDecile}{$decile} / $counts{$increment}{$indexDecile};
                    print $fileH "$freq,";  
                }
                print $fileH "\n";
            }        
        
        }

    }
}


1;



