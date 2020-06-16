#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %samFields %samFlags));

#callable parameters and commands in this script
defined $command{reduceNoise2} or $command{reduceNoise2} = ['multiThread', '48:00:00', 2000, 0];

sub reduceNoise2 { #this version only applies maxNFrags AND irreversibly deletes all frags from offending pair!
    my ($sample) = @_; 
    my $fragsTable = getTableName('Frags', $sample);

    #provide feedback on effectiveness of noise reduction
    status("anomalous fragment count\n");        
    my $startFragCount = countPersistentAnomalies($fragsTable);
    status("  before noise reduction\t$startFragCount\n"); 
    
    #DELETE all anomalous fragments from pairs with more than threshold number of total fragments
    reduceNoise_parameter2($fragsTable);
    my $endFragCount = countPersistentAnomalies($fragsTable);
    my $percentRemaining = int((($endFragCount / $startFragCount)*100)+0.5);
    status("  after maxNFrags\t$endFragCount ($percentRemaining"."\%)\n");         

}
sub reduceNoise_parameter2 { #irreversibly delete all anomalous fragments where source pairID violates maxNFrags
    my ($fragsTable) = @_;
    getRNParameter($fragsTable, 'maxNFrags', 'nFrags');
    runSQL("DELETE FROM $fragsTable
            WHERE fragmenttype > $types{Frags}{Normal}
              AND nFrags > $param{maxNFrags}");        
}

1;

