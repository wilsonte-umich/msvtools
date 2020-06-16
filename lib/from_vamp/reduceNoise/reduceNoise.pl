#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %samFields %samFlags));

###############################################################################
#'reduceNoise.pl' reduces the "noise" of low quality anomalous Fragments by 
#deleting anomalous fragments that have a high likelihood of being falsely anomalous.  
#It is tolerable if a few true anomalous Fragments are lost if coverage is high enough, 
#with the net result of improving Set finding.  
#----------------------------------------------------------------------------------
#Approach was motivated mainly by the need to look for SVs in non-barcoded pooled samples,
#where the noise is determined by the sample pool, but signal is determined by just one sample, 
#creating a situtation of low signal to noise ratio that greatly benefits from 
#a reduction of the noise of false anomalies.
#----------------------------------------------------------------------------------
#The logic of the analysis thus goes:
#----------------------------------------------------------------------------------
#existing steps:
#  1) keep initial maxDisc high during mapReads to capture 'all' potential Expected
#  2) find Expected based on these lenient criteria (parseFragments)
#  3) purge false anomalies with Expected from same pair (parseFragments)
#----------------------------------------------------------------------------------
#new reduceNoise steps:
#  4) determine the nFrags and nDiscs profiles of Expected Fragments
#     (Expected Fragments give robust assessment of bona fide read pair quality parameters)
#  5) select maxNFrags and maxNFragDiscs based on above, desiring to retain only high quality pairs
#  6) mask (don't delete) anomalous Fragments with more than maxNFrags and maxNFragDiscs 
#  7) the only anomalies that remain have similar mapping quality as Expected
#     (this necessarily limits anomalies to unique genome space!) 
###############################################################################


#callable parameters and commands in this script
defined $param{maxNFrags} or $param{maxNFrags} = -1; #max number of frags allowed per pair
defined $param{maxNFragDiscs} or $param{maxNFragDiscs} = -1; #max number of Discs allowed in both reads combined
defined $param{rnFractionExpected} or $param{rnFractionExpected} = 0.96; #% of Expected to account for if above not provided
defined $command{reduceNoise} or $command{reduceNoise} = ['multiThread', '48:00:00', 2000, 0];

sub reduceNoise {
    my ($sample) = @_; 
    my $fragsTable = getTableName('Frags', $sample);
    
    #revese any prior noise reduction
    reverseNoiseReduction($fragsTable); 

    #provide feedback on effectiveness of noise reduction
    status("anomalous fragment count\n");        
    my $startFragCount = countPersistentAnomalies($fragsTable);
    status("  before noise reduction\t$startFragCount\n"); 
    
    #mask all anomalous fragments from pairs with more than threshold number of total fragments
    reduceNoise_parameter($fragsTable, 'maxNFrags', 'nFrags');
    my $endFragCount = countPersistentAnomalies($fragsTable);
    my $percentRemaining = int((($endFragCount / $startFragCount)*100)+0.5);
    status("  after maxNFrags\t$endFragCount ($percentRemaining"."\%)\n");     
     
    #remove anomalous fragments with more than threshold number of total discrepancies
    my $nFragDiscsSQL = "ceil((length(discrepancies1)/4)-0.25) + ceil((length(discrepancies2)/4)-0.25)";
    reduceNoise_parameter($fragsTable, 'maxNFragDiscs', $nFragDiscsSQL); 
    $endFragCount = countPersistentAnomalies($fragsTable);
    $percentRemaining = int((($endFragCount / $startFragCount)*100)+0.5);
    status("  after maxNFragDiscs\t$endFragCount ($percentRemaining"."\%)\n");      

}
sub reverseNoiseReduction { #revese any prior noise reduction by re-negating any negative anomalous fragmentid
    my ($fragsTable) = @_;
    runSQL("UPDATE $fragsTable
            SET fragmentid = -fragmentid
            WHERE fragmenttype > $types{Frags}{Normal}
              AND fragmentid < 0");
}
sub countPersistentAnomalies { #simple anomaly count
    my ($fragsTable) = @_;
    my $nChrom = nChrom();
    my $maxFragmentType = $types{Frags}{DiffChrom};
    $param{noDiffChrom} and $maxFragmentType = $types{Frags}{Duplication};
    runSQL("SELECT count(*)
            FROM $fragsTable
            WHERE fragmenttype > $types{Frags}{Normal}
              AND fragmenttype <= $maxFragmentType
              AND fragmentid > 0",\my($count));
    fetchRow();
    $count or $count = 0;
    return $count;
}
sub reduceNoise_parameter { #mask low quality based on input parameter by negating fragmentid
    my ($fragsTable, $paramKey, $valField) = @_;
    getRNParameter($fragsTable, $paramKey, $valField);
    runSQL("UPDATE $fragsTable
            SET fragmentid = -fragmentid
            WHERE fragmenttype > $types{Frags}{Normal}
              AND fragmentid > 0
              AND $valField > $param{$paramKey}");        
}
sub getRNParameter {
    my ($fragsTable, $paramKey, $valField) = @_;
    unless($param{$paramKey} > -1){ #use the value set by user if present
        my $nChrom = nChrom(); #else use rnFractionExpected to dynamically determine maxNFragDiscs from Expected Fragments
        my $expectedSQL = "SELECT $valField val
                           FROM $fragsTable
                           WHERE chromosome1 <= $nChrom
                             AND fragmenttype = $types{Frags}{Normal}";  
        my $countSQL = "SELECT val, count(*) N
                        FROM ($expectedSQL)
                        GROUP BY val";
        runSQL($countSQL, \my($val,$count));
        my ($nFrags, %val, $cumFraction);
        while (fetchRow()){ #collect counts by nFragDiscs
            $val{$val} = $count;
            $nFrags += $count;
        }
        foreach my $val(keys %val){ $val{$val} /= $nFrags } #convert to a fraction
        foreach my $val(sort {$a <=> $b} keys %val){ #determine lowest nFragDiscs that satisfies rnFractionExpected
            $cumFraction += $val{$val};
            $cumFraction >= $param{rnFractionExpected} and $param{$paramKey} = $val and last;
        }
    }
    status("  using $paramKey = $param{$paramKey}\n");
}

1;




