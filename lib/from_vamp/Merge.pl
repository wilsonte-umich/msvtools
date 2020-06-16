#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#Merge.pl contains the subs to merge Reads, Pairs or Fragments
#of a list of input samples into a single output sample of a new name.
#########################################################################

use vars(qw(%param %types %fields %refSeqs));

sub mergeReads {
    my($sampleOut, @samplesIn) = @_;
    
    unless (scalar @samplesIn) {die "no samples to merge\n"}
    
    status("creating output directories...\n");
        my %dirsOut;
        getDirectories($sampleOut, \%dirsOut);
        unless(-d $dirsOut{sample}) {mkdir $dirsOut{sample}}
        unless(-d $dirsOut{read1}) {mkdir $dirsOut{read1}}
        unless(-d $dirsOut{read2}) {mkdir $dirsOut{read2}}
    
    status("collecting input read files...\n");
        my @read1Files;
        my @read2Files;
        foreach my $sampleIn (@samplesIn) {
            my %readInFiles;
            getReadFiles($sampleIn, $param{readType}, \%readInFiles);
            push @read1Files, $readInFiles{read1};
            push @read2Files, $readInFiles{read2};   
        }
    
    status("merging input read files...\n");
        my %readOutFiles;
        getReadFiles($sampleOut, $param{readType}, \%readOutFiles);
        system("cat @read1Files > $readOutFiles{read1}");
        system("cat @read2Files > $readOutFiles{read2}");
    
    if ($param{dropInput}){
        status("deleting input read files...\n");
            unlink @read1Files;
            unlink @read2Files; 
    }
}

sub mergePairs{
    my ($sampleOut, @samplesIn) = @_;
    
    unless (scalar @samplesIn) {die "no samples to merge\n"}
    
    status("creating output table...\n");
        my $pairsTableOut = newTable('Pairs', $sampleOut);
    
    status("merging input pairs tables...\n");
        foreach my $sampleIn (@samplesIn){
            my $pairsTableIn = getTableName('Pairs', $sampleIn); 
            runSQL("INSERT INTO $pairsTableOut SELECT * FROM $pairsTableIn");
            if ($param{dropInput}){
                status("dropping $sampleIn...\n");
                    dropTable($pairsTableIn);
                    dropTable(getTableName('Stats', $sampleIn));                            
            }  
        }
    
    status("creating pair size histogram...\n");
    my $histogramTable = newTable('h', $pairsTableOut);
    createPairSizeHistogram($histogramTable, $pairsTableOut);

    status("calculating pair statistics...\n");
    my $statsTableOut = newTable('Stats', $sampleOut);
	calculateNormalStats($histogramTable, $statsTableOut);
}

sub mergeFragments{
    my ($sampleOut, @samplesIn) = @_;
       
    unless (scalar @samplesIn) {die "no samples to merge\n"}
 
    status("creating output tables...\n");
        my $fragsTableOut = newTable('Frags', $sampleOut);    
        my $statsTableOut = newTable('Stats', $sampleOut);    
        my ($hFragsTableOut, $fMapTableOut, $rMapTableOut);
        unless($param{rooted}){
            $hFragsTableOut = newTable('h', $fragsTableOut);
            $fMapTableOut = newTable('FMap', $sampleOut);
            unless($param{noDisc}){ $rMapTableOut = newTable('RMap', $sampleOut) }        
        }

    status("merging input fragments tables...\n");
        my $normalSum = 0;
        my $reverseNormalSum = 0;
        my $fragCoverageSum = 0;
        my $readCoverageSum = 0;
        my $minNormal = 1E9;
        my $maxNormal = 0;
        my $minInsertion = 1E9;
        my $minDeletion = 1E9;
        foreach my $sampleIn (@samplesIn){
            my $fragsTableIn = getTableName('Frags', $sampleIn);
            runSQL("INSERT INTO $fragsTableOut SELECT * FROM $fragsTableIn");
            my (%stats, $statsTableIn, $hFragsTableIn, $fMapTableIn, $rMapTableIn);            
            unless($param{rooted}){            
                $statsTableIn = getTableName('Stats', $sampleIn);
                getStatistics($statsTableIn, \%stats);
                $normalSum += $stats{normalCount};
                $reverseNormalSum += $stats{reverseNormalCount};
                $fragCoverageSum += $stats{fragCoverage};
                $readCoverageSum += $stats{readCoverage};
                $minNormal <= $stats{minNormal} or $minNormal = $stats{minNormal};
                $maxNormal >= $stats{maxNormal} or $maxNormal = $stats{maxNormal}; 
                $minInsertion <= $stats{minInsertion} or $minInsertion = $stats{minInsertion};
                $minDeletion <= $stats{minDeletion} or $minDeletion = $stats{minDeletion};
                $hFragsTableIn = getTableName('h', $fragsTableIn);            
                runSQL("INSERT INTO $hFragsTableOut SELECT * FROM $hFragsTableIn");
                $fMapTableIn = getTableName('FMap', $sampleIn);
                runSQL("INSERT INTO $fMapTableOut SELECT * FROM $fMapTableIn");
                if ($rMapTableOut){
                    $rMapTableIn = getTableName('RMap', $sampleIn);
                    runSQL("INSERT INTO $rMapTableOut SELECT * FROM $rMapTableIn");
                }            
            }
            if ($param{dropInput}){ #drop inputs
                status("dropping $sampleIn...\n");
                    dropTable($fragsTableIn);                
                    dropTable($statsTableIn);
                    dropTable($hFragsTableIn);
                    dropTable($fMapTableIn);
                    dropTable($rMapTableIn);
            }
        }

    if ($param{rooted}){
        status("copying stats from reference merged stats...\n");
            $sampleOut =~ m/(.+)\_$param{rooted}$/;
            $1 or return;
            my $refStatsTable = getTableName('Stats', $1);
            tableExists($refStatsTable) or die "could not find table $refStatsTable";
            updateTable('Stats', $statsTableOut, "SELECT * FROM $refStatsTable");
    } else {
        updateStat($statsTableOut, 'normalCount', $normalSum);
        updateStat($statsTableOut, 'reverseNormalCount', $reverseNormalSum);
        updateStat($statsTableOut, 'minNormal', $minNormal);
        updateStat($statsTableOut, 'maxNormal', $maxNormal);
        updateStat($statsTableOut, 'minInsertion', $minInsertion);
        updateStat($statsTableOut, 'minDeletion', $minDeletion);
        updateStat($statsTableOut, 'fragCoverage', $fragCoverageSum);
        updateStat($statsTableOut, 'readCoverage', $readCoverageSum);
        status("  $normalSum Normal Fragments\n");
        status("  $reverseNormalSum ReverseNormal Fragments\n");
        status("  combined Normal peaks range from $minNormal to $maxNormal bp\n");
        status("  minimum size of detectable Insertion = $minInsertion\n");
        status("  minimum size of detectable Deletion = $minDeletion\n");        
        status("  Fragment coverage = $fragCoverageSum\n");
        status("  Read coverage = $readCoverageSum\n");    

        status("combining size histogram and maps...\n");
            updateTable('h', $hFragsTableOut, 
                            "SELECT SERIES, X, Sum(Y) Y 
                             FROM $hFragsTableOut
                             GROUP BY SERIES, X"); 
            updateTable('FMap', $fMapTableOut, 
                        "SELECT CHROMOSOME, POSITION, Sum(COVERAGE) COVERAGE, 0 NORMALIZEDCOVERAGE
                         FROM $fMapTableOut
                         GROUP BY CHROMOSOME, POSITION");              
            $rMapTableOut and 
                updateTable('RMap', $rMapTableOut, 
                        "SELECT CHROMOSOME, POSITION, Sum(COVERAGE) COVERAGE, 0 NORMALIZEDCOVERAGE
                         FROM $rMapTableOut
                         GROUP BY CHROMOSOME, POSITION"); 
                         
        status("calculating map histograms...\n");
            my $hFMapTableOut = newTable('h', $fMapTableOut);
            createMapHistogram($hFMapTableOut, $fMapTableOut);
            calculateCoverageStats($hFMapTableOut, 'FMap', $statsTableOut);
            if ($rMapTableOut){
                my $hRMapTableOut = newTable('h', $rMapTableOut);
                createMapHistogram($hRMapTableOut, $rMapTableOut);
                calculateCoverageStats($hRMapTableOut, 'RMap', $statsTableOut);   
            }
        
        status("calculating normalized coverages...\n");
            $fMapTableOut and $fragCoverageSum and calculateNormalizedCoverage($fMapTableOut, 'FMap', $fragCoverageSum);
            $rMapTableOut and $readCoverageSum and calculateNormalizedCoverage($rMapTableOut, 'RMap', $readCoverageSum);                                               
    }   
}

sub copyFragments{
    my ($sampleOut, $sampleIn) = @_;
    copyTable('Frags', $sampleOut, $sampleIn);
    copyTable('Stats', $sampleOut, $sampleIn);
    copyTable('FMap', $sampleOut, $sampleIn); 
    #add RMap too?
}

sub copyTable{
    my ($tableType, $sampleOut, $sampleIn) = @_;
    status("copying $tableType table...\n");
        my $tableIn = getTableName($tableType, $sampleIn);
        my $tableOut = newTable($tableType, $sampleOut);
        runSQL("INSERT INTO $tableOut SELECT * FROM $tableIn"); 
}

1;
