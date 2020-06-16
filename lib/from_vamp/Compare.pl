#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#Compare.pl finds Sets by comparing sample1 to sample2
#########################################################################

use vars(qw(%param %types %fields %refSeqs $uid));

sub compare{
    my ($sample1, $sample2) = @_;

    newTable('Marks', $sample1);
    newTable('Marks', $sample2);
    
    if($param{unpaired}){ processDiscrepancies2($sample1, $sample2); exit }
             
    status("getting statistics...\n");
        my $statsTable1 = getTableName('Stats', $sample1);
        my $statsTable2 = getTableName('Stats', $sample2);
        getStatistics($statsTable1, \my %stats1);
        getStatistics($statsTable2, \my %stats2);
        #user can override default set coverage if desired
        $param{minFragsSet1} or $param{minFragsSet1} = $stats1{minCoverageFMap} / 2; #/2 to account for mono-allelic events
        $param{minFragsSet1} >= 2 or $param{minFragsSet1} = 2; 
        $param{maxFragsSet1} or $param{maxFragsSet1} = $stats1{maxCoverageFMap}; 
        $param{minFragsSet2} or $param{minFragsSet2} = $stats2{minCoverageFMap} / 2; 
        $param{minFragsSet2} >= 2 or $param{minFragsSet2} = 2; 
        $param{maxFragsSet2} or $param{maxFragsSet2} = $stats2{maxCoverageFMap};  
        $param{minFragsSet} = $param{minFragsSet1};
        $param{minFragsSet} <= $param{minFragsSet2} or $param{minFragsSet} = $param{minFragsSet2};        
        $param{maxFragsSet} = ($param{maxFragsSet1} + $param{maxFragsSet2});  

    status("finding sets...\n");
        $param{fragsTable1} = getTableName('Frags', $sample1);
        $param{fragsTable2} = getTableName('Frags', $sample2);
        $param{strataFilter} = '';
        runCompare($sample1, $sample2);

    status("establishing mismatch strata of fragments in sets...\n");
        establishFragmentStrata($param{fragsTable1}, 1);
        establishFragmentStrata($param{fragsTable2}, 2);
        
    status("re-finding sets within lowest mismatch strata...\n");
        $param{strataFilter} = ' AND NSETSFRAG = -1 ';
        runCompare($sample1, $sample2);

    status("loading sets...\n");
        loadFSData(1);
        loadFSData(2);
        
    status("calculating fraction of promiscuous fragments in sets...\n");
        calculateFracBadFrags(1);
        calculateFracBadFrags(2);
    
    status("filling set counts...\n");
        fillSetCounts($param{fragsTable1}, 1);
        fillSetCounts($param{fragsTable2}, 2);
        
        
    if ($param{rooted}){
        
        status("estimating fraction overlap...\n");
            calcRootedFracOverlap($sample1);
            calcRootedFracOverlap($sample2);

        status("finding nearest neighbor rooted sets\n    Chr: ");
            findRootedNeighbors1($sample1); #this maybe should find neighbors on combined sample frags?
            findRootedNeighbors1($sample2);

        # status("marking neighbors sets which contain root elements within overlap...\n");
            # markKnownRoots($sample);    
        
    } else {
          
        status("adjusting set limits to reflect only non-promiscuous fragments...\n");
            fixSetLimits($sample1); 
            fixSetLimits($sample2);           
            
        status("creating set size histograms...\n");
            createSetSizeHistogram(1); 
            createSetSizeHistogram(2);         

        status("finding Z sets...\n");
            findZSets($sample1);           
            findZSets($sample2);    

        ##this works but is only used for making a total genome plot, never accessed by downstream code
        ##the visualizer provides a more useful dynamic ratio plot anyway
        #status("calculating coverage ratio...\n");
        #    my $mapTable1 = getTableName('FMap', $sample1);
        #    my $mapTable2 = getTableName('FMap', $sample2);
        #    calculateCoverageRatio($mapTable1, $mapTable2, $sample1, $sample2);
        #    calculateCoverageRatio($mapTable2, $mapTable1, $sample2, $sample1);     

    #    status("writing BED files...\n");
    #        makeSetsBED($sample1);
    #        makeSetsBED($sample2);
        
        unless($param{noDisc}){
            status("processing discrepancies...\n");
                processDiscrepancies2($sample1, $sample2);      
        }  
    }
}

sub runCompare{
    my ($sample1, $sample2) = @_;
    initializeTables($sample1, 1);
    initializeTables($sample2, 2);
    $param{setID} = 1; #setID incremented over all fragment types
    if($param{rooted}){
        findSetsRooted2();
        #die "compare not yet implemented for rooted analysis";
    } else {
        findSets2($types{Sets}{Deletion});
        $param{skipInsertions} or findSets2($types{Sets}{Insertion});
        findSets2($types{Sets}{Inversion});
        findSets2($types{Sets}{Duplication});
        unless($param{noDiffChrom}){findTranslocations2()}
        #not currently finding singleReads since there are so many
    }
    closeFiles(1);
    closeFiles(2);
}

#sub calculateCoverageRatio{
#    my ($mapTable1, $mapTable2, $sample1, $sample2) = @_;
#    my $rowsPF1E4 = int((1E4/$param{binSize})/2);
#    my $rowsPF1E5 = int((1E5/$param{binSize})/2);
#    my $ratioTable = newTable('FRatio', "$sample1\_$sample2");
#    runSQL("INSERT INTO $ratioTable
#            SELECT m2.CHROMOSOME, m2.POSITION,
#                (nvl(m1.NORMALIZEDCOVERAGE,0)/m2.NORMALIZEDCOVERAGE) COVERAGERATIO,
#                0 COVERAGERATIO1E4, 0 COVERAGERATIO1E5
#            FROM $mapTable1 m1, $mapTable2 m2
#            WHERE m1.CHROMOSOME(+) = m2.CHROMOSOME
#                AND m1.POSITION(+) = m2.POSITION");
#    updateTable('FRatio', $ratioTable,
#                "SELECT CHROMOSOME, POSITION, COVERAGERATIO, 
#                    Avg(Sum(COVERAGERATIO))
#                        OVER (ORDER BY CHROMOSOME, POSITION 
#                                ROWS BETWEEN $rowsPF1E4 PRECEDING
#                                AND $rowsPF1E4 FOLLOWING) COVERAGERATIO1E4,
#                    Avg(Sum(COVERAGERATIO))
#                        OVER (ORDER BY CHROMOSOME, POSITION 
#                                ROWS BETWEEN $rowsPF1E5 PRECEDING
#                                AND $rowsPF1E5 FOLLOWING) COVERAGERATIO1E5                    
#                FROM $ratioTable
#                GROUP BY CHROMOSOME, POSITION, COVERAGERATIO
#                ORDER BY CHROMOSOME, POSITION");
#    #this moving average is NOT perfect since it goes by rows,
#    #and all positions are not represented in rows
#    #(only bin positions covered in sample2 have ratio table rows)
#    #could average over a position range but this is very slow
#}

1;

