#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#Find.pl finds Sets in a single lane or merged sample
#########################################################################

use vars(qw(%param %types %fields %refSeqs $uid));

sub find{
    my ($sample) = @_;

    newTable('Marks', $sample);
    
    if($param{unpaired}){ processDiscrepancies1($sample); exit }

    status("getting statistics...\n");
        my $statsTable = getTableName('Stats', $sample);
        getStatistics($statsTable, \my %stats);
        #user can override default set coverage if desired
        $param{minFragsSet1} or $param{minFragsSet1} = $stats{minCoverageFMap} / 2; #/2 to account for mono-allelic events
        $param{minFragsSet1} >= 2 or $param{minFragsSet1} = 2; 
        $param{maxFragsSet1} or $param{maxFragsSet1} = $stats{maxCoverageFMap};
        ($param{minFragsSet2}, $param{maxFragsSet2}, $param{minFragsSet}, $param{maxFragsSet}) =
            (999999, 0, $param{minFragsSet1}, $param{maxFragsSet1}); 

    status("finding sets...\n");
        $param{fragsTable1} = getTableName('Frags', $sample);
        $param{strataFilter} = '';
        runFind($sample);

    status("establishing mismatch strata of fragments in sets...\n");
        establishFragmentStrata($param{fragsTable1}, 1);
        
    status("re-finding sets within lowest mismatch strata...\n");
        $param{strataFilter} = ' AND NSETSFRAG = -1 ';
        runFind($sample);

    status("loading sets...\n");
        loadFSData(1);        

    status("calculating fraction of promiscuous pairs in sets...\n");
        calculateFracBadFrags(1); 

    status("filling set counts...\n");
        fillSetCounts($param{fragsTable1}, 1);  

    if ($param{rooted}){
        
        status("estimating fraction overlap...\n");
            calcRootedFracOverlap($sample);
    
        status("finding nearest neighbor rooted sets\n    Chr: ");
            findRootedNeighbors1($sample);
            
        # status("marking neighbors sets which contain root elements within overlap...\n");
            # markKnownRoots($sample);    
        
    } else {
    
        status("adjusting set limits to reflect only non-promiscuous fragments...\n");
            fixSetLimits($sample);    
    
        status("creating set size histogram...\n");
            createSetSizeHistogram(1); 
            
        unless($param{exome}){
            status("finding Z sets...\n");
                findZSets($sample);  
        }   
                 
        # status("writing BED file...\n");
            # makeSetsBED($sample);

        unless ($param{noDisc}){
            status("processing discrepancies...\n");
            processDiscrepancies1($sample);          
        }    
    }  
}

sub runFind{
    my ($sample) = @_;
    initializeTables($sample, 1);         
    $param{setID} = 1; #setID incremented over all fragment types
    if($param{rooted}){
        findSetsRooted1();
    } else {
        findSets1($types{Sets}{Deletion});
        $param{skipInsertions} or findSets1($types{Sets}{Insertion});
        findSets1($types{Sets}{Inversion});
        findSets1($types{Sets}{Duplication});
        unless($param{noDiffChrom}){findTranslocations1()}  
        #not currently finding singleReads since there are so many
    }
    closeFiles(1);
}

sub initializeTables{
    my ($sample, $sampleN) = @_;
    initializeTable('Sets', $sampleN, '', $sample);
    initializeTable('IDs', $sampleN, 'F', $sample);
    initializeTable('IDs', $sampleN, 'P', $sample);
}

sub initializeTable{
    my ($tableType, $sampleN, $suffix, $sample) = @_;
    my ($pT, $pTF, $pTFH) = getFSPNames($tableType, $sampleN, $suffix);
    $param{$pT} = newTable($tableType, "$sample$suffix");  
    $param{$pTF} = "$param{$pT}.csv";
    open $param{$pTFH}, ">", $param{$pTF};
}

sub getFSPNames{
    my ($tableType, $sampleN, $suffix) = @_;
    return (getFSPName($tableType, $sampleN, $suffix, 'T'),
            getFSPName($tableType, $sampleN, $suffix, 'TF'),
            getFSPName($tableType, $sampleN, $suffix, 'TFH'))
}

sub closeFiles{
    my ($sampleN) = @_;
    foreach my $tfh (getFSFileHs($sampleN)){ close $tfh }
}

sub getFSTables{
    my ($sampleN) = @_;
    return ($param{getFSPName('Sets', $sampleN, '', 'T')},
            $param{getFSPName('IDs', $sampleN, 'F', 'T')},
            $param{getFSPName('IDs', $sampleN, 'P', 'T')} )   
}

sub getFSFiles{
    my ($sampleN) = @_;
    return ($param{getFSPName('Sets', $sampleN, '', 'TF')},
            $param{getFSPName('IDs', $sampleN, 'F', 'TF')},
            $param{getFSPName('IDs', $sampleN, 'P', 'TF')} )   
}

sub getFSFileHs{
    my ($sampleN) = @_;
    return ($param{getFSPName('Sets', $sampleN, '', 'TFH')},
            $param{getFSPName('IDs', $sampleN, 'F', 'TFH')},
            $param{getFSPName('IDs', $sampleN, 'P', 'TFH')} )   
}

sub getFSPName{
    my ($tableType, $sampleN, $suffix, $nameType) = @_;
    return "$tableType\_$nameType\_$sampleN\_$suffix";
}

sub createSetSizeHistogram{
    my ($sampleN) = @_;
    my ($setsTable) =  getFSTables($sampleN);    
    my $histogramTable = newTable('h', $setsTable);   
    runSQL("INSERT INTO $histogramTable
            (SELECT SETTYPE SERIES, Abs(EVENTMEAN) X, Count(SETID) Y
                FROM $setsTable
                GROUP BY SETTYPE, EVENTMEAN)  ");
}


1;


