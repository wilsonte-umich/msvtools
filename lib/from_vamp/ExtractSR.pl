#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#ExtractSR.pl pulls just the reads (without any pair information) from
#map files as "single reads".  A lookup is made to prefer reads that correspond
#to "preferred" (e.g. Normal) Fragments.  Other reads (mainly true single reads)
#are then stratified to keep the lowest number of discrepancies.  
#Finally read map tables are generated for subsequent discrepancy finding.
#########################################################################

use vars(qw(%param %types %fields %refSeqs));

sub extractSR{
    my ($sample) = @_;
    my ($inputSample, $outputSample) = splitSampleName($sample, 1);
    $outputSample or die "extractSR requires different input and output sample names in form inputSample::outputSample";
    #alternative: force name of $outputsample to $inputSampleSR
    
	getMapFiles($inputSample, \my%mapFilesIn);   
	my $fragsTableIn = getTableName('Frags', $inputSample);           
	my $fragsTableOut = newTable('Frags', $outputSample);     

    status("extracting mapped reads...\n");    
        loadMapsSR(\%mapFilesIn, $fragsTableOut);       

	status("finding best reads...\n");  
        markSRPreferred($fragsTableIn, $fragsTableOut);
        purgeSRWithPreferred($fragsTableOut);    
        pickBestNonPreferred($fragsTableOut, $outputSample);
		                            
    status("creating read map...\n"); 
        my $readMapTable = newTable('RMap', $outputSample);    
        createReadMapSR($fragsTableOut, $readMapTable);      
        my $histogramTable = newTable('h', $readMapTable);
        createMapHistogram($histogramTable, $readMapTable);          			
}    

#sub loadMapsSR{
#    my ($mapFilesRef, $fragsTable) = @_;
#	my $fragsFile = "$fragsTable.csv";           
#	open my $fragsFileH, ">", $fragsFile or die "could not open $fragsFile"; 
#    my $timeStamp = getTimeStamp(); 
#    loadMapsSR_Read($$mapFilesRef{read1}, $fragsFileH, 1, $timeStamp);  
#    loadMapsSR_Read($$mapFilesRef{read2}, $fragsFileH, 2, $timeStamp);    	
#	close $fragsFileH; 
#	loadReadsSR($fragsFile, $fragsTable);  
#}

#sub loadMapsSR_Read{
#    my ($mapFile, $fragsFileH, $readN, $timeStamp) = @_;
#    status("  extracting read $readN\n");
#        my $getLineSub = "getLine\_$param{mapType}";
#        my $getLineSubRef = \&{$getLineSub};
#        open my $mapFileH, "<", $mapFile or return;
#        my ($pairID, $map, $prevPairID, @prevMaps);
#        my $line = <$mapFileH>;
#        ($prevPairID, @prevMaps) = &$getLineSubRef($line);         
#        while (<$mapFileH>){     
#            ($pairID, $map) = &$getLineSubRef($_);
#            if ($pairID eq $prevPairID){
#                push @prevMaps, $map,
#            } else {
#                printReadsSR($fragsFileH, $timeStamp, $prevPairID, \@prevMaps, $readN);              
#                @prevMaps = ($map);
#            }
#            $prevPairID = $pairID;      
#        }
#        printReadsSR($fragsFileH, $timeStamp, $prevPairID, \@prevMaps, $readN);
#        close $mapFileH; 
#}

#sub printReadsSR{
#    my ($fragsFileH, $timeStamp, $pairID, $mapsRef, $readN) = @_;
#    $pairID = ($pairID * 1E8) + $timeStamp; #pairid same as extractPairs, allows cross correlation of usual to SR Frags tables
#    my $readID = ($pairID * 10) + $readN; #fragid carries readid which include read#
#    my $nReads = scalar @$mapsRef;   
#    foreach my $readMap(@$mapsRef){  
#        my ($chrom, $pos, $length, $strand, $discs) = split(/:/, $readMap);    
#        ##############################
#        #these start as type=singleRead => some will stay this way (others will be overridden with Normal etc)
#        #must make sure that collateDiscrepancies, findHomo, LODplot will accept singleRead!
#        printReadSR($fragsFileH, $readID, $types{Frags}{SingleRead},
#                                 $chrom, $pos, $length, $strand, $discs, 
#                                 $pairID, $nReads);
#    }  
#}

#sub printReadSR{
#    my ($fragsFileH, $readID, $fragType,
#                     $chrom, $pos, $length, $strand, $discs, 
#                     $pairID, $nReads) = @_;
#    print $fragsFileH join(",", $readID, $fragType,
#                                $chrom, $pos, $length, $strand, $discs, 
#                                $pairID, $nReads )."\n"; 
#}

#sub loadReadsSR{
#    my ($fragsFile, $fragsTable) = @_;
#    status("  uploading read data...\n");
#    loadData($fragsFile, $fragsTable, ",", "FRAGMENTID, FRAGMENTTYPE, 
#                                                  CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1,
#					                              PAIRID, NFRAGS");
#}

#sub markSRPreferred{
#    my ($fragsTableIn, $fragsTableOut) = @_;
#    status("  updating fragment type in preferred pairs (Normal, ReverseNormal, in Set)...\n");  
#    my $parentPair = "Trunc(PAIRID/1E8)";
#    my $inSQL1 = getSRPreferredInSQL($fragsTableIn, 1, $parentPair);
#    my $inSQL2 = getSRPreferredInSQL($fragsTableIn, 2, $parentPair);
#    my $inSQL = " $inSQL1 UNION ALL $inSQL2 ";   
#    my $outSQL = "SELECT $parentPair PARENTPAIR, FRAGMENTID, FRAGMENTTYPE,
#                         CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1,
#                         CHROMOSOME2, POSITION2, LENGTH2, STRAND2, DISCREPANCIES2,
#                         FRAGMENTSIZE, PAIRID, NFRAGS,
#                         EVENTSIZE, STDEVNORMAL, ENDTOLERANCE, 
#                         NSETSFRAG, NSETSPAIR
#                  FROM $fragsTableOut";                    
#    my $joinSQL = "SELECT o.FRAGMENTID, nvl(i.FRAGMENTTYPE, o.FRAGMENTTYPE) FRAGMENTTYPE,
#                           o.CHROMOSOME1, o.POSITION1, o.LENGTH1, o.STRAND1, o.DISCREPANCIES1,
#                           o.CHROMOSOME2, o.POSITION2, o.LENGTH2, o.STRAND2, o.DISCREPANCIES2,
#                           o.FRAGMENTSIZE, o.PAIRID, o.NFRAGS,
#                           o.EVENTSIZE, o.STDEVNORMAL, o.ENDTOLERANCE, 
#                           o.NSETSFRAG, o.NSETSPAIR
#                    FROM ($inSQL) i, ($outSQL) o
#                    WHERE i.PARENTPAIR(+) = o.PARENTPAIR 
#                      AND i.CHROMOSOME(+) = o.CHROMOSOME1
#                      AND i.POSITION(+) = o.POSITION1
#                      AND i.LENGTH(+) = o.LENGTH1
#                      AND i.STRAND(+) = o.STRAND1 ";
#    updateTable('Frags', $fragsTableOut, $joinSQL);    
#}

#sub getSRPreferredInSQL{
#    my ($fragsTableIn, $readN, $parentPair) = @_;
#    #returns unique preferred fragment for each PARENTPAIR
#    my $inSQL1 = "SELECT $parentPair PARENTPAIR, FRAGMENTTYPE,
#                         CHROMOSOME$readN CHROMOSOME, POSITION$readN POSITION, LENGTH$readN LENGTH, STRAND$readN STRAND
#                  FROM $fragsTableIn
#                  WHERE FRAGMENTTYPE = $types{Frags}{Normal}
#                     OR FRAGMENTTYPE = $types{Frags}{ReverseNormal} 
#                     OR (NSETSPAIR = 1 AND NSETSFRAG = 1)";
#                     #be careful never to ADD 'OR FRAGMENTTYPE = $types{Frags}{SingleRead}' to this!!!!
#                     #since this is where single read status is getting assigned!
#}

#sub purgeSRWithPreferred{
#    my ($fragsTable) = @_;
#    #alternative: purge ALL including those in SR preferred
#    status("  deleting non-preferred reads from same pair as a preferred read...\n");  
#    my $hasPreferredSQL = "SELECT FRAGMENTID
#                           FROM $fragsTable
#                           WHERE FRAGMENTTYPE != $types{Frags}{SingleRead}
#                           GROUP BY FRAGMENTID";
#    runSQL("DELETE FROM $fragsTable
#            WHERE FRAGMENTTYPE = $types{Frags}{SingleRead}
#              AND FRAGMENTID IN ($hasPreferredSQL) ");
#}

#sub pickBestNonPreferred{
#    my ($fragsTable, $sample) = @_;
#    status("  picking best read from remaining non-preferred reads\n");    
#    my $fragsTableTMP = newTable('Frags', "$sample\_TMP");     
#    my $fragsFileTMP = "$fragsTableTMP.csv";   
#    open my $fragsFileTMPH, ">", $fragsFileTMP; 
#    runSQL("SELECT FRAGMENTID, FRAGMENTTYPE, 
#                   CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1, 
#                   PAIRID, NFRAGS 
#            FROM $fragsTable
#            WHERE FRAGMENTTYPE = $types{Frags}{SingleRead}
#            ORDER BY FRAGMENTID"); 
#    my @readRefs;
#    my $prevReadID = 0; #readID = fragID       
#    while (my $readRef = fetchRowHashRef()){     
#        if ($$readRef{FRAGMENTID} > $prevReadID and $prevReadID){
#             processNonPreferred($fragsFileTMPH, \@readRefs);           
#            @readRefs = ();
#        }                     
#        push @readRefs, $readRef;
#        $prevReadID = $$readRef{FRAGMENTID};
#    } 
#    processNonPreferred($fragsFileTMPH, \@readRefs);  
#    close $fragsFileTMPH;
#    loadReadsSR($fragsFileTMP, $fragsTableTMP);
#    runSQL("DELETE FROM $fragsTable
#            WHERE FRAGMENTTYPE = $types{Frags}{SingleRead} ");  
#    runSQL("INSERT INTO $fragsTable SELECT * FROM $fragsTableTMP");      
#    dropTable($fragsTableTMP);    
#}         

#sub processNonPreferred{
#    my ($fragsFileH, $readRefsRef) = @_;
#    scalar @$readRefsRef == 1 and return printBestReadsSR($fragsFileH, $readRefsRef);
#    $readRefsRef = stratifyReadsSR(\&getNDisc, $readRefsRef);
#    scalar @$readRefsRef == 1 and return printBestReadsSR($fragsFileH, $readRefsRef);
#    $readRefsRef = stratifyReadsSR(\&getNInternalDiscs, $readRefsRef);     
#    printBestReadsSR($fragsFileH, $readRefsRef);   #is possible to have two "frags" from the same read! (but should be rare)      
#}

#sub printBestReadsSR{
#    my ($fragsFileH, $readRefsRef) = @_;
#    foreach my $readRef(@$readRefsRef){  
#        printReadSR($fragsFileH, $$readRef{FRAGMENTID}, $$readRef{FRAGMENTTYPE},
#                                 $$readRef{CHROMOSOME1}, $$readRef{POSITION1}, $$readRef{LENGTH1}, $$readRef{STRAND1}, $$readRef{DISCREPANCIES1}, 
#                                 $$readRef{PAIRID}, $$readRef{NFRAGS});
#    }  
#}

#sub stratifyReadsSR{
#    my ($nDiscsSubRef, $readRefsRef) = @_;
#    my %strata; 
#    foreach my $readRef(@$readRefsRef){ 
#        my $nDiscs = &$nDiscsSubRef($$readRef{DISCREPANCIES1}, $$readRef{LENGTH1});
#        push @{$strata{$nDiscs}}, $readRef;
#    }
#    my @strata = sort {$a <=> $b} keys %strata;
#    my $minDiscs = $strata[0]; 
#    return $strata{$minDiscs};
#}

#sub createReadMapSR{
#    my ($fragsTable, $readMapTable) = @_;
#    foreach my $chrom (1..$refSeqs{$param{refSeqBase}}{nChrom}){
#        my %coverage = ();
#        runSQL("SELECT ((trunc(POSITION1/$param{binSize})*$param{binSize})+$param{binSize}) AS LOWBIN1,
#                        (trunc((POSITION1 + LENGTH1 - 1)/$param{binSize})*$param{binSize}) AS HIGHBIN1
#                FROM $fragsTable
#                WHERE CHROMOSOME1 = $chrom",
#                \my($lowBin1, $highBin1));
#        while (fetchRow()){ fillRead(\%coverage, $lowBin1, $highBin1) }
#        my $readMapFile = "$readMapTable.csv";
#        open my $readMapFileH, ">", $readMapFile;
#        foreach my $pos (keys %coverage){ print $readMapFileH join(",", ($chrom, $pos, $coverage{$pos}))."\n" }
#        close $readMapFileH;
#        loadData($readMapFile, $readMapTable, ",", "CHROMOSOME, POSITION, COVERAGE");
#    }
#}

1;


