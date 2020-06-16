#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#ExtractPairs.pl coordinates the conversion of text files of mapped reads
#into Oracle Pairs data tables, and then calculates Pair statistics.
#########################################################################

use vars(qw(%param %types %fields %refSeqs));

sub extractPairsSplit {
    my ($sample) = @_;
    my ($inputSample, $outputSample) = splitSampleName($sample);
     
    getDirectories($inputSample, \my%dirs);    
    my $paramFile = "$dirs{sample}/$outputSample\_extract\_parameters.txt";    
    printParameters($paramFile, 'extractPairsSplit', $sample);    

    status("extracting pairs...\n");
        my (@mapFiles1, @mapFiles2);
    	my %ext = (pass => 'gff', bowtie => 'bowtie');
        my $prepDir1 = $dirs{prepared_read1};
        my $prepDir2 = $dirs{prepared_read2};
        foreach my $groupDir1 (<$prepDir1/*>){
            -d $groupDir1 or next;
            $groupDir1 =~ m/$prepDir1\/(.*)/;
            my $group = $1;
            push @mapFiles1, "$prepDir1/$group/$group.fa.$ext{$param{mapType}}";
            push @mapFiles2, "$prepDir2/$group/$group.fa.$ext{$param{mapType}}";
        }
        my $pairsTable;
        if ($param{addToExisting}) {
            $pairsTable = getTableName('Pairs', $outputSample);
        } else {
            $pairsTable = newTable('Pairs', $outputSample);
        }     
        runExtraction(\@mapFiles1, \@mapFiles2, $pairsTable);

    unless($param{extractDataOnly}){
        status("creating pair size histogram...\n");
            my $histogramTable = newTable('h', $pairsTable);
            createPairSizeHistogram($histogramTable, $pairsTable);

        status("calculating pair statistics...\n");
            my $statsTable = newTable('Stats', $outputSample);    
            calculateNormalStats($histogramTable, $statsTable);
    }
    
    return $pairsTable;
    
}

sub extractPairs{
    my ($sample) = @_;
    my ($inputSample, $outputSample) = splitSampleName($sample);
     
    getDirectories($inputSample, \my%dirs);    
    my $paramFile = "$dirs{sample}/$outputSample\_extract\_parameters.txt";    
    printParameters($paramFile, 'extractPairs', $sample);    

    status("extracting pairs...\n");
    	getMapFiles($inputSample, \my%mapFiles);
    	my $pairsTable;
        if ($param{addToExisting}) {
            $pairsTable = getTableName('Pairs', $outputSample);
        } else {
            $pairsTable = newTable('Pairs', $outputSample);
        }  
        if ($param{mapType} eq 'bwa'){
            pairBWAMaps($sample);
            require "$param{vampPath}/bin/sam/extractPairedSam.pl";
            extractPairedSam($sample); 
        } else {
            runExtraction([$mapFiles{read1}], [$mapFiles{read2}], $pairsTable); 
        }

    unless($param{extractDataOnly}){
        status("creating pair size histogram...\n");
            my $histogramTable = newTable('h', $pairsTable);
            createPairSizeHistogram($histogramTable, $pairsTable);

        status("calculating pair statistics...\n");
            my $statsTable = newTable('Stats', $outputSample);    
            calculateNormalStats($histogramTable, $statsTable);
    }
    
    return $pairsTable;
    
}

my ($readsTable, %loadNs, %loadFileHs);
sub runExtraction{
    my ($inputFiles1, $inputFiles2, $pairsTable) = @_; #$inputFiles1 and 2 are array refs
    
    status("  parsing map data...\n");
        my $readsTable1 = getTableName('Reads', "$pairsTable\_1");
        my $readsTable2 = getTableName('Reads', "$pairsTable\_2");
        %loadNs = ();
        parseMapData($inputFiles1, $readsTable1, 1);
        parseMapData($inputFiles2, $readsTable2, 2);
        
    my $np = 0;
    my $outputFile = "$pairsTable.csv";
    my $timeStamp = getTimeStamp();
    
    foreach my $loadN(keys %loadNs){
        my $million = $loadN + 1;
        status("   million pairs: $million...\n");
    
        status("    loading map data...\n");
            my $readsTable1 = newTable('Reads', "$pairsTable\_1");
            my $readsTable2 = newTable('Reads', "$pairsTable\_2");
            loadData("$readsTable1\_$loadN\.csv", $readsTable1, ",", "PAIRID, MAPS CHAR(4000)");  
            loadData("$readsTable2\_$loadN\.csv", $readsTable2, ",", "PAIRID, MAPS CHAR(4000)");  
            open my $outputFileH, ">", $outputFile or die "could not open $outputFile";    
    
        status("    retrieving correlated sample maps...\n");
            runSQL("SELECT nvl(t1.PAIRID, t2.PAIRID) PAIRID, 
                           nvl(t1.MAPS, 0) MAPS1, nvl(t2.MAPS, 0) MAPS2
                   FROM $readsTable1 t1 FULL OUTER JOIN $readsTable2 t2
                   ON t1.PAIRID = t2.PAIRID",
    	           \my($pairID, $maps1, $maps2) );      
    
        status("    determining pair types...\n");
        MAP_: while (fetchRow()){
        	$np++;
        	$pairID = ($pairID * 1E8) + $timeStamp;
            if($maps1 xor $maps2){
                my $maps = $maps1;
                $maps or $maps = $maps2;
                my @maps = split(/::/, $maps);
                commitSingletonMaps($outputFileH, $pairID, \@maps);
            } elsif($maps1 and $maps2){
                my @maps1 = split(/::/, $maps1);
                my @maps2 = split(/::/, $maps2);
                commitPairedMaps($outputFileH, $pairID, \@maps1, \@maps2);
            }
        }  
        close $outputFileH;
        dropTable($readsTable1);
        dropTable($readsTable2);
        
        status("    uploading pairs and discrepancies...\n");
            loadData($outputFile, $pairsTable, ",", "PAIRTYPE, 
                        CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1, NHITS1,
    					CHROMOSOME2, POSITION2, LENGTH2, STRAND2, DISCREPANCIES2, NHITS2,
    					FRAGMENTSIZE, PAIRID, NFRAGS");

    }

	status("    $np two-read pairs processed\n");
    				
}

sub splitReadMap {
    my ($map) = @_;
    if(ref($map)){ #$map is a passed array reference (used by extractPairedBam)
        return @$map;
    } else { #$map is a joined scalar
        return split(/:/, $map);
    }
}

sub commitSingletonMaps {
    my ($outputFileH, $pairID, $maps1) = @_;
    my $nHits1 = scalar @$maps1;
    my $nFrags = $nHits1;
    my $pairType = $types{Pairs}{SingleRead};
    foreach my $map1 (@$maps1){
        my ($chrom1, $pos1, $length1, $strand1, $discs1) = splitReadMap($map1);
        my $size = 0;
        print $outputFileH join(",", ($pairType, $chrom1, $pos1, $length1, $strand1, $discs1, $nHits1,
                                	             0, 0, 0, 0, 0, 0,
                                		         $size, $pairID, $nFrags))."\n";     
    }
}

sub commitPairedMaps {
    my ($outputFileH, $pairID, $maps1, $maps2) = @_;
    my $nHits1 = scalar @$maps1;
    my $nHits2 = scalar @$maps2;
    my $nFrags = $nHits1 * $nHits2;
	foreach my $map1 (@$maps1){
	    my ($chrom1, $pos1, $length1, $strand1, $discs1) = splitReadMap($map1);
	    foreach my $map2 (@$maps2){
            my ($chrom2, $pos2, $length2, $strand2, $discs2) = splitReadMap($map2);    	    
	        my ($chrom1, $pos1, $length1, $strand1, $discs1) = ($chrom1, $pos1, $length1, $strand1, $discs1); #temp hold for read1 since may invert
            my ($pairType, $size);
            if ($param{rooted}){ #do NOT re-orient rooted pairs, so read1 is always the partner/unrooted read
        		if ($chrom1 == $chrom2){                
                    $size = abs($pos2 - $pos1) + $param{readLength};   
                    $chrom2 = 0; #override chrom2 to 0 as way of tracking sameChrom pairs
        		    if ($strand1 == $strand2){
                        $pairType = $types{Pairs}{Colinear}; 
        		    } else {
            			if ($strand1 == 1 and $pos2 > $pos1){
            			    $pairType = $types{Pairs}{Convergent};
            			} else {
            				$pairType = $types{Pairs}{Divergent};				    
            			}
        		    }
        		} else {                           
                    $size = 0;  #size not meaningful for DiffChrom
        		    $pairType = $types{Pairs}{DiffChrom};
        		}
            } else {
        		if ($chrom1 == $chrom2){
        		    #re-orient all sameChrom pairs so that smallest position is called position1
        		    if ($pos1 > $pos2){ invertReads(\($chrom1, $pos1, $length1, $strand1, $discs1, $nHits1,
                                                      $chrom2, $pos2, $length2, $strand2, $discs2, $nHits2)) }    
                    $size = $pos2 - $pos1 + $length2;
                    $chrom2 = 0; #override chrom2 to 0 as way of tracking sameChrom pairs
        		    if ($strand1 == $strand2){
                        $pairType = $types{Pairs}{Colinear};          
        		    } else {
            			if ($strand1 == 1){
            			    $pairType = $types{Pairs}{Convergent};
            			} else {
            				$pairType = $types{Pairs}{Divergent};				    
            			}

        		    }
        		} else {
        		    ($param{noDiffChrom} or $nFrags > $param{maxHits} + 1) and next; #suppress DiffChrom for highly promiscuous pairs
        		    #re-orient all DiffChrom read pairs so that lowest number chrom is called chrom1        		
        		    if ($chrom1 > $chrom2){ invertReads(\($chrom1, $pos1, $length1, $strand1, $discs1, $nHits1,
                                                          $chrom2, $pos2, $length2, $strand2, $discs2, $nHits2)) }                             
                    $size = 0;  #size not meaningful for DiffChrom
        		    $pairType = $types{Pairs}{DiffChrom};
        		}
                                     
            }
            print $outputFileH join(",", ($pairType, $chrom1, $pos1, $length1, $strand1, $discs1, $nHits1,
                                    	             $chrom2, $pos2, $length2, $strand2, $discs2, $nHits2,
                                    		         $size, $pairID, $nFrags))."\n";                                                                                                        			
	    }
	} 
}

sub parseMapData{  #data are split into chunks of pairIDs for efficiency
    my ($inputFiles, $readsTable_, $read) = @_;
    status("    parsing read$read maps\n");
    $readsTable = $readsTable_;
    %loadFileHs = ();
    foreach my $inputFile(@$inputFiles){
        my $getLineSub = "getLine\_$param{mapType}";
        $param{forceGetLineSub} and $getLineSub = $param{forceGetLineSub};
        my $getLineSubRef = \&{$getLineSub};
        open my $inputFileH, "<", $inputFile or die "could not open $inputFile";
        my ($pairID, $map, $prevPairID);
        my $line = <$inputFileH>;
        #($prevPairID, $prevMaps) = &$getLineSubRef(\$line);
        ($prevPairID, $map) = &$getLineSubRef(\$line);
        my $loadFileH = getLoadFileH($prevPairID);
        print $loadFileH "\n$prevPairID,$$map";
        while (<$inputFileH>){  
            ($pairID, $map) = &$getLineSubRef(\$_);
            if ($pairID){
                if ($pairID eq $prevPairID){
                    #$prevMaps .= "::$$map";
                    print $loadFileH "::$$map";
                } else {
                    $loadFileH = getLoadFileH($pairID);
                    print $loadFileH "\n$pairID,$$map";
                    $prevPairID = $pairID; 
                    #my $loadFileH = getLoadFileH($readsTable, $pairID);
                    #print $loadFileH "$prevPairID,$prevMaps\n";
                    #$prevMaps = $$map;  
                }             
            }  
        }
        close $inputFileH;
    }
    foreach my $loadN (keys %loadFileHs) {
        my $loadFileH = $loadFileHs{$loadN};
        close $loadFileH;
    }       
}

sub getLoadFileH{
    my ($pairID) = @_;
    my $loadN = int($pairID / 1E6);
    $loadFileHs{$loadN} and return $loadFileHs{$loadN};
    $loadNs{$loadN} = $loadN;
    my $loadFile = "$readsTable\_$loadN\.csv";
    open my $loadFileH, ">", $loadFile or die  "could not open $loadFile";
    $loadFileHs{$loadN} = $loadFileH;
    return $loadFileHs{$loadN};
}

sub getFields{
    my ($line) = @_;
    chomp $$line;
    my @fields = split(/\t/, $$line);
    return \@fields;
}

sub addDiscrepancy{
    my ($discs, $relPos, $type, $extra, $encounteredRef) = @_;
    my $disc = ($relPos * 1E2) + ($type * 1E1) + $extra;
    unless ($$encounteredRef{$disc}) {$discs = ($discs * 1E4) + $disc}
    $$encounteredRef{$disc}++;
    return $discs;
    #use of the encountered hash purges duplicate discrepancies within a read
    #necessary since pass seems to have a bug
    #evidenced by the apparent duplicate gap reported in this actual gff line:
    #   Note="M:0,G:3 -> Q10/C Q10/C R26/A/A";
    #or is the duplicate Q10/C somehow meaningful??
}

sub invertReads{
    my ($chrom1, $pos1, $length1, $strand1, $discs1, $nHits1,
        $chrom2, $pos2, $length2, $strand2, $discs2, $nHits2) = @_; #all as scalar Refs
    ($$chrom1, $$chrom2) = ($$chrom2, $$chrom1);
    ($$pos1, $$pos2) = ($$pos2, $$pos1);
    ($$length1, $$length2) = ($$length2, $$length1);
    ($$strand1, $$strand2) = ($$strand2, $$strand1);
    ($$discs1, $$discs2) = ($$discs2, $$discs1);
    ($$nHits1, $$nHits2) = ($$nHits2, $$nHits1); 
}

sub getTimeStamp{
    my ($sec, $min, $hr, $day, $month, $year) = localtime(time);
    return (($day * 1E6) +
	    ($hr * 1E4) +
	    ($min * 1E2) +
	    $sec);
}

sub createPairSizeHistogram{
    my ($histogramTable, $pairsTable) = @_;
    runSQL("INSERT INTO $histogramTable
            (SELECT PAIRTYPE SERIES, FRAGMENTSIZE X, Count(PAIRID) Y
                FROM $pairsTable
                GROUP BY PAIRTYPE, FRAGMENTSIZE)  ");
}

1;

