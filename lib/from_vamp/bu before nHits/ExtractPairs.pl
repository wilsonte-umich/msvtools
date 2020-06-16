#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#ExtractPairs.pl coordinates the conversion of text files of mapped reads
#into Oracle Pairs data tables, and then calculates Pair statistics.
#########################################################################

use vars(qw(%param %types %fields %refSeqs));

sub extractPairs{
    my ($sample) = @_;
    my ($inputSample, $outputSample) = split("::", $sample);    
    $outputSample or $outputSample = $inputSample;   
     
    getDirectories($inputSample, \my%dirs);    
    my $paramFile = "$dirs{sample}/$outputSample\_extract\_parameters.txt";    
    printParameters($paramFile, 'extractPairs', $sample);    

    status("extracting pairs...\n");
    	getMapFiles($inputSample, \my%mapFiles);
    	my $pairsTable = newTable('Pairs', $outputSample);
      	runExtraction($mapFiles{read1}, $mapFiles{read2}, $pairsTable); 

    status("creating pair size histogram...\n");
        my $histogramTable = newTable('h', $pairsTable);
        createPairSizeHistogram($histogramTable, $pairsTable);

    status("calculating pair statistics...\n");
    	my $statsTable = newTable('Stats', $outputSample);    
        calculateNormalStats($histogramTable, $statsTable);
}

sub runExtraction{
    my ($inputFile1, $inputFile2, $pairsTable) = @_;

    status("  loading map data...\n");
        my $table1 = loadMaps($inputFile1, $pairsTable, 1);
        my $table2 = loadMaps($inputFile2, $pairsTable, 2);

        my $outputFile = "$pairsTable.csv";
        open my $outputFileH, ">", $outputFile or die "could not open $outputFile";    

#    status("  retrieving correlated sample maps...\n");
#        runSQL("SELECT $table1.PAIRID AS PAIRID, $table1.MAPS AS MAPS1, $table2.MAPS AS MAPS2
#               FROM $table1, $table2
#               WHERE $table1.PAIRID = $table2.PAIRID",
#	           \my($pairID, $maps1, $maps2) );
    status("  retrieving correlated sample maps...\n");
        runSQL("SELECT nvl(t1.PAIRID, t2.PAIRID) PAIRID, 
                       nvl(t1.MAPS, 0) MAPS1, nvl(t2.MAPS, 0) MAPS2
               FROM $table1 t1 FULL OUTER JOIN $table2 t2
               ON t1.PAIRID = t2.PAIRID",
	           \my($pairID, $maps1, $maps2) );      

    status("  determining pair types...\n");
        my $timeStamp = getTimeStamp();
    
    my $np = 0;
    MAP_: while (fetchRow()){
    	$np++;
    	$pairID = ($pairID * 1E8) + $timeStamp;
    	if(!($maps1 and $maps2)){ #the SingleReads block
    	   my $maps = $maps1;
    	   $maps or $maps = $maps2;
    	   my @maps = split(/::/, $maps);
    	   my $nFrags = scalar @maps;
    	   my $pairType = $types{Pairs}{SingleRead};
           foreach my $map (@maps){
                my ($chrom1, $pos1, $length1, $strand1, $discs1) = split(/:/, $map);
                my $size = 0;
                print $outputFileH join(",", ($pairType, $chrom1, $pos1, $length1, $strand1, $discs1,
                                        	             0, 0, 0, 0, 0,
                                        		         $size, $pairID, $nFrags))."\n";     
           }
           next MAP_;
        }
    	my @maps1 = split(/::/, $maps1);
    	my @maps2 = split(/::/, $maps2);
    	my $nFrags = @maps1 * @maps2;
    	foreach my $map1 (@maps1){
    	    my ($chrom1, $pos1, $length1, $strand1, $discs1) = split(/:/, $map1);
    	    foreach my $map2 (@maps2){
                my ($chrom2, $pos2, $length2, $strand2, $discs2) = split(/:/, $map2);    	    
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
            		    if ($pos1 > $pos2){ invertReads(\($chrom1, $pos1, $length1, $strand1, $discs1,
                                                          $chrom2, $pos2, $length2, $strand2, $discs2)) }    
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
            		    ($param{noDiffChrom} or $nFrags > $param{maxHits}) and next; #suppress DiffChrom for highly promiscuous pairs
            		    #re-orient all DiffChrom read pairs so that lowest number chrom is called chrom1        		
            		    if ($chrom1 > $chrom2){ invertReads(\($chrom1, $pos1, $length1, $strand1, $discs1,
                                                              $chrom2, $pos2, $length2, $strand2, $discs2)) }                             
                        $size = 0;  #size not meaningful for DiffChrom
            		    $pairType = $types{Pairs}{DiffChrom};
            		}
                                         
                }
                print $outputFileH join(",", ($pairType, $chrom1, $pos1, $length1, $strand1, $discs1,
                                        	             $chrom2, $pos2, $length2, $strand2, $discs2,
                                        		         $size, $pairID, $nFrags))."\n";                                                                                                        			
    	    }
    	}
    }  
    close $outputFileH;
    dropTable($table1);
    dropTable($table2);
    
    status("  uploading final pairs and discrepancies...\n");
        loadData($outputFile, $pairsTable, ",", "PAIRTYPE, 
                    CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1,
					CHROMOSOME2, POSITION2, LENGTH2, STRAND2, DISCREPANCIES2,
					FRAGMENTSIZE, PAIRID, NFRAGS");
					
	status("    $np two-read pairs processed\n");				
}

sub loadMaps{
    my ($inputFile, $pairsTable, $read) = @_;
    status("    loading read $read\n");
        my $getLineSub = "getLine\_$param{mapType}";
        my $getLineSubRef = \&{$getLineSub};
        my $table = newTable('Reads', "$pairsTable$read");
        my $outputFile = "$table\.csv";
        open my $outputFileH, ">", $outputFile or die "could not open $outputFile";
        open my $inputFileH, "<", $inputFile or die "could not open $inputFile";
        my ($pairID, $map, $prevPairID, $prevMaps);
        my $line = <$inputFileH>;
        ($prevPairID, $prevMaps) = &$getLineSubRef($line);
        while (<$inputFileH>){  
            ($pairID, $map) = &$getLineSubRef($_);
            if ($pairID){
                if ($pairID eq $prevPairID){
                    $prevMaps .= "::$map";
                } else {
                    print $outputFileH "$prevPairID,$prevMaps\n";
                    $prevMaps = $map;  
                }
                $prevPairID = $pairID;             
            }  
        }
        print $outputFileH "$prevPairID,$prevMaps\n";
        close $outputFileH;
        close $inputFileH;
        loadData($outputFile, $table, ",", "PAIRID, MAPS CHAR(4000)"); 
        return $table;
}

sub getFields{
    my ($line) = @_;
    chomp $line;
    my @fields = split(/\t/,$line);
    return @fields;
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
    my ($chrom1, $pos1, $length1, $strand1, $discs1,
        $chrom2, $pos2, $length2, $strand2, $discs2) = @_; #all as scalar Refs
    ($$chrom1, $$chrom2) = ($$chrom2, $$chrom1);
    ($$pos1, $$pos2) = ($$pos2, $$pos1);
    ($$length1, $$length2) = ($$length2, $$length1);
    ($$strand1, $$strand2) = ($$strand2, $$strand1);
    ($$discs1, $$discs2) = ($$discs2, $$discs1); 
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



