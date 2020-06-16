#!/usr/bin/perl
use strict;
use warnings;

##################################################
#ExtractUnpaired extracts and allows for the import
#of 'unpaired'/single-end read (i.e. not paired or mate-pair)
#into VAMP for analyes such as discrepancy finding
##################################################

use vars(qw(%param %types %fields %refSeqs));

my @binSizes = (100,1000);

sub extractUnpaired{
    my ($sample) = @_;
    my ($inputSample, $outputSample) = splitSampleName($sample);
    
	getMapFiles($inputSample, \my%mapFilesIn);            
	my $fragsTableOut;
	if ($param{addToExisting}) {
	    $fragsTableOut = getTableName('Frags', $outputSample); 	   
	} else {
	    $fragsTableOut = newTable('Frags', $outputSample); 
	}

    status("extracting mapped reads...\n");    
        loadMapsUnpaired(\%mapFilesIn, $fragsTableOut);

    unless($param{extractDataOnly}){  
        #this is ideal for unique mapped read counting, e.g. Chip-Seq...
        status("creating bin counts and histograms...\n");
            createUnpairedBinCounts($fragsTableOut, $outputSample);
                
        status("creating read map...\n"); 
            my $readMapTable = newTable('RMap', $outputSample);    
    #        createReadMapUnpaired($fragsTableOut, $readMapTable);      
    #        my $histogramTable = newTable('h', $readMapTable);
    #        createMapHistogram($histogramTable, $readMapTable);  

            createReadMapUnpairedAllBins($fragsTableOut, $readMapTable);
    }
    
    return $fragsTableOut;
        			
}    

sub loadMapsUnpaired{
    my ($mapFilesRef, $fragsTable) = @_;
	my $fragsFile = "$fragsTable.csv";           
	open my $fragsFileH, ">", $fragsFile or die "could not open $fragsFile"; 
    my $timeStamp = getTimeStamp(); 
    my $getLineSub = "getLine\_$param{mapType}";
    $param{forceGetLineSub} and $getLineSub = $param{forceGetLineSub};
    my $getLineSubRef = \&{$getLineSub};
    open my $mapFileH, "<", $$mapFilesRef{read1} or return;
    my ($pairID, $map, $prevPairID, @prevMaps);
    my $line = <$mapFileH>; 
    ($prevPairID, $map) = &$getLineSubRef(\$line);   
    push @prevMaps, $map;
    while (<$mapFileH>){     
        ($pairID, $map) = &$getLineSubRef(\$_);
        if ($pairID) {
            if ($pairID eq $prevPairID){
                push @prevMaps, $map;
            } else {    
                printReadsUnpaired($fragsFileH, $timeStamp, $prevPairID, \@prevMaps);              
                @prevMaps = ($map);
                $prevPairID = $pairID; 
            } 
        }    
    }
    printReadsUnpaired($fragsFileH, $timeStamp, $prevPairID, \@prevMaps);
    close $mapFileH;     
	close $fragsFileH; 
    status("  uploading read data...\n");
    loadData($fragsFile, $fragsTable, ",", "FRAGMENTID, FRAGMENTTYPE, 
                                                  CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1, NHITS1,
					                              PAIRID, NFRAGS");	      
}

sub printReadsUnpaired{
    my ($fragsFileH, $timeStamp, $pairID, $mapsRef) = @_;
    $pairID = ($pairID * 1E8) + $timeStamp; 
    my $nReads = scalar @$mapsRef;   
    foreach my $readMap(@$mapsRef){  
        my ($chrom, $pos, $length, $strand, $discs) = split(/:/, $$readMap); 
        print $fragsFileH join(",", $pairID, $types{Frags}{SingleRead},
                                    $chrom, $pos, $length, $strand, $discs, $nReads,
                                    $pairID, $nReads )."\n";         
    }  
}
  
sub createReadMapUnpaired{
    my ($fragsTable, $readMapTable) = @_;
    foreach my $chrom (1..nChrom()){
        my %coverage = ();
        runSQL("SELECT ((trunc(POSITION1/$param{binSize})*$param{binSize})+$param{binSize}) AS LOWBIN1,
                        (trunc((POSITION1 + LENGTH1 - 1)/$param{binSize})*$param{binSize}) AS HIGHBIN1
                FROM $fragsTable
                WHERE CHROMOSOME1 = $chrom",
                \my($lowBin1, $highBin1));
        while (fetchRow()){ fillRead(\%coverage, $lowBin1, $highBin1) }
        my $readMapFile = "$readMapTable.csv";
        open my $readMapFileH, ">", $readMapFile;
        foreach my $pos (keys %coverage){ print $readMapFileH join(",", ($chrom, $pos, $coverage{$pos}))."\n" }
        close $readMapFileH;
        loadData($readMapFile, $readMapTable, ",", "CHROMOSOME, POSITION, COVERAGE");
    }
}

sub createUnpairedBinCounts{
    my ($fragsTable, $sample) = @_;
    runSQL("SELECT Sum(END_) GENOMESIZE
            FROM CHROMINFO_$param{refSeqBase}
            WHERE CHROMOSOME <= $refSeqs{$param{refSeqBase}}{nChrom}",
            \my $genomeSize);
    fetchRow();  
    status("  refSeq $param{refSeqBase} genome size = $genomeSize\n");
    runSQL("SELECT Count(*) 
            FROM $fragsTable
            WHERE CHROMOSOME1 <= $refSeqs{$param{refSeqBase}}{nChrom}",
            \my$nReads);
    fetchRow();  
    status("  sample $sample number of reads = $nReads\n");      
    foreach my $binSize(@binSizes){
        my $expectedCoverage = ($binSize / $genomeSize) * $nReads;
        status("    binSize $binSize yields expected random bin coverage = $expectedCoverage\n");
        my $readMapTable = getTableName('RMap', "$sample\_$binSize");
        my $sql = "SELECT CHROMOSOME1 CHROMOSOME, (Round(POSITION1 / $binSize) * $binSize) POSITION
                   FROM $fragsTable
                   WHERE CHROMOSOME1 <= $refSeqs{$param{refSeqBase}}{nChrom}";
        $sql = "SELECT CHROMOSOME, POSITION, 
                       Count(*) COVERAGE, 
                       Round((Count(*) / $expectedCoverage)/0.1)*0.1 NORMALIZEDCOVERAGE
                FROM ($sql)
                GROUP BY CHROMOSOME, POSITION";
        dropTable($readMapTable);
        runSQL("CREATE TABLE $readMapTable AS $sql");
        my $histogramTable = getTableName('h', $readMapTable);
        $sql = "SELECT 1 SERIES, NORMALIZEDCOVERAGE X, Count(*) Y
                FROM $readMapTable
                GROUP BY NORMALIZEDCOVERAGE";
        dropTable($histogramTable);
        runSQL("CREATE TABLE $histogramTable AS $sql");
    }
}

sub createReadMapUnpairedAllBins{
    my ($fragsTable, $readMapTable) = @_;
    runSQL("SELECT Sum(END_) GENOMESIZE
            FROM CHROMINFO_$param{refSeqBase}
            WHERE CHROMOSOME <= $refSeqs{$param{refSeqBase}}{nChrom}",
            \my $genomeSize);
    fetchRow();  
    status("  refSeq $param{refSeqBase} genome size = $genomeSize\n");
    runSQL("SELECT Count(*) 
            FROM $fragsTable
            WHERE CHROMOSOME1 <= $refSeqs{$param{refSeqBase}}{nChrom}",
            \my$nReads);
    fetchRow();  
    status("  $fragsTable number of reads = $nReads\n"); 
    my $expectedCoverage = ($param{binSize} / $genomeSize) * $nReads;
    status("  binSize $param{binSize} yields expected random bin coverage = $expectedCoverage\n");
    foreach my $chrom (1..nChrom()){
        runSQL("SELECT (Round(END_ / $param{binSize}) * $param{binSize}) AS HIGHBIN
                FROM CHROMINFO_$param{refSeqBase}
                WHERE CHROMOSOME = $chrom",
                \my $highBin);
        fetchRow();
        my %coverage = ();
        for (my $bin = 0; $bin <= $highBin; $bin += $param{binSize}){ $coverage{$bin} = 0 }
        runSQL("SELECT (Round(POSITION1 / $param{binSize}) * $param{binSize}) BIN
                   FROM $fragsTable
                   WHERE CHROMOSOME1 = $chrom",
                \my $bin);
        while (fetchRow()){ $coverage{$bin}++ }
        my $readMapFile = "$readMapTable.csv";
        open my $readMapFileH, ">", $readMapFile;
        foreach my $pos (keys %coverage){ 
            print $readMapFileH join(",", ($chrom, $pos, $coverage{$pos}, $coverage{$pos}/$expectedCoverage))."\n" 
        }
        close $readMapFileH;
        loadData($readMapFile, $readMapTable, ",", "CHROMOSOME, POSITION, COVERAGE, NORMALIZEDCOVERAGE");
    }
}

1;


