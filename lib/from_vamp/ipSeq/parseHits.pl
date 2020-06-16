#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $param{stranded} or $param{stranded} = 0; #whether or not the IP hits are strand-specific
defined $param{reverseStrands} or $param{reverseStrands} = 0; #whether or not to reverse strand assignments to make top strand = 1 (determined empirically)
defined $param{binSizes} or $param{binSizes} = 0; #comma-delimited list, uses binSize if not specified
#defined $param{fragmentLength} or $param{fragmentLength} = 200; #the estimated average size of library fragments
#defined $param{fMapBinSize} or $param{fMapBinSize} = 50; #binSize applied to the fragment map
defined $command{parseHits} or $command{parseHits} = ['multiThread', '4:00:00', 2000, 0];

#defined $command{createHitFMap} or $command{createHitFMap} = ['singleThread', '4:00:00', 2000, 0];

#my ($halfBinSize, $halfFragLength);

sub parseHits { 
    my ($sample) = @_;
    my $fragsTable = getTableName('Frags', $sample);
    my $hitsTable = newTable('Hits', $sample);
    runHitParse($fragsTable, $hitsTable);
    #known possible limitation:
    #  HMAP does not carry strand information at present (although HITS table does)
    #  anything derivative from HMAP will NOT have strand-specific bin counts
    createHitMapsAllBins($hitsTable, $sample);
    calcHitStats($sample);
}
sub runHitParse { 
    my ($fragsTable, $hitsTable) = @_;
    status("parsing sequence reads to genome hits...\n");
    my $strand = "strand1";
    $param{reverseStrands} and $strand = "decode(strand1,1,2,1)";
    $param{stranded} or $strand = "1";
    my $strandSQL = "SELECT chromosome1 chromosome, position1 position, $strand strand
                     FROM $fragsTable ";           
    my $countSQL = "SELECT chromosome, position, strand, count(*) count_
                    FROM ($strandSQL)
                    GROUP BY chromosome, position, strand";              
    updateTable('Hits', $hitsTable, $countSQL);   
}
sub createHitMapsAllBins{
    my ($hitsTable, $sample) = @_;
    status("creating hit coverage map...\n");
    my $genomeSize = getGenomeSize();
    my $nHits = getGenomeHitCount($hitsTable);
    my $expectedDensity = $nHits / $genomeSize;
    status("  average random hit density = $expectedDensity\n");
    $param{binSizes} or $param{binSizes} = $param{binSize};
    my @binSizes = split(",", $param{binSizes});
    foreach my $binSize(@binSizes){  
        my $hMapTable = newTable('HMap', "$sample\_$binSize");  
        createHitMap($hitsTable, $hMapTable, $binSize, $expectedDensity); 
        my $histogramTable = newTable('h', $hMapTable);
        createHitMapHistogram($histogramTable, $hMapTable);
    }
}
sub createHitMap{
    my ($hitsTable, $hMapTable, $binSize, $expectedDensity) = @_;
    my $expectedCoverage = $expectedDensity * $binSize;
    status("  binSize $binSize yields random bin coverage = $expectedCoverage\n");
    foreach my $chrom (1..nChrom()){
        my $highBin = getHighBin($chrom, $binSize);
        my %coverage = ();
        for (my $bin = 0; $bin <= $highBin; $bin += $binSize){ $coverage{$bin} = 0 }
        my $val = "1";
        $param{keepDups} and $val = "count_";
        my $binSQL = "SELECT (Round(POSITION / $binSize) * $param{binSize}) bin, $val count_
                      FROM $hitsTable
                      WHERE chromosome = $chrom";
        runSQL($binSQL, \my($bin, $count));
        while (fetchRow()){ $coverage{$bin} += $count }
        my $hMapFile = "$hMapTable.csv";
        open my $hMapFileH, ">", $hMapFile;
        foreach my $bin (keys %coverage){ 
            my $coverage = $coverage{$bin}; #number of hits in the bin
            my $density = $coverage / $binSize; #number of hits per bp in the bin
            my $normCoverage = $coverage / $expectedCoverage; 
            my $normDensity = $density / $expectedDensity; #norm values are actually the same...
            print $hMapFileH join(",", ($chrom, $bin, $coverage, $density, $normCoverage, $normDensity))."\n" 
        }
        close $hMapFileH;
        loadData($hMapFile, $hMapTable, ",", "CHROMOSOME, POSITION, COVERAGE, DENSITY, NORMALIZEDCOVERAGE, NORMALIZEDDENSITY");
    }
}


#sub createHitFMap{
#    my ($fragsTable, $sample) = @_;
#    my $fMapTable = newTable('FMap', "hits_$sample");    
#    $param{fMapBinSize} % 2 and $param{fMapBinSize}++;
#    $halfBinSize = $param{fMapBinSize} / 2; 
#    $param{fragmentLength} % 2 and $param{fragmentLength}++;
#    $halfFragLength = $param{fragmentLength} / 2; 
#
#    #foreach my $chrom (1..nChrom()){
#    foreach my $chrom (21..21){
#    
#        my %coverage = ();    
#        #my $highBin = getHighBin($chrom, $param{fMapBinSize});
#        #for (my $bin = 0; $bin <= $highBin; $bin += $param{fMapBinSize}){ $coverage{$bin} = 0 }
#        my $hardSql =  "SELECT strand1,
#                               decode(strand1, 1, position1, position1 + length1 - $param{fragmentLength}) hardStart,  
#                               decode(strand1, 1, position1 + $param{fragmentLength} - 1, position1 + length1 - 1) hardEnd
#                         FROM $fragsTable
#                         WHERE CHROMOSOME1 = $chrom
#                         
#                           AND position1 >= 22300000
#                           AND position1 <= 22950000
#                         
#                         ";    
#        my $fuzzySql =  "SELECT strand1, hardStart, hardEnd,
#                               decode(strand1, 1, hardEnd + 1, hardStart - $halfFragLength) fuzzyStart,  
#                               decode(strand1, 1, hardEnd + $halfFragLength, hardStart - 1) fuzzyEnd           
#                         FROM ($hardSql)";                            
#        runSQL($fuzzySql, \my($strand, $hardStart, $hardEnd, $fuzzyStart, $fuzzyEnd));               
#        while (fetchRow()){ 
#            my $lowBin =  getHitFMapBin($hardStart);
#            my $highBin = getHitFMapBin($hardEnd);
#            for (my $bin = $lowBin; $bin <= $highBin; $bin += $param{fMapBinSize}){ 
#                my ($countStart, $countEnd) = getHitFMapCountEnds($bin, $hardStart, $hardEnd);  
#                my $scaledCount = ($countEnd - $countStart + 1) / $param{fMapBinSize}; #can be fractional
#                $coverage{$bin} += $scaledCount;
#            }
#            $lowBin =  getHitFMapBin($fuzzyStart);
#            $highBin = getHitFMapBin($fuzzyEnd);
#            for (my $bin = $lowBin; $bin <= $highBin; $bin += $param{fMapBinSize}){ 
#                my ($countStart, $countEnd) = getHitFMapCountEnds($bin, $fuzzyStart, $fuzzyEnd);
#                my $startIntensity = ($countStart - $fuzzyStart) / $halfFragLength;
#                my $endIntensity = ($fuzzyEnd - $countEnd) / $halfFragLength;
#                $strand == 1 and $startIntensity = 1 - $startIntensity;
#                $strand == 2 and $endIntensity = 1 - $endIntensity;
#                my $averageIntensity = ($startIntensity + $endIntensity) / 2;   
#                my $scaledCount = $averageIntensity * ($countEnd - $countStart + 1) / $param{fMapBinSize}; #can be fractional
#                $coverage{$bin} += $scaledCount;  
#            }
#        }
#        my $fMapFile = "$fMapTable.csv";
#        open my $fMapFileH, ">", $fMapFile;
#        foreach my $bin (keys %coverage){ 
#            my $coverage = int($coverage{$bin} + 0.5); 
#            print $fMapFileH join(",", ($chrom, $bin, $coverage))."\n" 
#        }
#        close $fMapFileH;
#        loadData($fMapFile, $fMapTable, ",", "CHROMOSOME, POSITION, COVERAGE");
#    }   
#}
#sub getHitFMapBin { 
#    my ($position) = @_;
#    return int(($position/$param{fMapBinSize}) + 0.5) * $param{fMapBinSize};
#}
#sub getHitFMapCountEnds { 
#    my ($bin, $fragStart, $fragEnd) = @_;
#    my $binStart = $bin - $halfBinSize;
#    my $binEnd =   $bin + $halfBinSize - 1;
#    my $countStart = $fragStart >= $binStart ? $fragStart : $binStart;
#    my $countEnd = $fragEnd <= $binEnd ? $fragEnd : $binEnd;
#    return ($countStart, $countEnd);
#}



sub getGenomeSize {
    runSQL("SELECT Sum(END_) GENOMESIZE
            FROM CHROMINFO_$param{refSeqBase}
            WHERE CHROMOSOME <= $refSeqs{$param{refSeqBase}}{nChrom}",
            \my $genomeSize);
    fetchRow();  
    status("  refSeq $param{refSeqBase} genome size = $genomeSize\n"); 
    return $genomeSize;
}
sub getGenomeHitCount {
    my ($hitTable) = @_;
    my $agg = "count(*)";
    $param{keepDups} and $agg = "sum(count_)";
    runSQL("SELECT $agg FROM $hitTable",
            \my $nHits);
    fetchRow();  
    status("  $hitTable number of hits = $nHits\n"); 
    return $nHits;
}
sub getHighBin {
    my ($chrom, $binSize) = @_;
    runSQL("SELECT (Round(END_ / $binSize) * $binSize) AS HIGHBIN
            FROM CHROMINFO_$param{refSeqBase}
            WHERE CHROMOSOME = $chrom",
            \my $highBin);
    fetchRow(); 
    return $highBin;
}
sub createHitMapHistogram {
    my ($histogramTable, $mapTable) = @_;
    my $roundSQL = "SELECT Round(NORMALIZEDDENSITY,2) NORMALIZEDDENSITY FROM $mapTable";
    runSQL("INSERT INTO $histogramTable
            (SELECT 0 SERIES, NORMALIZEDDENSITY X, Count(*) Y
                FROM ($roundSQL)
                GROUP BY NORMALIZEDDENSITY)  ");
}

1;
