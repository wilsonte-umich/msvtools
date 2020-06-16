#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %fieldNames));

###############################################################################################
#establishes the fraction of all hits within a bin that correspond to UV hits vs. UV + non-UV hits
#fractions over minFracUV threshold identify transcription start sites, iff sufficient UV hits are present
#----------------------------------------------------------------------------------------------
#although HITS table does carry strand information, TSS are determined WITHOUT respect to strand
#probably more sensitive for true blocks, can work out strand content later
###############################################################################################

#callable parameters and commands in this script
defined $param{minFracUV} or $param{minFracUV} = 0.65; #minimum fraction of UV hits for a bin to be called as part of a TSS
defined $param{minUVND} or $param{minUVND} = 1; #minimum UV normalizedDensity for calling a SINGLE bin as a TSS
defined $command{findTSS} or $command{findTSS} = ['singleThread', '12:00:00', 2000, 0];

my ($uvExpectedDensity, $tssFileH, $tssID, $nCountSum, $uCountSum, $nBins, $blockStart, $prevBin);
my $halfBinSize = int($param{binSize}/2);
my ($binSize, $fragLength, $halfFragLength, $effFragLength, $maNBins, $maNAdj, %sampleInfo, %coverage, $uvProb);
my $mapSourceIDSql = "fragmentID-trunc(fragmentID/1E8)*1E8";

sub findTSS { 
    my ($nonUVsample, $uvSample) = @_;
    status("finding transcription start sites...\n");
    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $uvHitCount) = getTSSScalars($nonUVsample, $uvSample);
    status("  nonUV scalar = $nonUVScalar\n  UV scalar = $uvScalar\n");
    my $genomeSize = getGenomeSize();
    $uvExpectedDensity = ($uvHitCount*$uvScalar)/$genomeSize;
    my $tssTable = newTable("HB", "$uvSample\_tss"); #output is a hit block table, modified to include the fracUV information
    runSQL("ALTER TABLE $tssTable ADD (FRACUVTHRESHOLD NUMBER, FRACUV NUMBER(*,5))"); 
    my $tssFile = "$tssTable.csv";
    open $tssFileH, ">", $tssFile or die "could not open $tssFile: $!\n";
    my ($nAgg, $uAgg) = ("decode(nCount,null,0,1)", "decode(uCount,null,0,1)");
    $param{keepDups} and ($nAgg, $uAgg) = ("nCount", "uCount");
    foreach my $chrom(1..nChrom()){
        status(" $chrom");
        #this is where strand tracking would have to be added, if desired; currently grouping by position
        my $nonUVHitsSQL = "SELECT position, sum(count_) count_ FROM $nonUVHitsTable WHERE chromosome = $chrom GROUP BY position";
        my $uvHitsSQL =    "SELECT position, sum(count_) count_ FROM $uvHitsTable    WHERE chromosome = $chrom GROUP BY position";
        my $hitJoinSQL = "SELECT round(nvl(n.position,u.position) / $param{binSize}) * $param{binSize} bin, 
                                 n.count_ * $nonUVScalar nCount, 
                                 u.count_ * $uvScalar uCount
                          FROM ($nonUVHitsSQL) n FULL OUTER JOIN ($uvHitsSQL) u
                            ON (n.position = u.position)";
        my $joinCountSQL = "SELECT bin, nvl(sum($nAgg),0) nCount, nvl(sum($uAgg),0) uCount
                            FROM ($hitJoinSQL)
                            GROUP BY bin
                            ORDER BY bin";            
        runSQL($joinCountSQL, \my($bin,$nCount,$uCount));
        ($nCountSum, $uCountSum, $nBins, $blockStart, $prevBin) = (0,0,0,0,0);    
        while(fetchRow()){ 
            my $fracUV = $uCount/($uCount + $nCount);
            if ($fracUV >= $param{minFracUV}){
                ($blockStart and $bin != $prevBin + $param{binSize}) and commitTssBlock($chrom);        
                $blockStart or $blockStart = $bin;
                $nCountSum += $nCount;
                $uCountSum += $uCount;
                $nBins++;
            } else {
                $blockStart and commitTssBlock($chrom);
            }  
            $prevBin = $bin;
        }
        $blockStart and commitTssBlock($chrom);
    }    
    status("\n");  
    close $tssFileH;
    my $fieldNames = "HBID, CHROMOSOME, START_, END_, BINSIZE, NBINS, THRESHOLD, NORMALIZEDDENSITY, FRACUVTHRESHOLD, FRACUV, INGENE";
    loadData($tssFile, $tssTable, ",", $fieldNames);
    scoreHitBlockGenes($tssTable);
}
sub commitTssBlock {
    my ($chrom) = @_;
    #note that wider TSS including more bins have a lower effective ND threshold
    #this approach in meangingful ONLY because regions were raised as candidates by an independent method
    my $uvDensity = $uCountSum / ($nBins * $param{binSize});
    my $uvNormDens = $uvDensity / $uvExpectedDensity;
    if ($uvNormDens * $nBins >= $param{minUVND}){ 
        my $fracUV = $uCountSum/($uCountSum + $nCountSum);
        $tssID++;
        $blockStart -= $halfBinSize;
        my $blockEnd = $prevBin + $halfBinSize;
        print $tssFileH "$tssID,$chrom,$blockStart,$blockEnd,$param{binSize},$nBins,$param{minUVND},$uvNormDens,$param{minFracUV},$fracUV,0\n";  
    }
    ($nCountSum, $uCountSum, $nBins, $blockStart) = (0,0,0,0);    
}
sub getTSSScalars {
    my ($nonUVsample, $uvSample) = @_;
    my $nonUVHitsTable = getTableName("Hits", $nonUVsample);
    my $uvHitsTable = getTableName("Hits", $uvSample);
    my $nonUVHitCount = getGenomeHitCount($nonUVHitsTable);
    my $uvHitCount = getGenomeHitCount($uvHitsTable);
    my $scalarNum = $nonUVHitCount;
    $scalarNum < $uvHitCount and $scalarNum = $uvHitCount;
    my $nonUVScalar = ($scalarNum/$nonUVHitCount);
    my $uvScalar = ($scalarNum/$uvHitCount);
    return ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $uvHitCount);   
}

1;

