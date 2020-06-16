#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %fieldNames));

###############################################################################################
#find UV enhancement peaks in the genome from paired UV and non-UV data sets
#----------------------------------------------------------------------------------------------
#method:
#  1) scale the samples so that each sample's genome hit count = the average hit count of the combined samples
#       this ensures that counts are not artificially inflated for Poisson calculation
#  2) project the mapped read position over an estimated RNA fragment length
#       projection includes a "hard" call = fragmentLenth, where intensity = 1,
#       as well as a "fuzzy" call = halfFragmentLenth at distal end only, 
#       over which intensity scales from 1 to 0
#  3) fractionally distribute the combined effective fragment length into bins crossed by the fragment
#       bins can and will be filled with fractional portions of the fragment,
#       the total in all bins for the fragment = 1, to maintain accurate numbers for Poisson
#  4) apply a moving average to smooth the coverage curve
#  5) find blocks of contiguous genome bins for which each bin shows UV > non-UV
#       these entire blocks will be kept as UVE peaks if they pass the next step
#       blocks are determined independently for the two strands if stranded=TRUE
#  6) check every block for statistically significant UV enrichment
#       statistical test is the cumulative Poisson,
#       comparing rounded UV count to lambda value = average(UV,non-UV)
#       this is synonymous with the null hypothesis that UV = non-UV, i.e. fracUV = 0.5
#       blocks checked on total span AND at bin with highest UV count; EITHER can pass the block
#  7) keep the entire block as a UVE peak when p <= pThreshold
#  8) merge adjacent blocks that are unnaturally split by zero coverage bins or rare bins where UV <= non-UV
#       merge when: p(summed adjacent blocks) > p(summed adjacent blocks plus gap)
#       apply maxGapLength to minimize improper merging
###############################################################################################

#callable parameters and commands in this script
defined $param{fragmentLength} or $param{fragmentLength} = 500; #the estimated average size of library fragments
defined $param{maNBins} or $param{maNBins} = 3; #the number of bins in the moving average sliding window
defined $param{pThreshold} or $param{pThreshold} = 0.05; #the number of bins in the moving average sliding window
defined $param{maxGapLength} or $param{maxGapLength} = 10000; #the number of bins in the moving average sliding window
defined $command{findUVE} or $command{findUVE} = ['singleThread', '24:00:00', 5000, 0];

my ($uvExpectedDensity, $uveFileH, $uveID);
my ($nonUVsample, $uvSample, %sampleInfo, %coverage, $uvProb, @bins, $minBin, $maxBin, $minI, $maxI, $newMinI) = @_;
my ($binSize, $halfBinSize, $fragLength, $halfFragLength, $effFragLength, $maNBins, $maNAdj, $maOffset, @blocks);
my ($nCountSum, $uCountSum, $maxStrand, $maxUCount, $bestCounts);
my $mapSourceIDSql = "fragmentID-trunc(fragmentID/1E8)*1E8";

sub findUVE { 
    ($nonUVsample, $uvSample) = @_;
    scaleSamplesUve();
    my $genomeSize = getGenomeSize();
    adjustInputParametersUve();
    $uvExpectedDensity = ($sampleInfo{$uvSample}{hitCount} * $sampleInfo{$uvSample}{scalar}) / $genomeSize; 
    status("uvExpectedDensity = $uvExpectedDensity\n");
    $uvProb = "$param{vampPath}\/bin\/bruSeq\/uvProb";
    $maxStrand = $param{stranded} ? 2 : 1;
    my $uveTable = newTable("HB", "$uvSample\_uve"); #output is a modified hit block table
    runSQL("ALTER TABLE $uveTable ADD (STRAND NUMBER, PVALUE NUMBER, FRACUV NUMBER(*,5))"); 
    my $uveFile = "$uveTable.csv";
    open $uveFileH, ">", $uveFile or die "could not open $uveFile: $!\n";    
    status("finding UV enhancement peaks...\n");
    foreach my $chrom (1..nChrom()){
    #foreach my $chrom (19..19){
        status("chrom $chrom\n");
        %coverage = ();
        my $sql = "SELECT end_ FROM chrominfo_$param{refSeq} WHERE chromosome = $chrom";
        runSQL($sql, \my($end));
        fetchRow();
        ($minBin, $maxBin) = (0, getTMapBin($end));  
        ($minI, $maxI) = (0, $maxBin / $binSize);
        @bins = map {$_*$binSize} $minI..$maxI;
        extractTMapCoverage($nonUVsample, $chrom); #read strands determine the fragment projetion in genome
        extractTMapCoverage($uvSample, $chrom);
        status("  finding blocks of adjacent bins where UV > non-UV\n");
        foreach my $orientation(1..$maxStrand){ #call blocks separately on each strand for stranded libraries
            @blocks = ();
            my $i_; #$x_ == xPrev
            ($nCountSum, $uCountSum, $maxUCount, $i_) = (0, 0, 0, undef);
            foreach my $i ($minI..$maxI){ 
                my $uCount = $coverage{$orientation}{$bins[$i]}{$uvSample};
                my $nCount = $coverage{$orientation}{$bins[$i]}{$nonUVsample};
                if ($uCount > $nCount){ #peaks at this stage represent adjacent bin blocks where UV > non-UV in all bins
                    $i_ or $i_ = $i;
                    $nCountSum += $nCount;
                    $uCountSum += $uCount;
                    if($maxUCount < $uCount){
                        $maxUCount = $uCount;
                        $bestCounts = [$uCount, $nCount];
                    } 
                } else {
                    defined $i_ and checkUveBlock($orientation, $i_, $i - 1);
                    ($nCountSum, $uCountSum, $maxUCount, $i_) = (0, 0, 0, undef);
                }
            }
            defined $i_ and checkUveBlock($orientation, $i_, $maxI);
            mergeBlocks($chrom, $orientation); 
        }
    }   
    close $uveFileH;
    my $fieldNames = "HBID, CHROMOSOME, START_, END_, BINSIZE, NBINS, THRESHOLD, NORMALIZEDDENSITY, STRAND, PVALUE, FRACUV, INGENE";
    loadData($uveFile, $uveTable, ",", $fieldNames);
    scoreHitBlockGenes($uveTable);
}    
sub scaleSamplesUve {
    status("scaling samples so that summed hit count is the same before and after scaling\n");
    foreach my $sample ($nonUVsample, $uvSample){ getUveFragTableInfo($sample) }   
    my $sampleAverage = ($sampleInfo{$nonUVsample}{hitCount} + $sampleInfo{$uvSample}{hitCount}) / 2;              
    $sampleInfo{$nonUVsample}{scalar} = $sampleAverage / $sampleInfo{$nonUVsample}{hitCount};
    $sampleInfo{$uvSample}{scalar} = $sampleAverage / $sampleInfo{$uvSample}{hitCount};
    status("  nonUV scalar = $sampleInfo{$nonUVsample}{scalar}\n  UV scalar = $sampleInfo{$uvSample}{scalar}\n");
}
sub getUveFragTableInfo {
    my ($sample) = @_;
    my $fragsTable = getTableName("Frags", $sample);
    my $mapSourceSql = "SELECT $mapSourceIDSql mapSourceID FROM $fragsTable";
    $mapSourceSql = "SELECT mapSourceID, count(*) hitCount FROM ($mapSourceSql) GROUP BY mapSourceID";
    $mapSourceSql = "SELECT mapSourceID, hitCount FROM ($mapSourceSql) ORDER BY hitCount DESC";  
    runSQL($mapSourceSql, \my($mapSourceID, $hitCount));
    fetchRow();    
    status("  using mapSourceID $mapSourceID for sample $sample, hit count = $hitCount\n");
    $sampleInfo{$sample}{fragsTable} = $fragsTable;
    $sampleInfo{$sample}{mapSourceID} = $mapSourceID;
    $sampleInfo{$sample}{hitCount} = $hitCount; 
} 
sub adjustInputParametersUve {
    status("adjusting input parameters so that binSize and fragmentLength are even, maNBins is odd\n");
    $binSize = $param{binSize};
    $binSize % 2 and $binSize++; 
    $halfBinSize = $binSize / 2; 
    $fragLength = $param{fragmentLength};
    $fragLength % 2 and $fragLength++;
    $halfFragLength = $fragLength / 2; 
    $effFragLength = $fragLength + ($halfFragLength / 2); #hard plus fuzzy call lengths with intensity adjustment 
    $maNBins = $param{maNBins};
    $maNBins % 2 or $maNBins++; 
    $maNAdj = ($maNBins - 1) / 2;
    $maOffset = $maNAdj * $binSize;
}
sub extractTMapCoverage {
    my ($sample, $chrom) = @_;
    status("  generating fragment map for sample $sample\n");
    initializeUveBins($sample);
    my $fragSql = "SELECT fragmentID, position1, strand1, length1
                   FROM $sampleInfo{$sample}{fragsTable}
                   WHERE CHROMOSOME1 = $chrom";
    runSQL($fragSql, \my($fragID, $pos, $strand, $length));
    my $mapSourceID = $sampleInfo{$sample}{mapSourceID};
    status("    processing fragments\n");
    while (fetchRow()){ 
        $fragID =~ m/$mapSourceID$/ or next;
        my ($hardStart, $hardEnd, $fuzzyStart, $fuzzyEnd);
        if($strand == 1){
            $hardStart = $pos;
            $hardEnd = $pos + $fragLength - 1;
            $fuzzyStart = $hardEnd + 1;
            $fuzzyEnd = $hardEnd + $halfFragLength;
        } else {
            $hardStart = $pos + $length - $fragLength;
            $hardEnd = $pos + $length - 1;
            $fuzzyStart = $hardStart - $halfFragLength;
            $fuzzyEnd = $hardStart - 1;
        }
        my $lowBin =  getTMapBin($hardStart);
        my $highBin = getTMapBin($hardEnd);
        for (my $bin = $lowBin; $bin <= $highBin; $bin += $binSize){ 
            my ($countStart, $countEnd) = getTMapCountEnds($bin, $hardStart, $hardEnd);
            countUveBinHits($sample, $countStart, $countEnd, 1, $strand, $bin); #intensity over stretch = 1 in hard call
        }
        $lowBin =  getTMapBin($fuzzyStart);
        $highBin = getTMapBin($fuzzyEnd);
        for (my $bin = $lowBin; $bin <= $highBin; $bin += $binSize){ 
            my ($countStart, $countEnd) = getTMapCountEnds($bin, $fuzzyStart, $fuzzyEnd);
            my $startIntensity = ($countStart - $fuzzyStart) / $halfFragLength;
            my $endIntensity = ($fuzzyEnd - $countEnd) / $halfFragLength;
            $strand == 1 and $startIntensity = 1 - $startIntensity;
            $strand == 2 and $endIntensity = 1 - $endIntensity; 
            my $averageIntensity = ($startIntensity + $endIntensity) / 2;  
            countUveBinHits($sample, $countStart, $countEnd, $averageIntensity, $strand, $bin); #intensity scaled from 1 to 0 over fuzzy call
        }
    }
}
sub initializeUveBins {
    status("    initializing bins\n");
    my ($sample) = @_;
    foreach my $orientation(1..$maxStrand){
        for (my $bin = $minBin; $bin <= $maxBin; $bin += $binSize){ 
            $coverage{$orientation}{$bin}{$sample} = 0;
        }
    }
}
sub getTMapBin { 
    my ($position) = @_;
    return int(($position/$binSize) + 0.5) * $binSize;
}
sub getTMapCountEnds { 
    my ($bin, $fragStart, $fragEnd) = @_;
    my $binStart = $bin - $halfBinSize;
    my $binEnd =   $bin + $halfBinSize - 1;
    my $countStart = $fragStart >= $binStart ? $fragStart : $binStart;
    my $countEnd = $fragEnd <= $binEnd ? $fragEnd : $binEnd;
    return ($countStart, $countEnd);
}
sub countUveBinHits {
    my ($sample, $countStart, $countEnd, $averageIntensity, $strand, $bin) = @_;
    my $effLength = ($countEnd - $countStart + 1) * $averageIntensity;
    my $scaledCount = (($effLength / $effFragLength) * $sampleInfo{$sample}{scalar}) / $maNBins;
    my $orientation = getTranscriptionOrientation($strand);
    for (my $maBin = $bin-$maOffset; $maBin <= $bin+$maOffset; $maBin += $binSize){
          $coverage{$orientation}{$maBin}{$sample} += $scaledCount;
    }
}
sub getTranscriptionOrientation { 
    #strand reversal applies to predicted transcription orientation, NOT to fragment projection from read position!
    my ($strand) = @_;
    $param{stranded} or return 1;
    return $param{reverseStrands} ? ($strand % 2) + 1 : $strand;
}
sub checkUveBlock {
    my ($orientation, $j, $k) = @_;
    #perform two checks for block significance
    my $p = getUveP($uCountSum, $nCountSum); #first check on entire block (highest overall counts)
    $p <= $param{pThreshold} or $p = getUveP(@$bestCounts); #second check on bin with highest UV count
    $p <= $param{pThreshold} and push @blocks, [$j, $k, $uCountSum, $nCountSum, $p];
#    #check all possible combinations of inner adjacent bin blocks for significant enrichment p-value 
#    foreach my $i($j..$k){
#        my ($uCount, $nCount) = (0, 0);
#        foreach my $l($i..$k){
#            $uCount += $coverage{$orientation}{$bins[$l]}{$uvSample};
#            $nCount += $coverage{$orientation}{$bins[$l]}{$nonUVsample};
#            my $p = getUveP($uCount, $nCount); #p-value on the candidate entry set of adjacent bins
#            if($p <= $param{pThreshold}){
#                $p = getUveP($uCountSum, $nCountSum); #final p-value on the entire block
#                push @blocks, [$j, $k, $uCountSum, $nCountSum, $p]; #any significant inner p-value keeps the entire block
#                return;
#            }
#        }  
#    }
}
sub getUveP { #cumulative Poisson p-value calculated externally using fortran executable
    my ($uCount, $nCount) = @_;
    $uCount or $uCount = 0;
    $nCount or $nCount = 0;
    $uCount = int($uCount + 0.5); #round values to be consistent with Poisson method
    $nCount = int($nCount + 0.5);
    $uCount > $nCount or return 1; #don't bother if UV <= non-UV, obviously won't be significant
    my $pCommand = "$uvProb $uCount $nCount";
    my $p = qx/$pCommand/;
    $p =~ m/NaN/ and return 1;
    $p =~ s/\s//g;
    return $p;
}
sub mergeBlocks {
    my($chrom, $orientation) = @_;
    status("  attempting to merge adjacent blocks\n");
    my ($j_, $k_, $uCount_, $nCount_, $p_);
    foreach my $block(@blocks){ #blocks were stored in bin order
        my ($j, $k, $uCount, $nCount, $p) = @$block;
        if($k_){
            my ($j__, $k__) = ($k_ + 1, $j - 1); #$x__ = xGap
            if( (($k__ - $j__ + 1) * $binSize) <= $param{maxGapLength} ){ #merge criterion #1: maxGapLength
                my $uCount__ = getUveBlockCount($uvSample, $orientation, $j__, $k__);
                my $nCount__ = getUveBlockCount($nonUVsample, $orientation, $j__, $k__);
                my $uCountS = $uCount_  + $uCount; #S = separate
                my $nCountS = $nCount_  + $nCount;
                my $uCountM = $uCount_ + $uCount__ + $uCount; #M = merged
                my $nCountM = $nCount_ + $nCount__ + $nCount;
                my $pS = getUveP($uCountS, $nCountS);                  
                my $pM = getUveP($uCountM, $nCountM);
                if($pM < $pS) { #merge criterion #2: merged p value must improve the combined significance
                    ($j_, $k_, $uCount_, $nCount_, $p_) = ($j_, $k, $uCountM, $nCountM, $pM);
                } else {
                    commitUveBlock($chrom, $orientation, $j_, $k_, $uCount_, $nCount_, $p_);
                    ($j_, $k_, $uCount_, $nCount_, $p_) = ($j, $k, $uCount, $nCount, $p);
                }
            } else {
                commitUveBlock($chrom, $orientation, $j_, $k_, $uCount_, $nCount_, $p_);
                ($j_, $k_, $uCount_, $nCount_, $p_) = ($j, $k, $uCount, $nCount, $p);
            }
        } else {
            ($j_, $k_, $uCount_, $nCount_, $p_) = ($j, $k, $uCount, $nCount, $p);
        }  
    } 
    $k_ and commitUveBlock($chrom, $orientation, $j_, $k_, $uCount_, $nCount_, $p_);
}
sub getUveBlockCount {
    my ($sample, $orientation, $j, $k) = @_;
    my $count = 0;
    foreach my $i($j..$k){ $count += $coverage{$orientation}{$bins[$i]}{$sample} }
    return $count;
}
sub commitUveBlock {
    my ($chrom, $orientation, $j, $k, $uCount, $nCount, $p) = @_;
    my $nBins = $k - $j + 1;
    my $uvDensity = $uCount / ($nBins * $binSize); 
    my $uvNormDens = $uvDensity / $uvExpectedDensity;
    my $fracUV = $uCount/($uCount + $nCount);
    my $blockStart = $bins[$j] - $halfBinSize;
    my $blockEnd = $bins[$k] + $halfBinSize;
    $uveID++;    
    print $uveFileH "$uveID,$chrom,$blockStart,$blockEnd,$param{binSize},$nBins,$param{pThreshold},$uvNormDens,$orientation,$p,$fracUV,0\n";
}

1;

