#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %fieldNames));

#callable parameters and commands in this script
defined $param{tssThresholds} or $param{tssThresholds} = '0.65'; #the minimum fraction of UV hits for calling
defined $param{tssMinBinFactor} or $param{tssMinBinFactor} = 1; #multiplier for the average UV hit density required for consideration
defined $command{findTSS} or $command{findTSS} = ['singleThread', '24:00:00', 2000, 0];

sub findTSS { 
    my ($nonUVsample, $uvSample) = @_;
    #establishes the fraction of all hits within a bin that correspond
    #to UV hits vs. non-UV hits
    #fractions over some threshold identify transcription start sites
    status("finding transcription start sites...\n");
    my @tssThresholds = split(",", $param{tssThresholds});
    $param{binSizes} or $param{binSizes} = $param{binSize};
    my @binSizes = split(",", $param{binSizes});    
    my ($nonUVHitsTable, $uvHitsTable, $nonUVScalar, $uvScalar, $uvHitCount) = getTSSScalars($nonUVsample, $uvSample);
    status("  nonUV scalar = $nonUVScalar\n  UV scalar = $uvScalar\n");
    my $genomeSize = getGenomeSize();
    my $tssTable = newTable("HB", "$nonUVsample\_tss"); #output is a hit block table
    my $tssFile = "$tssTable.csv";
    open my $tssFileH, ">", $tssFile or die "could not open $tssFile: $!\n";
    my $tssID = 0;
    my ($nAgg, $uAgg) = ("decode(nCount,null,0,1)", "decode(uCount,null,0,1)");
    $param{keepDups} and ($nAgg, $uAgg) = ("nCount", "uCount");
    foreach my $tssThreshold(@tssThresholds){ 
        status("  frac UV threshold = $tssThreshold\n");
        foreach my $binSize(@binSizes){ 
            status("    bin size = $binSize\n"); 
            my $minBinHits = int((($uvHitCount*$uvScalar)/$genomeSize)*$binSize*$param{tssMinBinFactor});
            status("    minBinHits = $minBinHits\n      chrom"); 
            foreach my $chrom(1..nChrom()){
                status(" $chrom");
                my $nonUVHitsSQL = "SELECT position, count_ FROM $nonUVHitsTable WHERE chromosome = $chrom";
                my $uvHitsSQL = "SELECT position, count_ FROM $uvHitsTable WHERE chromosome = $chrom";
                my $hitJoinSQL = "SELECT round(nvl(n.position,u.position) / $binSize) * $binSize bin, 
                                         n.count_ * $nonUVScalar nCount, 
                                         u.count_ * $uvScalar uCount
                                  FROM ($nonUVHitsSQL) n FULL OUTER JOIN ($uvHitsSQL) u
                                    ON (n.position = u.position)";
                my $joinCountSQL = "SELECT bin, nvl(sum($nAgg),0) nCount, nvl(sum($uAgg),0) uCount
                                    FROM ($hitJoinSQL)
                                    GROUP BY bin";
                my $hitFracSQL = "SELECT bin, uCount/(uCount + nCount) fracUV
                                  FROM ($joinCountSQL)  
                                  WHERE uCount > $minBinHits
                                  ORDER BY bin"; #only need to enforce uCount since looking for high UV fractions                                  
                runSQL($hitFracSQL, \my($bin,$fracUV));
                my ($sum, $count, $blockStart, $prevBin);    
                while(fetchRow()){ 
                    if ($fracUV >= $tssThreshold){
                        if($blockStart and $bin != $prevBin + $binSize){
                            commitHitBlock($tssFileH, $tssThreshold, $binSize, 
                                           $chrom, $blockStart, $prevBin, 
                                           $sum, $count, \$tssID);
                            $blockStart = 0;
                            $sum = 0;
                            $count = 0;                 
                        }  
                        $blockStart or $blockStart = $bin;
                        $sum += $fracUV;
                        $count++;
                    } else {
                        $blockStart and commitHitBlock($tssFileH, $tssThreshold, $binSize, 
                                                       $chrom, $blockStart, $prevBin, 
                                                       $sum, $count, \$tssID);
                        $blockStart = 0;
                        $sum = 0;
                        $count = 0;
                    }  
                    $prevBin = $bin;
                }
                $blockStart and commitHitBlock($tssFileH, $tssThreshold, $binSize, 
                                               $chrom, $blockStart, $prevBin, 
                                               $sum, $count, \$tssID);
            }    
            status("\n");  
             
        }       
    }    
    close $tssFileH;
    loadData($tssFile, $tssTable, ",", $fieldNames{HB});
    scoreHitBlockGenes($tssTable);
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
