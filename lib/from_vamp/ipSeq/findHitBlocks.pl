#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %fieldNames %refSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $param{ndThresholds} or $param{ndThresholds} = 0; 
defined $command{findHitBlocks} or $command{findHitBlocks} = ['multiThread', '12:00:00', 2000, 0];

sub findHitBlocks { 
    #known possible limitation:
    #  HMAP does not carry strand information at present (although HITS table does)
    #  therefore hit blocks are NOT strand-specific
    my ($sample) = @_;
    status("finding transcribed genome regions...\n");
    $param{ndThresholds} or die "must specify parameter ndThresholds\n";
    my @ndThresholds = split(",", $param{ndThresholds});
    $param{binSizes} or $param{binSizes} = $param{binSize};
    my @binSizes = split(",", $param{binSizes});   
    my $hbTable = newTable("HB", $sample);
    my $hbFile = "$hbTable.csv";
    open my $hbFileH, ">", $hbFile;
    my $hbID = 0;
    foreach my $ndThreshold(@ndThresholds){    
        foreach my $binSize(@binSizes){ 
            my $hMapTable = getTableName('HMap', "$sample\_$binSize");  
            foreach my $chrom(1..nChrom()){ 
                runSQL("SELECT POSITION, NORMALIZEDDENSITY
                        FROM $hMapTable
                        WHERE chromosome = $chrom
                        ORDER BY POSITION",
                        \my($pos,$normDens));
                my ($ndSum, $ndCount, $blockStart);    
                while(fetchRow()){
                    if ($normDens >= $ndThreshold){
                        $blockStart or $blockStart = $pos;
                        $ndSum += $normDens;
                        $ndCount++;
                    } else {
                        $blockStart and commitHitBlock($hbFileH, $ndThreshold, $binSize, 
                                                    $chrom, $blockStart, $pos - $binSize,
                                                    $ndSum, $ndCount, \$hbID);
                        $blockStart = 0;
                        $ndSum = 0;
                        $ndCount = 0;
                    }  
                }
                $blockStart and commitHitBlock($hbFileH, $ndThreshold, $binSize, 
                                            $chrom, $blockStart, $pos,
                                            $ndSum, $ndCount, \$hbID);
            }    
        }       
    }    
    close $hbFileH;
    loadData($hbFile, $hbTable, ",", $fieldNames{HB});
    scoreHitBlockGenes($hbTable);
}

sub commitHitBlock {
    my ($fileH, $threshold, $binSize, $chrom, $blockStart, $blockEnd, $sum, $count, $id) = @_;
    my $mean = $sum / $count;
    $$id++;
    my $halfBinSize = int($binSize/2);
    $blockStart -= $halfBinSize;
    $blockEnd += $halfBinSize;
    print $fileH "$$id,$threshold,$binSize,$chrom,$blockStart,$blockEnd,$count,$mean,0\n";
}

sub scoreHitBlockGenes {
    my ($hbTable) = @_;
    status("determining which transcribed regions cross known genes...\n");
    my $refGeneTable = getRefGeneTableName('All');
    runSQL("UPDATE $hbTable SET ingene = 0");
    runSQL("UPDATE $hbTable 
            SET ingene = 1 
            WHERE hbid IN 
            (SELECT b.hbid
            FROM $hbTable b, $refGeneTable g
            WHERE b.chromosome = g.chromosome
              AND b.start_ <= g.end_
              AND g.corrstart_ <= b.end_)");
}

1;


