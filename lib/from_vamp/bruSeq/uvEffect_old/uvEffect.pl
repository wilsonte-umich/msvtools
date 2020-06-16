#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#binSize is required but already defined by VAMP
defined $param{upstreamBins} or $param{upstreamBins} = 5; #number of bins preceding the gene start
defined $param{minBinHits} or $param{minBinHits} = 3; #3 hits/bin = 95% chance that a bin will have >= 1 hit
defined $param{maxGeneDistance} or $param{maxGeneDistance} = 100000; #3 hits/bin = 95% chance that a bin will have >= 1 hit
defined $command{trackUV} or $command{trackUV} = ['singleThread', '10:00:00', 2000, 0];

#derivative and internal variables
my $upstreamSize = $param{binSize} * $param{upstreamBins};
my ($ctlSample, $uvSample, %hitsTables, $gunoTable);
my (%samples, %scalars, %bins, %genes, %strands, %sizes, %densities, %files);

sub trackUV {
    ($ctlSample, $uvSample) = @_;  
    $hitsTables{$ctlSample} = getTableName('Hits', $ctlSample);
    $hitsTables{$uvSample} = getTableName('Hits', $uvSample);
    prepareGenesUVE();
    countsReadsUVE();
    scaleSamplesUVE(); 
    scoreExpressionUVE();
    printCountsUVE();
    my @fileRoots = calcStatsUVE();
    plotCurvesUVE(@fileRoots);
    reportGeneNumbers();
}
sub prepareGenesUVE { #find genes with that are uniquely defined and do not overlap other genes
    status("finding non-overlapping gene set...\n");
    my $regGeneTable = "refGene_$param{refSeq}";
    $gunoTable = "$regGeneTable\_guno"; #guno = grouped, unique, non-overlapping
    dropTable($gunoTable);
    my $groupSQL = "SELECT name2, chromosome, min(start_) start_, max(end_) end_, strand
                    FROM $regGeneTable
                    GROUP BY name2, chromosome, strand";
    my $countSQL = "SELECT name2, count(*) N FROM ($groupSQL) GROUP BY name2";
    my $singleSQL = "SELECT name2 FROM ($countSQL) WHERE N = 1";      
    my $uniqSQL  = "SELECT name2, chromosome, strand, end_ - start_ size_,
                           decode(strand,1,start_,end_) tss, 
                           decode(strand,1,end_ - start_,$upstreamSize) lt,
                           decode(strand,1,$upstreamSize,end_ - start_) gt,
                           decode(strand,1,start_-$upstreamSize,start_) paddedStart, 
                           decode(strand,1,end_,end_+$upstreamSize) paddedEnd
                    FROM ($groupSQL)
                    WHERE name2 IN ($singleSQL)";      
    my $overlapSQL = "SELECT g1.name2
                      FROM ($uniqSQL) g1, ($uniqSQL) g2
                      WHERE g1.chromosome = g2.chromosome
                        and g1.paddedStart <= g2.paddedEnd
                        and g2.paddedStart <= g1.paddedEnd
                        and g1.name2 != g2.name2";     
    my $nonOverlapSQL = "SELECT name2, strand, size_, chromosome, tss, lt, gt
                         FROM ($uniqSQL) WHERE name2 NOT IN ($overlapSQL)";                            
    runSQL("CREATE TABLE $gunoTable AS $nonOverlapSQL");
}
sub countsReadsUVE { #loop chromosomes and count hits per bin, gene, etc.
    status("counting hits per gene bin...\n");
    foreach my $chrom(1..nChrom()) {
    #foreach my $chrom(21..21) {  
        foreach my $sample ($ctlSample, $uvSample){
            my $distSQL = "SELECT g.name2, g.strand, g.size_, g.lt, g.gt, r.position - g.tss distance  
                           FROM $hitsTables{$sample} r, $gunoTable g
                           where r.chromosome = g.chromosome
                             AND r.chromosome = $chrom";   
            my $inGeneSQL = "SELECT name2, strand, size_, round(distance/ $param{binSize}) * $param{binSize} distance,
                                    case when distance * decode(strand,1,1,-1) > 0 then 1 else 0 end inGene
                             FROM ($distSQL)
                             WHERE distance <= lt 
                               AND distance >= -gt";       
            runSQL($inGeneSQL, \my($name, $strand, $size, $distance, $inGene));
            while (fetchRow()){ 
                $bins{$sample}{$name}{$distance}++;
                $bins{combined}{$name}{$distance}++;
                $bins{header}{$distance}++;
                $strands{$name} = $strand;
                $sizes{$name} = $size;   
                $inGene and $genes{$sample}{$name}++;
            }    
        }
    }  
}
sub scaleSamplesUVE { #normalize samples onto the same scale for inter-sample comparisons
    status("scaling samples to common scale...\n");
    $samples{$ctlSample} = sampleCountUVE($ctlSample);
    $samples{$uvSample} = sampleCountUVE($uvSample);
    my $maxSampleCount = $samples{$ctlSample}; #scale to sample with greatest hit count
    $maxSampleCount >= $samples{$uvSample} or $maxSampleCount = $samples{$uvSample};  
    $scalars{$ctlSample} = $maxSampleCount / $samples{$ctlSample};
    $scalars{$uvSample} = $maxSampleCount / $samples{$uvSample};  
    scaleSampleUVE($ctlSample);
    scaleSampleUVE($uvSample);   
} 
sub sampleCountUVE { #scaling is done based on total chromosome hits per sample
    my ($sample) = @_;  
    runSQL("SELECT count(*) N from $hitsTables{$sample}", \my$sampleCount);
    fetchRow();
    return $sampleCount;
}
sub scaleSampleUVE {
    my ($sample) = @_;
    foreach my $name(keys %sizes) { 
        $genes{$sample}{$name} or $genes{$sample}{$name} = 0; #null bins assigned zero counts
        $genes{$sample}{$name} *= $scalars{$sample};
        foreach my $distance(keys %{$bins{combined}{$name}}) {
            $bins{$sample}{$name}{$distance} or $bins{$sample}{$name}{$distance} = 0;
            $bins{$sample}{$name}{$distance} *= $scalars{$sample};  
        }   
    }
}
sub scoreExpressionUVE { #expression scored as the hit density of the control non-UV sample
    status("scoring relative expression of control sample...\n");
    foreach my $name(keys %sizes) { 
        my $count = $genes{$ctlSample}{$name};   
        $densities{$name} = $count / $sizes{$name}; #density = hits per bp
        my $hitsPerBin = $densities{$name} * $param{binSize};
        $hitsPerBin < $param{minBinHits} and $sizes{$name} = undef; #block poor expressors as uninformative
    }  
}
sub printCountsUVE { #output file giving the relative bin density for all distances (rows) by genes (cols)
    status("printing output...\n");
        my $ctlH1 = getFileHUVE($ctlSample, 1);  #generate output files for each sample and gene orientation
        my $uvH1 = getFileHUVE($uvSample, 1);
        my $ctlH2 = getFileHUVE($ctlSample, 2);
        my $uvH2 = getFileHUVE($uvSample, 2);  
        foreach my $distance(sort {$a<=>$b} keys %{$bins{header}}) {
            print $ctlH1 $distance;
            print $uvH1 $distance;
            print $ctlH2 $distance;
            print $uvH2 $distance;
            foreach my $name(sort {$a cmp $b} keys %sizes) { 
                $sizes{$name} or next; #block poor expressors as uninformative
                my ($ctlH, $uvH) = ($ctlH1, $uvH1); 
                $strands{$name} == 2 and ($ctlH, $uvH) = ($ctlH2, $uvH2);
                my $ctlRatio = getRatioHUVE($ctlSample, $name, $distance);
                my $uvRatio = getRatioHUVE($uvSample, $name, $distance);
                print $ctlH ",$ctlRatio";
                print $uvH ",$uvRatio";
            }   
            print $ctlH1 "\n";
            print $uvH1 "\n";
            print $ctlH2 "\n";
            print $uvH2 "\n";        
        }  
        close $ctlH1;
        close $uvH1;
        close $ctlH2;
        close $uvH2;       
}
sub getFileHUVE {
    my ($sample, $strand) = @_;
    my $file = getFileNameUVE($sample, $strand);
    open my $fileH, ">", $file or die "could not open $file\n";  
    print $fileH "DISTANCE";
    foreach my $name(sort {$a cmp $b} keys %sizes) { 
        $sizes{$name} or next;
        $strands{$name} == $strand or next;
        print $fileH ",$name" 
    }
    print $fileH "\n";    
    return $fileH;   
}
sub getFileNameUVE {
    my ($sample, $strand) = @_;
    my $fileRoot = getFileRootUVE($sample, $strand);
    return "$fileRoot.csv";
}
sub getFileRootUVE {
    my ($sample, $strand) = @_;
    return "$param{inputPath}/$sample\_$strand\_$param{binSize}\_UVE";
}
sub getRatioHUVE { #bind ratio, or relative density = bin hit density / gene hit density for the control sample
    my ($sample, $name, $distance) = @_;
    my $count = $bins{$sample}{$name}{$distance};
    my $ratio = "";
    if (defined $count) {
        my $binDensity = $count / $param{binSize}; #density = hits per bp
        $ratio =  $binDensity / $densities{$name}; #denom determined previously for control sample
    } 
    return $ratio;
}
sub calcStatsUVE { #use R to calculate the mean, median and stdev of all genes in each sample/strand combination
    status("calculating mean and median density per bin...\n");
    my $rScript = "$param{vampPath}/bin/bru/uvEffectStats.R";
    my @fileRoots;
    foreach my $strand(1..2) {
        foreach my $sample ($ctlSample,$uvSample) {
            my $fileRoot = getFileRootUVE($sample, $strand);
            push @fileRoots, $fileRoot;
            #R script takes since file root as args
            system("Rscript $rScript $fileRoot"); #script generates sample_UVE_stats.csv output
        }
    }  
    return @fileRoots
}
sub plotCurvesUVE { #use R to plot all sample/strand combinations to a single plot
    my (@fileRoots) = @_;
    status("creating final plots...\n");
    my $rScript = "$param{vampPath}/bin/bru/uvEffectPlot.R";
    my $fileRoots = join(" ", @fileRoots);
    #R script takes -maxX, maxX, binSize, outputPath and a list of files roots as args
    system("Rscript $rScript -$param{maxGeneDistance} $param{maxGeneDistance} $param{binSize} $param{inputPath} $fileRoots");  
}
sub reportGeneNumbers { 
    my %geneCounts;
    foreach my $name(keys %sizes) { $sizes{$name} and $geneCounts{$strands{$name}}++ }
    foreach my $strand(1..2){ status("\n$geneCounts{$strand} genes on strand $strand\n") }
}

# this version finds genes with divergent promoters
# later code is NOT yet set to deal with this information
# did not pursue this because
# 1) the divergent promoter pattern quickly disappears as two start sites
# are separated by even a few kb, consistent with < 1 kb map location
# 2) the pattern strongly depends on having ~equal exprssion of both genes
# thus, in the end, divergent pattern best shown anecdotally
# sub prepareGenesUVE {
    # my $regGeneTable = "refGene_$param{refSeq}";
    # $gunoTable = "$regGeneTable\_guno"; #guno = grouped, unique, non-overlapping
    # dropTable($gunoTable);
    # my $groupSQL = "SELECT name2, chromosome, min(start_) start_, max(end_) end_, strand
                    # FROM $regGeneTable
                    # GROUP BY name2, chromosome, strand";
    # my $countSQL = "SELECT name2, count(*) N FROM ($groupSQL) GROUP BY name2";
    # my $singleSQL = "SELECT name2 FROM ($countSQL) WHERE N = 1";      
    # my $uniqSQL  = "SELECT name2, chromosome, strand, end_ - start_ size_,
                           # decode(strand,1,start_,end_) tss, 
                           # decode(strand,1,end_ - start_,$upstreamSize) lt,
                           # decode(strand,1,$upstreamSize,end_ - start_) gt,
                           # decode(strand,1,start_-$upstreamSize,start_) paddedStart, 
                           # decode(strand,1,end_,end_+$upstreamSize) paddedEnd
                    # FROM ($groupSQL)
                    # WHERE name2 IN ($singleSQL)";     
    # my $overlapSQL = "SELECT g1.name2
                      # FROM ($uniqSQL) g1, ($uniqSQL) g2
                      # WHERE g1.chromosome = g2.chromosome
                        # AND g1.paddedStart <= g2.paddedEnd
                        # AND g2.paddedStart <= g1.paddedEnd
                        # AND g1.name2 != g2.name2
                        # AND ((g1.strand = 2 AND g2.strand = 1 
                             # AND g1.paddedStart < g2.paddedStart
                             # AND g1.paddedEnd - g2.paddedStart < 5000)
                             # OR(g2.strand = 2 AND g1.strand = 1 
                             # AND g2.paddedStart < g1.paddedStart
                             # AND g2.paddedEnd - g1.paddedStart < 5000))";  
    # my $divergentSQL = "SELECT name2, strand, size_, chromosome, tss, lt, gt
                         # FROM ($uniqSQL) WHERE name2 IN ($overlapSQL)";                            
    # runSQL("CREATE TABLE $gunoTable AS $divergentSQL");
# }



1;
