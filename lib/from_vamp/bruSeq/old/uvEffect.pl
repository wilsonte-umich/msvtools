#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs));

my $binSize = 1000; #bp
my $upstreamBins = 3;
my $upstreamSize = $binSize * $upstreamBins;
my $minHitsPerBin = 3; #3 hits/bin density = 95% chance that a bin will have >= 1 hit
my ($ctlSample, $uvSample, %readTables, $gunoTable);
my (%samples, %scalars, %bins, %genes, %strands, %sizes, %densities, %files);

sub uvEffect {
    ($ctlSample, $uvSample) = @_;  
    prepareReadsUVE();
    prepareGenesUVE();
    # countsReadsUVE();
    # scaleSamplesUVE(); 
    # scoreExpressionUVE();
    # printCountsUVE();
    # runRUVE();
}
sub prepareReadsUVE {
    foreach my $sample ($ctlSample, $uvSample){
        my $fragsTable = "FRAGS_$sample";
        $readTables{$sample} = "READS_$sample\_UVE";
        dropTable($readTables{$sample});
        runSQL("CREATE TABLE $readTables{$sample} AS
                SELECT chromosome1 chromosome, position1 position
                FROM $fragsTable
                GROUP BY chromosome1, position1, strand1");  
    }   
}
# #this version finds genes with NON-divergent promoters
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
                        # and g1.paddedStart <= g2.paddedEnd
                        # and g2.paddedStart <= g1.paddedEnd
                        # and g1.name2 != g2.name2";     
    # my $nonOverlapSQL = "SELECT name2, strand, size_, chromosome, tss, lt, gt
                         # FROM ($uniqSQL) WHERE name2 NOT IN ($overlapSQL)";                            
    # runSQL("CREATE TABLE $gunoTable AS $nonOverlapSQL");
# }
#this version finds genes with divergent promoters
#later code is NOT set to deal with this information
#did not pursue this because
# 1) the divergent promoter pattern quickly disappears at two start sites
#    are separated by even a few kb, consistent with < 1 kb map location
# 2) the pattern strongly depends on having ~equal exprssion of both genes
#thus, in the end, divergent pattern best shown anecdotally
sub prepareGenesUVE {
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
                        AND g1.paddedStart <= g2.paddedEnd
                        AND g2.paddedStart <= g1.paddedEnd
                        AND g1.name2 != g2.name2
                        AND ((g1.strand = 2 AND g2.strand = 1 
                             AND g1.paddedStart < g2.paddedStart
                             AND g1.paddedEnd - g2.paddedStart < 5000)
                             OR(g2.strand = 2 AND g1.strand = 1 
                             AND g2.paddedStart < g1.paddedStart
                             AND g2.paddedEnd - g1.paddedStart < 5000))";  
    my $divergentSQL = "SELECT name2, strand, size_, chromosome, tss, lt, gt
                         FROM ($uniqSQL) WHERE name2 IN ($overlapSQL)";                            
    runSQL("CREATE TABLE $gunoTable AS $divergentSQL");
}
sub countsReadsUVE {
    foreach my $chrom(1..nChrom()) {
    #foreach my $chrom(21..21) {  
        foreach my $sample ($ctlSample, $uvSample){
            my $distSQL = "SELECT g.name2, g.strand, g.size_, g.lt, g.gt, r.position - g.tss distance  
                           FROM $readTables{$sample} r, $gunoTable g
                           where r.chromosome = g.chromosome
                             AND r.chromosome = $chrom";   
            my $inGeneSQL = "SELECT name2, strand, size_, round(distance/$binSize) * $binSize distance,
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
sub scaleSamplesUVE {
    $samples{$ctlSample} = sampleCountUVE($ctlSample);
    $samples{$uvSample} = sampleCountUVE($uvSample);
    my $maxSampleCount = $samples{$ctlSample};
    $maxSampleCount >= $samples{$uvSample} or $maxSampleCount = $samples{$uvSample};  
    $scalars{$ctlSample} = $maxSampleCount / $samples{$ctlSample};
    $scalars{$uvSample} = $maxSampleCount / $samples{$uvSample};  
    scaleSampleUVE($ctlSample);
    scaleSampleUVE($uvSample);   
} 
sub sampleCountUVE {
    my ($sample) = @_;  
    runSQL("SELECT count(*) N from $readTables{$sample}", \my$sampleCount);
    fetchRow();
    return $sampleCount;
}
sub scaleSampleUVE {
    my ($sample) = @_;
    foreach my $name(keys %sizes) { 
        $genes{$sample}{$name} or $genes{$sample}{$name} = 0;
        $genes{$sample}{$name} *= $scalars{$sample};
        foreach my $distance(keys %{$bins{combined}{$name}}) {
            $bins{$sample}{$name}{$distance} or $bins{$sample}{$name}{$distance} = 0;
            $bins{$sample}{$name}{$distance} *= $scalars{$sample};  
        }   
    }
}
sub scoreExpressionUVE {
    foreach my $name(keys %sizes) { 
        my $count = $genes{$ctlSample}{$name};   
        $densities{$name} = $count / $sizes{$name}; #hits per bp for gene in control sample
        my $hitsPerBin = $densities{$name} * $binSize;
        $hitsPerBin < $minHitsPerBin and $sizes{$name} = undef; #block poor expressors to limit divide by zeros
    }  
}
sub printCountsUVE {  
        my $ctlH1 = getFileHUVE($ctlSample, 1);
        my $uvH1 = getFileHUVE($uvSample, 1);
        my $ctlH2 = getFileHUVE($ctlSample, 2);
        my $uvH2 = getFileHUVE($uvSample, 2);  
        foreach my $distance(sort {$a<=>$b} keys %{$bins{header}}) {
            print $ctlH1 $distance;
            print $uvH1 $distance;
            print $ctlH2 $distance;
            print $uvH2 $distance;
            foreach my $name(sort {$a cmp $b} keys %sizes) { 
                $sizes{$name} or next;
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
    return "$param{inputPath}/$sample\_$strand\_UVE.csv";
}
sub getRatioHUVE {
    my ($sample, $name, $distance) = @_;
    my $count = $bins{$sample}{$name}{$distance};
    my $ratio = "";
    if (defined $count) {
        my $binDensity = $count / $binSize;
        $ratio =  $binDensity / $densities{$name};
    } 
    return $ratio;
}
sub runRUVE {
    foreach my $strand(1..2) {
        foreach my $sample ($ctlSample,$uvSample) {
            my $scriptFile = "$param{inputPath}/uve.R";
            open my $scriptH, ">", $scriptFile or die "could not open $scriptFile\n";
            my $file = getFileNameUVE($sample, $strand);
            print $scriptH "data <- read.table('$file',header=TRUE,sep=',')", "\n";  
            print $scriptH "nCols <- length(data)", "\n";  
            print $scriptH "values <- data[,2:nCols]", "\n";    
            print $scriptH "output <- data.frame(data\$DISTANCE)", "\n";       
            print $scriptH "output\$MEAN <- apply(values, 1, mean, na.rm=TRUE)", "\n";   
            print $scriptH "output\$MEDIAN <- apply(values, 1, median, na.rm=TRUE)", "\n";   
            print $scriptH "output\$STDEV <- apply(values, 1, sd, na.rm=TRUE)", "\n"; 
            print $scriptH "write.table(output,file='$file.summary.csv', sep=',')", "\n";  
            # print $scriptH "bitmap(file='$file.jpg',type='jpeg',width=1000,height=1000,units='px')", "\n";
            # print $scriptH "plot(output\$DISTANCE,output\$MEAN,main='$sample',
                            # xlab='distance from TSS (bp)',
                            # ylab='bin ratio')", "\n"; 
            close $scriptH;
            system("Rscript $scriptFile");  
        }
    }
}

# sub drawGeneGraphBrdu {
    # my ($name, $gene, $csvFile) = @_;   
    # my $orientation = "forward";
    # $$gene{strand} == 2 and $orientation = "reverse";
    # my $title = "$name (Chr$$gene{chrom}, $$gene{start}-$$gene{end}, $orientation)"; 
    # system ("cp $csvFile $dirs{jpg}/hits.csv");  
    # createExonsFileBrdu($gene); 
    # my $scriptFile = "$dirs{jpg}/jpg.R";
    # open my $scriptH, ">", $scriptFile or die "could not open $scriptFile\n";
    # print $scriptH "samples <- c('".join("','",@samples)."')", "\n";   
    # print $scriptH "setwd('$dirs{jpg}')", "\n";
    # print $scriptH "hits <- read.table('hits.csv',header=TRUE,sep=',')", "\n";
    # print $scriptH "exons <- read.table('exons.csv',header=TRUE,sep=',')", "\n";
    # print $scriptH "minX <- min(hits\$BIN)", "\n";
    # print $scriptH "maxX <- max(hits\$BIN)", "\n";  
    # print $scriptH "maxY <- 0", "\n";  
    # print $scriptH "for (sample in samples) {", "\n";
    # print $scriptH "    maxY <- max(maxY,max(hits[[sample]]))", "\n";   
    # print $scriptH "}", "\n";
    # print $scriptH "maxY <- round(maxY,2)", "\n";  
    # print $scriptH "maxY95 <- 0.95 * maxY", "\n";  
    # print $scriptH "bitmap(file='$name.jpg',type='jpeg',width=1000,height=1000,units='px')", "\n";
    # print $scriptH "plot(0,0,xlim=c(minX,maxX),ylim=c(0,maxY),main='$title',
                    # xlab='gene bin (bin 0 = tss, bin $nBins = gene end)',
                    # ylab='corrected read count in bin',pch='')", "\n";    
    # print $scriptH "for (i in 1:length(exons\$START)) {", "\n";
    # print $scriptH "    rect(xleft=exons\$START[i],ybottom=maxY95,xright=exons\$END[i],ytop=maxY,col=\"red\")", "\n";  ;  
    # print $scriptH "}", "\n";
    # print $scriptH "color <- 1", "\n"; 
    # print $scriptH "for (sample in samples) {", "\n";
    # print $scriptH "    lines(hits\$BIN,hits[[sample]],col=color,lwd=2)", "\n";  
    # print $scriptH "    color <- color + 1", "\n";  
    # print $scriptH "}", "\n";
    # print $scriptH "legend(maxX/2,maxY95,samples,lty=1,lwd=2,col=c(1,2,3))", "\n";
# 
    # close $scriptH;
    # system("Rscript $scriptFile");
# }

1;
