#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs));

my (@samples, %dirs, %genes, %sampleCounts, %readCounts, %density, %scalars);
my $genePadding = 1000;
my $nBins = 100;

sub brduSeq{
    (@samples) = @_;
    getDirsBrdu();
    getGeneInfoBrdu();
    #getSamplesCountsBrdu();
    getReadsBrdu();
    scaleSamplesBrdu();    
    printGenesBrdu();
    printGeneRanksBrdu();  
}
sub getDirsBrdu{
    $dirs{outDir} = "$param{inputPath}/brdu";
    -d $dirs{outDir} or mkdir $dirs{outDir};
    $dirs{csv} = "$dirs{outDir}/csv";
    -d $dirs{csv} or mkdir $dirs{csv};
    $dirs{jpg} = "$dirs{outDir}/jpg";
    -d $dirs{jpg} or mkdir $dirs{jpg};
    -d $dirs{csv} and system("rm $dirs{csv}/*");
    -d $dirs{jpg} and system("rm $dirs{jpg}/*");   
}
sub getGeneInfoBrdu {
    my $sql = "SELECT g.name2, g.chromosome, g.strand, g.start_, g.end_, c.start_, c.end_
               FROM refgene_$param{refSeq} g, cds_$param{refSeq} c      
               WHERE g.name2 = c.name2
               GROUP BY g.name2, g.chromosome, g.strand, g.start_, g.end_, c.start_, c.end_";
    runSQL($sql, \my($name, $chrom, $strand, $start, $end, $exonStart, $exonEnd));
    while (fetchRow()){
        $genes{$name}{chrom} = $chrom;
        $genes{$name}{strand} = $strand; 
        ($genes{$name}{start} and $genes{$name}{start} <= $start) or $genes{$name}{start} = $start;
        ($genes{$name}{end} and $genes{$name}{end} >= $end) or $genes{$name}{end} = $end;
        push @{$genes{$name}{exons}}, [$exonStart, $exonEnd];
    }  
    foreach my $name(keys %genes) {
        $genes{$name}{size} = $genes{$name}{end} - $genes{$name}{start} + 1;
        $genes{$name}{binSize} = int(($genes{$name}{size} / $nBins) + 0.5);
        $genes{$name}{paddedSize} = $genes{$name}{size} + (2 * $genePadding);      
    }   
}
# sub getSamplesCountsBrdu {
    # foreach my $sample (@samples){
        # my $fragsTable = getTableName("FRAGS", $sample);
        # my $uFragsSql = "SELECT chromosome1, position1, strand1
                         # FROM $fragsTable
                         # GROUP BY chromosome1, position1, strand1";
        # my $countSql = "SELECT count(*) N FROM ($uFragsSql)";
        # runSQL($countSql, \my$count);
        # fetchRow();
        # $sampleCounts{$sample} = $count;  
    # }    
# }
# sub scaleSamplesBrdu {
    # my $maxSampleCount = 0;
    # foreach my $sample (@samples){
        # my $sampleCount = $sampleCounts{$sample};
        # $maxSampleCount >= $sampleCount or $maxSampleCount = $sampleCount;  
    # }
    # foreach my $sample (@samples){
        # $scalars{$sample} = $maxSampleCount / $sampleCounts{$sample};  
    # }   
    # ####################
    # foreach my $sample (@samples){
        # print "$sample\t$sampleCounts{$sample}\t$scalars{$sample}\n"; 
    # } 
# }
sub getReadsBrdu {
    foreach my $name(keys %genes) {
        my %gene = %{$genes{$name}};
        foreach my $sample (@samples){
            my $fragsTable = getTableName("FRAGS", $sample);
            my $uFragsSql = "SELECT position1, strand1
                             FROM $fragsTable
                             WHERE chromosome1 = $gene{chrom}
                               AND position1 >= $gene{start} - $genePadding
                               AND position1 <= $gene{end} + $genePadding
                             GROUP BY position1, strand1";
            runSQL($uFragsSql, \my($pos, $strand)); #not using strand info since method does not appear strand specific
            while (fetchRow()){ 
                my $bin = pos2BinBrdu(\%gene, $pos);
                $readCounts{$name}{$bin}{$sample}++;
                $sampleCounts{$name}{$sample}++;
                $sampleCounts{inGene}{$sample}++;
            }    
        }
    }
}
sub scaleSamplesBrdu {
    my $maxSampleCount = 0;
    foreach my $sample (@samples){
        my $sampleCount = $sampleCounts{inGene}{$sample};
        $maxSampleCount >= $sampleCount or $maxSampleCount = $sampleCount;  
    }
    foreach my $sample (@samples){
        $scalars{$sample} = $maxSampleCount / $sampleCounts{inGene}{$sample};  
    }   
    ####################
    foreach my $sample (@samples){
        status("$sample\t$sampleCounts{inGene}{$sample}\t$scalars{$sample}\n"); 
    } 
}
sub printGenesBrdu {
    foreach my $name(keys %genes) {
        my %gene = %{$genes{$name}};
        my $csvFile = "$dirs{csv}/$name.csv";
        open my $csvH, ">", $csvFile or die "could not open $csvFile\n";
        my $samplesLabel = join(",", @samples);
        print $csvH "BIN,CHROMOSOME,STRAND,POSITION,$samplesLabel\n";
        my @bins = sort {$a<=>$b} keys %{$readCounts{$name}};
        my $minBin = $bins[0];
        $minBin or $minBin = 0;
        $minBin > 0 and $minBin = 0;
        my $maxBin = pop @bins;  
        $maxBin or $maxBin = $nBins;
        $maxBin < $nBins and $maxBin = $nBins;
        foreach my $bin ($minBin..$maxBin){ 
            my $pos = bin2PosBrdu(\%gene, $bin);
            print $csvH "$bin,$gene{chrom},$gene{strand},$pos";
            foreach my $sample (@samples){
                my $count = $readCounts{$name}{$bin}{$sample};
                $count or $count = 0;
                $scalars{$sample} and $count *= $scalars{$sample}; 
                print $csvH ",$count";      
            }
            print $csvH "\n";    
        }
        close $csvH;
        drawGeneGraphBrdu($name, \%gene, $csvFile);       
    }
}
sub drawGeneGraphBrdu {
    my ($name, $gene, $csvFile) = @_;   
    my $orientation = "forward";
    $$gene{strand} == 2 and $orientation = "reverse";
    my $title = "$name (Chr$$gene{chrom}, $$gene{start}-$$gene{end}, $orientation)"; 
    system ("cp $csvFile $dirs{jpg}/hits.csv");  
    createExonsFileBrdu($gene); 
    my $scriptFile = "$dirs{jpg}/jpg.R";
    open my $scriptH, ">", $scriptFile or die "could not open $scriptFile\n";
    print $scriptH "samples <- c('".join("','",@samples)."')", "\n";   
    print $scriptH "setwd('$dirs{jpg}')", "\n";
    print $scriptH "hits <- read.table('hits.csv',header=TRUE,sep=',')", "\n";
    print $scriptH "exons <- read.table('exons.csv',header=TRUE,sep=',')", "\n";
    print $scriptH "minX <- min(hits\$BIN)", "\n";
    print $scriptH "maxX <- max(hits\$BIN)", "\n";  
    print $scriptH "maxY <- 0", "\n";  
    print $scriptH "for (sample in samples) {", "\n";
    print $scriptH "    maxY <- max(maxY,max(hits[[sample]]))", "\n";   
    print $scriptH "}", "\n";
    print $scriptH "maxY <- round(maxY,2)", "\n";  
    print $scriptH "maxY95 <- 0.95 * maxY", "\n";  
    print $scriptH "bitmap(file='$name.jpg',type='jpeg',width=1000,height=1000,units='px')", "\n";
    print $scriptH "plot(0,0,xlim=c(minX,maxX),ylim=c(0,maxY),main='$title',
                    xlab='gene bin (bin 0 = tss, bin $nBins = gene end)',
                    ylab='corrected read count in bin',pch='')", "\n";    
    print $scriptH "for (i in 1:length(exons\$START)) {", "\n";
    print $scriptH "    rect(xleft=exons\$START[i],ybottom=maxY95,xright=exons\$END[i],ytop=maxY,col=\"red\")", "\n";  ;  
    print $scriptH "}", "\n";
    print $scriptH "color <- 1", "\n"; 
    print $scriptH "for (sample in samples) {", "\n";
    print $scriptH "    lines(hits\$BIN,hits[[sample]],col=color,lwd=2)", "\n";  
    print $scriptH "    color <- color + 1", "\n";  
    print $scriptH "}", "\n";
    print $scriptH "legend(maxX/2,maxY95,samples,lty=1,lwd=2,col=c(1,2,3))", "\n";

    close $scriptH;
    system("Rscript $scriptFile");
}
sub createExonsFileBrdu {
    my ($gene) = @_;  
    my $exonFile = "$dirs{jpg}/exons.csv";
    open my $exonH, ">", $exonFile or die "could not open $exonFile\n";
    print $exonH "START,END\n";
    my $exons = $$gene{exons};
    foreach my $exon(@$exons) {
        my ($start, $end) = @$exon;
        $start = pos2BinBrdu($gene, $start, 1);
        $end = pos2BinBrdu($gene, $end, 1);
        print $exonH "$start,$end\n";    
    }
    close $exonH;
}
sub printGeneRanksBrdu {
    open my $rankH, ">", "$dirs{outDir}/densityRank.csv";
    my $samplesLabel = join(",", @samples);
    print $rankH "GENE,$samplesLabel\n";   
    foreach my $name(keys %genes) {
        my %gene = %{$genes{$name}};
        print $rankH "$name";
        foreach my $sample (@samples){   
            my $count = $sampleCounts{$name}{$sample};
            $count or $count = 0;
            $scalars{$sample} and $count *= $scalars{$sample}; 
            my $density = $count / $gene{paddedSize} ; #hits per bp
            print $rankH ",$density";
        }
        print $rankH "\n";     
    }  
    close $rankH;
}
sub pos2BinBrdu {
    my ($gene, $pos, $suppressRound) = @_;
    my $bin;
    if ($$gene{strand} == 1) {
        $bin = ($pos - $$gene{start}) / $$gene{binSize};  
    } else {
        $bin = ($$gene{end} - $pos) / $$gene{binSize};
    }  
    $suppressRound or $bin = int($bin + 0.5);
    return $bin;
}
sub bin2PosBrdu {
    my ($gene, $bin) = @_;
    if ($$gene{strand} == 1) {
        return ($bin * $$gene{binSize}) + $$gene{start};
    } else {
        return $$gene{end} - ($bin * $$gene{binSize});
    }   
}

1;
