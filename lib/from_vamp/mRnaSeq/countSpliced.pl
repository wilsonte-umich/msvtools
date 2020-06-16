#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{xx} or $command{xx} = ['multiThread', '24:00:00', 5000, 0];
defined $param{minJunctionOverlap} or $param{minJunctionOverlap} = 3; 
defined $command{countSpliced} or $command{countSpliced} = ['singleThread', '24:00:00', 13000, 0];

my $delimiter = '@@';

sub countSpliced {
    my ($sample) = @_;
    my ($inSample, $outSample) = splitSampleName($sample);
    my $mapSample = "$outSample\_splice";
    $param{mapType} = 'pass';
    getDirsCountSpliced($inSample, $mapSample, \my%readFiles, \my%outDirs, \my%mapFiles);
    createIntronRefSeq($mapSample); #must do this fresh every time in order to collect %refSeq and %reverseRefSeq
    $param{refSeq} = $mapSample;
    $param{refSeqBase} = $param{refSeq}; 
    mapReadsCountSpliced(\%readFiles, \%mapFiles);
    $param{expandEnds} = 1;
    extractPairsSpliced(\%mapFiles, $outSample); 
}
sub getDirsCountSpliced {
    my ($inSample, $mapSample, $readFiles, $outDirs, $mapFiles) = @_;
    getReadFiles($inSample, 'purgedFasta', $readFiles);  
    getDirectories($mapSample, $outDirs);
    -d $$outDirs{sample} or mkdir $$outDirs{sample};
    -d $$outDirs{read1} or mkdir $$outDirs{read1};
    -d $$outDirs{read2} or mkdir $$outDirs{read2};
    getMapFiles($mapSample, $mapFiles);
}
sub createIntronRefSeq {
    my ($mapSample) = @_; 
    status("creating refSeq $mapSample...\n");
    status("  loading $param{refSeq} chromosome sequences...\n"); #subs found in findHomozygous.pl
    my $lineSize = getLineSize(); 
    my $chromLines = loadChromLines(); 
    getGeneSpliceInfo(\my%chroms, \my%intronStarts, \my%intronEnds);
    commitSliceJunctions($mapSample, \%chroms, \%intronStarts, \%intronEnds, $lineSize, $chromLines);
    $chromLines = undef;
    createDatabase($mapSample);  
}
sub getGeneSpliceInfo {
    my ($chroms, $intronStarts, $intronEnds) = @_;
    status("  retrieving exon and intron information for all genes...\n"); 
    my $intronTable = "refGeneIntrons_Unq_$param{refSeq}";
        runSQL("SELECT name2, chromosome, start_, corrend_ 
                FROM $intronTable",
                \my($gene, $chrom, $intronStart, $intronEnd));          
        while(fetchRow()){
            $$chroms{$gene} = $chrom;
            push @{$$intronStarts{$gene}}, $intronStart;
            push @{$$intronEnds{$gene}}, $intronEnd;
        }    
    #}
}
sub commitSliceJunctions {
    my ($mapSample, $chroms, $intronStarts, $intronEnds, $lineSize, $chromLines) = @_;
    my $refSeqFile = "$param{refSeqPath}/$mapSample/$mapSample.fa";
    open my $refSeqFileH, ">", $refSeqFile or die "could not open $refSeqFile\n";  
    my $i = 0;
    foreach my $gene(keys %$chroms){
        my $chrom = $$chroms{$gene};
        commitIntronStarts($refSeqFileH, $mapSample, \$i, $gene, $chrom, $intronStarts, $lineSize, $chromLines, \my%exonEnds);
        commitIntronEnds($refSeqFileH, $mapSample, \$i, $gene, $chrom, $intronEnds, $lineSize, $chromLines, \my%exonStarts);
        status("  committing exon splice junctions...\n"); 
        foreach my $intronStart(keys %exonEnds){
            foreach my $intronEnd(keys %exonStarts){
                if ($intronEnd > $intronStart){
                   commitSpliceJunction($refSeqFileH, $mapSample, \$i, $gene, $chrom, $intronStart, $intronEnd, $exonEnds{$intronStart}, $exonStarts{$intronEnd});
                }
            } 
        }
    }
    close $refSeqFileH;
}
sub commitIntronStarts {
    my ($refSeqFileH, $mapSample, $i, $gene, $chrom, $intronStarts, $lineSize, $chromLines, $exonEnds) = @_;
    foreach my $intronStart(@{$$intronStarts{$gene}}){
        my $seg1 = getDNASegment($chrom, $intronStart - $param{readLength} + $param{minJunctionOverlap}, $intronStart - 1, $lineSize, $chromLines); 
        my $seg2 = getDNASegment($chrom, $intronStart, $intronStart + $param{readLength} - 1 - $param{minJunctionOverlap}, $lineSize, $chromLines); 
        commitSpliceJunction($refSeqFileH, $mapSample, $i, $gene, $chrom, $intronStart, 0, $seg1, $seg2);    
        $$exonEnds{$intronStart} = $seg1;
    }
}
sub commitIntronEnds {
    my ($refSeqFileH, $mapSample, $i, $gene, $chrom, $intronEnds, $lineSize, $chromLines, $exonStarts) = @_;
    foreach my $intronEnd(@{$$intronEnds{$gene}}){
        my $seg3 = getDNASegment($chrom, $intronEnd - $param{readLength} + 1 + $param{minJunctionOverlap}, $intronEnd, $lineSize, $chromLines);     
        my $seg4 = getDNASegment($chrom, $intronEnd + 1, $intronEnd + $param{readLength} - $param{minJunctionOverlap}, $lineSize, $chromLines); 
        commitSpliceJunction($refSeqFileH, $mapSample, $i, $gene, $chrom, 0, $intronEnd, $seg3, $seg4);
        $$exonStarts{$intronEnd} = $seg4;
    }
}
sub commitSpliceJunction {
    my ($refSeqFileH, $mapSample, $i, $gene, $chrom, $intronStart, $intronEnd, $leftSeg, $rightSeg) = @_;
    $$i++;    
    my $junction = "$gene$delimiter$chrom$delimiter$intronStart$delimiter$intronEnd";
    print $refSeqFileH ">$junction\n$leftSeg$rightSeg\n";
    $refSeqs{$mapSample}{$junction} = $$i;
    $reverseRefSeqs{$mapSample}{$$i} = $junction;
}
sub mapReadsCountSpliced {
    my ($readFiles, $mapFiles) =  @_;
    status("mapping reads to $param{refSeq}...\n");
    my $maxReadN = 2;
    $param{unpaired} and $maxReadN = 1;
    foreach my $readN(1..$maxReadN){ #map the reads, one readN at a time
        my $read = "read$readN";
        status("    working on $read...\n");
        my $inputFile = $$readFiles{$read};  #reads as purged fasta
        (-e $inputFile) or die "$inputFile does not exist";    
        my $gffFile = $$mapFiles{$read};  #final catenated data destination  
        my $oGFF = "-o $gffFile";
        my $output  = "-gff -info_gff"; #output in gff format, including mismatch info 
        my $refSeqFile = "$param{refSeqPath}/$param{refSeq}/$param{refSeq}.fa";
        my $refdb = "$refSeqFile.$param{longWord}.db";
        my $R = "-R $refdb";  #reference sequence derived from inputs
        my $g = "-g 1";   #allow alignment gaps
        my $pstPath = "$param{passPath}PST/W$param{shortWord}M1m0G0X0.pst";
        my $iQuery = "-i $inputFile";
        my $matchThreshold = int((($param{readLength}-$param{maxDisc})/$param{readLength})*100);
        my $fid = "-fid $matchThreshold";  #allowed discrepancies communicated as % match
        my $pstThreshhold = (2*$param{shortWord})-$param{maxDisc};  #same as table provided by pass
        my $pst = "-pst $pstPath $pstThreshhold";
        unlink ($gffFile);  #clear any previous mapping output
        runPass(($iQuery, $R, $fid, $g, $pst, $output, $oGFF)); 
    } 
}
sub extractPairsSpliced {
    my ($mapFiles, $outSample) = @_;
    status("extracting read mappings...\n");    
    my $outTable = "SPLICE_$outSample";
    my $outTableTMP = "$outTable\_tmp";    
    dropTable($outTable);
    dropTable($outTableTMP);
    runSQL("CREATE TABLE $outTableTMP
            (NAME2 VARCHAR2(255), CHROMOSOME NUMBER, INTRONSTART NUMBER, INTRONEND NUMBER)");
    my $outFile = "$outTable.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile\n";
    my $maxReadN = 2;
    $param{unpaired} and $maxReadN = 1;
    foreach my $readN(1..$maxReadN){
        my $read = "read$readN";
        status("    working on $read...\n"); 
        open my $mapFileH, "<", $$mapFiles{$read} or die "could not open $$mapFiles{$read}";
        my $line = <$mapFileH>;
        my ($prevPairID, $map) = getLine_pass(\$line);  
        stratifyIntronMap($map, \my%strata);    
        while (my $line = <$mapFileH>){     
            my ($pairID, $map) = getLine_pass(\$line);
            if ($pairID) {
                unless ($pairID eq $prevPairID){
                    commitSplicedReads($outFileH, \%strata);
                    %strata = ();
                    $prevPairID = $pairID; 
                }
                stratifyIntronMap($map, \%strata);    
            }     
        }   
    }  
    close $outFileH;
    loadData($outFile, $outTableTMP, ',', "NAME2, CHROMOSOME, INTRONSTART, INTRONEND");    
    runSQL("CREATE TABLE $outTable AS
            SELECT NAME2, CHROMOSOME, INTRONSTART, INTRONEND, count(*) CROSSINGREADS
            FROM $outTableTMP
            GROUP BY NAME2, CHROMOSOME, INTRONSTART, INTRONEND");
    dropTable($outTableTMP);   
    my $posTable = "$outTable\_pos";
    dropTable($posTable);   
    runSQL("CREATE TABLE $posTable AS
            SELECT name2, chromosome, intronStart position, 1 isIntronStart,
                    sum(case when intronEnd  = 0 then crossingReads else 0 end) unspliced,
                    sum(case when intronEnd != 0 then crossingReads else 0 end) spliced,
                    0 fractionSpliced
             FROM $outTable
             WHERE intronStart > 0
             GROUP BY name2, chromosome, intronStart
             UNION ALL
             SELECT name2, chromosome, intronEnd position, 0 isIntronStart,
                    sum(case when intronStart  = 0 then crossingReads else 0 end) unspliced,
                    sum(case when intronStart != 0 then crossingReads else 0 end) spliced,
                    0 fractionSpliced
             FROM $outTable
             WHERE intronEnd > 0
             GROUP BY name2, chromosome, intronEnd");  
    runSQL("UPDATE $posTable SET fractionSpliced = round(spliced/nullif((spliced + unspliced),0),5)");      
}
sub stratifyIntronMap {
    my ($map, $strata) = @_;
    my ($junction, $pos, $length, $strand, $discs) = split(":", $$map);  
    my $nDiscs = int(length($discs/4) + 0.5);    
    push @{$$strata{$nDiscs}}, $junction;
}
sub commitSplicedReads {
    my ($outFileH, $strata) = @_;
    my @nDiscs = sort {$a <=> $b} keys %$strata;
    my $nDiscs = $nDiscs[0]; 
    scalar(@{$$strata{$nDiscs}}) > 1 and return; #enforce single best read, ignore reads with multiple hits in best stratum
    my $junction = ${$$strata{$nDiscs}}[0];
    $junction = $reverseRefSeqs{$param{refSeq}}{$junction};
    my ($gene, $chrom, $intronStart, $intronEnd) = split($delimiter, $junction);
    print $outFileH "$gene,$chrom,$intronStart,$intronEnd\n";                          
}

1;

