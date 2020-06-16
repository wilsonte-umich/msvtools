#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %gffFields %geneticCode));

my (%wtStats, %mutStats);

sub findCausative{
    my ($childRefSeq, $wtSample, $mutSample) = @_;

    status("getting parameters...\n");
        my $sample = "$wtSample\_$mutSample";         
        my $discsTable = getTableName('Discs', $sample); 
        tableExists($discsTable) or die "could not find table $discsTable";
        getStatistics(getTableName('Stats', $wtSample), \%wtStats);        
        getStatistics(getTableName('Stats', $mutSample), \%mutStats);  
        my $lineSize = getLineSize();

    status("retrieving heterozygous discrepancies...\n"); 
        my $hetRef = getHetMut($discsTable, $lineSize);
        my $candRef = getCandidates($discsTable, $lineSize);
        
    status("finding missing/extra bases consistent with heterozygosity...\n");        
        my $chromLinesRef = loadChromLines($param{refSeq});         
        checkOtherME($discsTable, $wtSample, $mutSample, $hetRef, $chromLinesRef, $lineSize, \&checkMEHeterozygous, \%wtStats, \%mutStats);      
        checkOtherME($discsTable, $wtSample, $mutSample, $candRef, $chromLinesRef, $lineSize, \&checkMECandidate, \%wtStats, \%mutStats); 
        
    status("creating mutated copy of $param{refSeq}...\n");
        mutateChromSeqs($childRefSeq, $hetRef);  
        rectifyRefSeqLines($childRefSeq, $lineSize);               
        my $childChromLinesRef = loadChromLines($childRefSeq); 
        printChromFiles($childChromLinesRef, $childRefSeq);         
        
    status("determining which heterozygous mutations are causative candidates\n");  
        setIsCandidate($hetRef, $candRef);

    status("printing table of mutations by position...\n");
        my $mutsByGeneRef = printMutations($childRefSeq, 'Heterozygous', $hetRef);          
        
    status("printing table of mutations by gene and generating CDS table...\n");
        printGenes($childRefSeq, $lineSize, $chromLinesRef, $childChromLinesRef, 'Heterozygous', $mutsByGeneRef);            
}

sub getHetMut{
    my ($discsTable, $lineSize) = @_;
    return getMutations($lineSize, getHetSQL($discsTable));    
}

sub getCandidates{
    my ($discsTable, $lineSize) = @_;
    my $hetSQL = getHetSQL($discsTable);
    #candidate logic = heterozygous in combined samples, homozygous in mut sample, absent in wt sample
    return getMutations($lineSize, "$hetSQL AND S2FRAC >= $param{minHomoF} AND S1FRAC <= (1 - $param{minHomoF}) "); 
}

sub getHetSQL{
    my ($discsTable) = @_;
    my $discsBaseSQL = getDiscsBaseSQL($discsTable);
    my $hetFilters = getHetFilters();
    return "SELECT CHROMOSOME, POSITION, DISCREPANCYTYPE, EXTRA,
                   S1_CONSISTENTCOUNT, S1_READCOUNT, S2_CONSISTENTCOUNT, S2_READCOUNT, CONSISTENTLOD
            FROM ($discsBaseSQL)
            WHERE $hetFilters";
}

sub getDiscsBaseSQL{
    my ($discsTable) = @_;
    return "SELECT CHROMOSOME, POSITION, DISCREPANCYTYPE, EXTRA, 
                    S1_CONSISTENTCOUNT, S1_READCOUNT, S2_CONSISTENTCOUNT, S2_READCOUNT, CONSISTENTLOD,  
                    NVL((S1_CONSISTENTCOUNT/NULLIF(S1_READCOUNT, 0)), 0) S1FRAC, 
                    (S2_CONSISTENTCOUNT/S2_READCOUNT) S2FRAC,
                    ((S1_CONSISTENTCOUNT + S2_CONSISTENTCOUNT)/(S1_READCOUNT + S2_READCOUNT)) S12FRAC                         
            FROM $discsTable";
}

sub getHetFilters{
    my $sampleFilters = getSampleFilters(\%wtStats, \%mutStats, 1); 
    #heterozygosity logic = NOT homozygous and IS real   
    return " NOT (S12FRAC >= $param{minHomoF} AND S1_CONSISTENTCOUNT > 0 AND S2_CONSISTENTCOUNT > 0)
             AND (S12FRAC >= $param{minRealF} OR (S1FRAC > $param{minHomoF} AND S1_READCOUNT > $wtStats{minCoverageRMap}) OR S2FRAC > $param{minHomoF})
             AND $sampleFilters ";
}

sub checkMEHeterozygous{
    my ($sample1, $sample2, $nReadsRef, $consistentRef) = @_;
    my ($passed, $s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD, $consistent, $nReads) = 
        checkMEHomozygous($sample1, $sample2, $nReadsRef, $consistentRef);
    $passed and return undef;
    $nReads or return undef;   
    my ($s1Frac, $s2Frac, $s12Frac) = (0, 0, 0);
    $s1Reads and $s1Frac = $s1Consistent / $s1Reads;
    $s2Reads and $s2Frac = $s2Consistent / $s2Reads;
    $s12Frac = $consistent / $nReads;
    $passed = ($s12Frac >= $param{minRealF} or
               ($s1Frac > $param{minHomoF} and $$nReadsRef{$sample1} > 2) or
               $s2Frac > $param{minHomoF} );  
    return ($passed, $s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD, $consistent, $nReads, $s1Frac, $s2Frac, $s12Frac);
}    

sub checkMECandidate{
    my ($sample1, $sample2, $nReadsRef, $consistentRef) = @_;
    my ($passed, $s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD, $consistent, $nReads, $s1Frac, $s2Frac, $s12Frac) = 
        checkMEHeterozygous($sample1, $sample2, $nReadsRef, $consistentRef);
    $passed or return undef;
    $passed = ($s2Frac >= $param{minHomoF} and $s1Frac <= (1 - $param{minHomoF}));   
    return ($passed, $s1Consistent, $s1Reads, $s2Consistent, $s2Reads, $LOD, $consistent, $nReads, $s1Frac, $s2Frac, $s12Frac);   
}

sub setIsCandidate{
    my ($hetRef, $candRef) = @_;
    foreach my $chrom(keys %$candRef){     
        foreach my $zrLine(keys %{$$candRef{$chrom}}){    
            foreach my $zrLinePos(keys %{$$candRef{$chrom}{$zrLine}}){     
                $$hetRef{$chrom}{$zrLine}{$zrLinePos}{isCandidate} = 1;
            }
        }
    }
}

1;


