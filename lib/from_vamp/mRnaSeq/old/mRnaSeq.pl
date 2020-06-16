#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %gffFields %bowtieFields));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{addMrnaMaps} or $command{addMrnaMaps} = ['multiThread', '24:00:00', 5000, 0];

my %strands = ('+' => 1, '-' => 2);
my %reverseStrands = (1 => '+', 2 => '-');
my ($inSample, $outSample, $refSeq, $extractTable);
my (%inDirs, %outDirs, %inMapFiles, %outMapFiles, %outPurgedFiles, $nullHits, %mRnaInfo);

sub addMrnaMaps { 
    #remaps reads that did not align to the genome to the transcriptome
    #merges data into the presumed already existing Pairs or Frags table, as appropriate
    #should only be run AFTER standard genome mapping and extraction (not parsing) steps
    #all coordinates end up as genome coordinates, not mRNA
    #may eventually want to have way of mapping ONLY RNA data, but this could be derived from here
    ($inSample) = @_;
    $outSample = "$inSample\_mRna";
    getDirsMrnaSeq();
    mapReadsMrnaSeq();
    loadExonCorrMrnaSeq();
    enforceMaxHitsMrnaSeq();
    extractMapsMrnaSeq();
}

sub getDirsMrnaSeq {
    getDirectories($inSample, \%inDirs);
    getDirectories($outSample, \%outDirs);
    -d $outDirs{sample} or mkdir $outDirs{sample};
    -d $outDirs{read1} or mkdir $outDirs{read1};
    -d $outDirs{read2} or mkdir $outDirs{read2};
    getMapFiles($inSample, \%inMapFiles);
    getMapFiles($outSample, \%outMapFiles);
    getReadFiles($outSample, 'purgedFasta', \%outPurgedFiles);   
}

sub mapReadsMrnaSeq {
    #MrnaSeq mapping algorithms ignore parameter maxHits and return ALL mappings to refMrna
    #this behavior is essential since gene positions may be represented by many mRna variants
    #application of maxHits occurs after genome positions are assigned and redundancy eliminated
    $refSeq = "refMrna_$param{refSeq}";
    status("mapping unaligned $inSample genome reads to $refSeq...\n");
    my $maxReadN = 2;
    $param{unpaired} and $maxReadN = 1;
    foreach my $readN(1..$maxReadN){ #map the reads, one readN at a time
        my $read = "read$readN";
        status("    working on $read...\n");
        my %notAligned = (pass => "$inDirs{$read}/not_aligned.fa", bowtie => "$inMapFiles{$read}.not_aligned");
        system("cp $notAligned{$param{mapType}} $outPurgedFiles{$read}"); #input = reads that failed alignment to genome
        my %mapReadSubs = (pass => \&mapReadMrnaSeq_pass, bowtie => \&mapReadMrnaSeq_bowtie);      
        &{$mapReadSubs{$param{mapType}}}($read);
    } 
}

sub mapReadMrnaSeq_pass {
    my ($read) = @_;
    my $inputFile = $outPurgedFiles{$read};  #reads as purged fasta
    (-e $inputFile) or die "$inputFile does not exist";    
    my $gffFile = $outMapFiles{$read};  #final catenated data destination
    my $oGFF = "-o $gffFile";
    my $output  = "-gff -info_gff"; #output in gff format, including mismatch info 
    my $refSeqFile = "$param{refSeqPath}/$refSeq.fa";
    my $refdb = "$refSeqFile.$param{longWord}.db";
    my $R = "-R $refdb";  #reference sequence derived from inputs
    my $g = "-g 1";   #allow alignment gaps
    my $pstPath = "$param{passPath}/PST/W$param{shortWord}M1m0G0X0.pst";
    my $iQuery = "-i $inputFile";
    my $matchThreshold = int((($param{readLength}-$param{maxDisc})/$param{readLength})*100);
    my $fid = "-fid $matchThreshold";  #allowed discrepancies communicated as % match
    my $pstThreshhold = (2*$param{shortWord})-$param{maxDisc};  #same as table provided by pass
    my $pst = "-pst $pstPath $pstThreshhold";
    unlink ($gffFile);  #clear any previous mapping output
    runPass(($iQuery, $R, $fid, $g, $pst, $output, $oGFF)); 
}

sub mapReadMrnaSeq_bowtie {
    my ($read) = @_;
    if ($param{maxDisc} > 3){
        status("Bowtie allows up to 3 mismatches, changing maxDisc to 3...\n");
        $param{maxDisc} = 3;
    }
    my $inputFile = $outPurgedFiles{$read};  #reads as purged fasta
    (-e $inputFile) or die "$inputFile does not exist";    
    my $outputFile = $outMapFiles{$read};  #final catenated data destinations
    my $notAligned = "$outputFile.not_aligned";
    unlink ($outputFile, $notAligned);  #clear any previous mapping output
    my $bowtieCommand = "$param{bowtiePath}/bowtie "; 
    $bowtieCommand .= "-f -B 1 ";  #fasta, report refSeq 1-referenced
    $bowtieCommand .= "--un $notAligned ";    
    $bowtieCommand .= " -v $param{maxDisc} ";
    $bowtieCommand .= " -a "; #report all valid alignments
    $bowtieCommand .= " $refSeq ";
    $bowtieCommand .= " $inputFile ";
    $bowtieCommand .= " $outputFile "; 
    $bowtieCommand .= " > $outputFile.stdout";
    runBowtie($bowtieCommand);          
}

sub loadExonCorrMrnaSeq {
    #store the information needed to correct mRna to genome coordinates
    status("loading exon correction values...\n");
    runSQL("SELECT name1, chromosome, strand, corrstart_, end_, mRnaStart
            FROM refgeneexons_$param{refSeq}
            ORDER BY name1, mRnaStart",
            \my($name1, $chrom, $strand, $corrStart, $end, $mrnaStart));
    while (fetchRow()){
        $mRnaInfo{$name1}{$mrnaStart} = [$chrom, $strand, $corrStart, $end];
        push @{$mRnaInfo{$name1}{mRnaStarts}}, $mrnaStart; #sorted low to high via SQL
    }
}

sub enforceMaxHitsMrnaSeq {
    #once mapping is completed, convert map file from mRna to genome coordinate space
    #keep only one hit per genome position per pairID, thereby discarding transcript redundancy
    status("correcting coordinates from mRna to genome space...\n");
    my $maxReadN = 2;
    $param{unpaired} and $maxReadN = 1;
    foreach my $readN(1..$maxReadN){ 
        my $read = "read$readN";
        status("    working on $read...\n");
        enforceMaxHitsMrnaSeqRead($read);
    }
}

sub enforceMaxHitsMrnaSeqRead {
    my ($read) = @_;
	my $corrFile = "$outMapFiles{$read}.corr"; 
	open my $corrFileH, ">", $corrFile or die "could not open $corrFile"; 
    open my $mapFileH, "<", $outMapFiles{$read} or die "could not open $outMapFiles{$read}\n";	
    my $getLineSub = "getLine\_$param{mapType}_mRna";
    my $getLineSubRef = \&{$getLineSub};
    my ($pairID, $prevPairID, %mapKeys);
    while (!$prevPairID){ #since many rRNA reads do not return a pairID
        my $line = <$mapFileH>; 
        $prevPairID = &$getLineSubRef(\$line, \%mapKeys) 
    } 
    while (my $line = <$mapFileH>){     
        $pairID = &$getLineSubRef(\$line, \%mapKeys); 
        if ($pairID) {
            unless ($pairID eq $prevPairID){  
                printCorrMapsMrnaSeq($corrFileH, $prevPairID, \%mapKeys);              
                delete $mapKeys{$prevPairID}; #just to keep memory requirements low  
                $prevPairID = $pairID;   
            } 
        }    
    } 
    printCorrMapsMrnaSeq($corrFileH, $prevPairID, \%mapKeys);  
    close $mapFileH;
    close $corrFileH;
    status("        $nullHits null gene hits, most presumably correspond to rRNA\n");
    system("mv $outMapFiles{$read} $outMapFiles{$read}.mRna"); #keep full mRna mapping just in case
    system("mv $corrFile $outMapFiles{$read}"); #finalize corrected and stratified map file
}

sub printCorrMapsMrnaSeq {
    #enforce maxHits, with full respect to discrepancy strata
    #commit allowable mappings to corrected map file
    my ($corrFileH, $pairID, $mapKeys) = @_;
    my %nDiscs = %{$$mapKeys{$pairID}};
    my @nDiscs = sort {$a <=> $b} keys %nDiscs; #give priority to lowest nDisc
    my @keepLines;
    foreach my $nDisc(@nDiscs){ #check against maxHits, working up through encountered strata
        my %mapKeys = %{$nDiscs{$nDisc}};
        my @mapKeys = keys %mapKeys;
        my $nMapKeys = scalar(@mapKeys); #number of _genome_ hits in strata
        scalar(@keepLines) + $nMapKeys > $param{maxHits} and last;
        foreach my $mapKey(@mapKeys){ push @keepLines, $mapKeys{$mapKey} }
    }
    foreach my $line(@keepLines){ print $corrFileH $$line } #line has newline
}

sub getLine_pass_mRna { #modified to get corrected chrom, position and strand info from exons table
    my ($line, $mapKeys) = @_;
    chomp $$line;
    my @in = split(/\t/, $$line);
    unless ($in[$gffFields{attributes}] =~ m/Name=(.+?);/){return 0}
    my $pairID = $1;
    my $name1 = $in[$gffFields{chromosome}]; #really just the target sequence name
    unless($mRnaInfo{$name1}){ $nullHits++; return undef }
    my $strand = $strands{$in[$gffFields{strand}]}; #the strand relative to the mRNA
    my $pos = $in[$gffFields{start}]; #the position relative to the mRNA    
    my ($chrom, $genomePos, $genomeStrand) = convertExonToGenome($name1, $pos, $strand);
    $genomeStrand or return undef;
    $in[$gffFields{chromosome}] = $reverseRefSeqs{$param{refSeq}}{$chrom};
    $in[$gffFields{start}] = $genomePos;
    $in[$gffFields{strand}] = $reverseStrands{$genomeStrand};
    my ($before, $discs) = getDiscrepancies_pass($in[$gffFields{attributes}]);    
    return getLine_common_mRna($line, $mapKeys, \@in, $pairID, $chrom, $genomePos, $genomeStrand, $discs);       
}

sub getLine_bowtie_mRna { #modified to get corrected chrom, position and strand info from exons table
    my ($line, $mapKeys) = @_;
    chomp $$line;
    my @in = split(/\t/, $$line);
    my $pairID = $in[$bowtieFields{name}];  
    my $name1 = $in[$bowtieFields{chromosome}]; #really just the target sequence name
    unless($mRnaInfo{$name1}){ $nullHits++; return undef }
    my $strand = $strands{$in[$bowtieFields{strand}]}; #the strand relative to the mRNA  
    my $pos = $in[$bowtieFields{position}]; #the position relative to the mRNA     
    my ($chrom, $genomePos, $genomeStrand) = convertExonToGenome($name1, $pos, $strand);
    $genomeStrand or return undef;
    $in[$bowtieFields{chromosome}] = $reverseRefSeqs{$param{refSeq}}{$chrom};
    $in[$bowtieFields{position}] = $genomePos;
    $in[$bowtieFields{strand}] = $reverseStrands{$genomeStrand};
    my ($discs) = getDiscrepancies_bowtie($in[$bowtieFields{mismatches}]);
    return getLine_common_mRna($line, $mapKeys, \@in, $pairID, $chrom, $genomePos, $genomeStrand, $discs);      
}

sub getLine_common_mRna {
    #reassemble the read line based on new genome coordinates and store by strata
    #it is assumed that discs are the same for reads mapping to the same genome coordinates
    my ($line, $mapKeys, $in, $pairID, $chrom, $genomePos, $genomeStrand, $discs) = @_;
    $$line = join("\t", @$in)."\n"; 
    my $nDiscs = getNDiscs($discs);
    my $mapKey = join(":", $chrom, $genomePos, $genomeStrand); #the line that reduces mRna redundancy to genome coordinates
    $$mapKeys{$pairID}{$nDiscs}{$mapKey} or $$mapKeys{$pairID}{$nDiscs}{$mapKey} = $line; #stratify, but only keep one line per map position
    return $pairID;   
}

sub convertExonToGenome {
    #returns  1) the chromosome the gene is on, in VAMP format
    #         2) the position corrected to genome coordinates
    #         3) the read strand corrected for reverse genes
    #uses corrstart_ to account for systematic UCSC refGene table error
    #otherwise get offset depending on how many exons preceded the mapped exon!!!
    my ($name1, $mRnaPos, $strand) = @_;
    my @mRnaStarts = @{$mRnaInfo{$name1}{mRnaStarts}};
    my $i = 0; #find the exon by comparing mRnaPos to mRnaStart of exons
    while ($mRnaStarts[$i + 1] and $mRnaPos >= $mRnaStarts[$i + 1]) { $i++ }
    my $mRnaStart = $mRnaStarts[$i]; 
    my $exonOffset = $mRnaPos - $mRnaStart;
    my ($chrom, $geneStrand, $corrStart, $end) =  @{$mRnaInfo{$name1}{$mRnaStart}};   
    my $genomePos;
    if($geneStrand == 1){ #place the mRnaPos into exon genome coordinates
        $genomePos = $corrStart + $exonOffset;
    } else {
        $genomePos = $end - $exonOffset;
        $strand = ($strand % 2) + 1; #correct the strand for reverse genes
    }
    return ($chrom, $genomePos, $strand);       
}

sub getNDiscs { #this should probably find a home somewhere else...
    my ($discs) = @_;
    $discs or return 0;
    return int( (length($discs) / 4) + 0.5 ); 
}

sub extractMapsMrnaSeq {
    #finally run standard extraction, but putting data into current table and skipping post-extraction processing
    $param{addToExisting} = 1;
    $param{extractDataOnly} = 1;
    my $sample = $outSample.'::'.$inSample; #reverse in/out order to put mRna reads (out) into existing table (in)
    if ($param{unpaired}){ 
        $extractTable = extractUnpaired($sample);
    } else { 
        $extractTable = extractPairs($sample);
    }
}

1;
