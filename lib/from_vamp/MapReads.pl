#!/usr/bin/perl -w
use strict;
use warnings;

#########################################################################
#mapReads.pl is a wrapper around PASS (http://pass.cribi.unipd.it/cgi-bin/pass.pl)
#It executes read mapping by a specific automated sequence
#in which an iteration of PASS is first run at high discrepancy tolerance.
#Reads with too many hits are then retested at progressively
#lower discrepancy tolerance.  Benchmarking shows this gives the highest
#density of useful read mapping with a reasonable time efficiency.
#########################################################################

use vars(qw(%param %types %fields %refSeqs $paramCat));

sub createDatabase{
    my ($refSeq) = @_;
    $refSeq or $refSeq = $param{refSeq};
    #stores pass database version of refSeq for later use
    #only ever called once for each refSeq and longWord combination
    status("creating PASS database...");
        my $refSeqFile = "$param{refSeqPath}/$refSeq/$refSeq.fa";
        my $refdb = "$refSeqFile.$param{longWord}.db";
        my $d = "-d $refSeqFile";    
        my $D = "-D $refdb";
        my $p = "-p $param{longWord}";  #long word for initial db search
        my $flc = "-flc 3";   #sets pass to use low linguistic complexity filter
        runPass(($d, $D, $p, $flc));
        return $refdb;
}

sub createBowtieIndex{
    my ($refSeq) = @_;
    $refSeq or $refSeq = $param{refSeq};
    my $bowtieCommand = "$param{bowtiePath}bowtie-build -f $param{refSeqPath}/$refSeq/$refSeq.fa $param{refSeqPath}/$refSeq/$refSeq"; 
    runBowtie($bowtieCommand);    
}

sub createBWAIndex{
    my ($refSeq) = @_;
    $refSeq or $refSeq = $param{refSeq};
    my $bwaCommand = "$param{bwaPath}bwa index -a bwtsw $param{refSeqPath}/$refSeq/$refSeq.fa";
    runBWA($bwaCommand);    
}

sub prepareReads{
    my ($sample, $getUserSequence) = @_; #$getUserSequence only used by internal vamp calls to override default
    my ($sampleIn, $sampleOut) = splitSampleName($sample);
    getDirectories($sampleOut, \my%dirsOut);
    mkdirBasic(\%dirsOut);
    $getUserSequence or $getUserSequence = \&getUserSequence;
    my $readType = createFastaFile($sampleIn, $sampleOut, $getUserSequence);
    if ($readType eq 'fasta') {  #true if readType was fasta or fastq/qseq was converted to fasta
        status("purging duplicate reads...\n");
            purgeDuplicates($sampleOut);
            $readType = 'purgedFasta';
    } else {
        die "'$readType' is not a recognized readType"
    }    
}

sub createFastaFile {
    my ($sampleIn, $sampleOut, $getUserSequence) = @_;
    my $readType = $param{readType};
    if ($readType eq 'solexa_fastq' or  #convert to fastq to fasta if needed
        $readType eq 'fastq') {         #'or' other fastq input types here as needed...
        status("converting fastq to fasta...\n");
            fastq_to_fastaSample($sampleIn, $sampleOut, $readType, $getUserSequence);
            $readType = 'fasta';
    } elsif ($readType eq 'qseq'){
        status("converting qseq.txt to fasta...\n");
            qseq_to_fastaSample($sampleIn, $sampleOut, $getUserSequence);
            $readType = 'fasta';
    } elsif ($readType eq 'fasta') {
        status("applying read restrictions to fasta reads...\n");
            fasta_to_fastaSample($sampleIn, $sampleOut, $getUserSequence);
    }
    return $readType;
}


sub fasta_to_fastaSample{  
    my ($sampleIn, $sampleOut, $getUserSequence) = @_;
    getReadFiles($sampleIn, 'fasta', \my%fastaFilesIn);
    getReadFiles($sampleOut, 'fasta', \my%fastaFilesOut);
    fasta_to_fastaRead($fastaFilesIn{read1}, $fastaFilesOut{read1}, $getUserSequence);
    $param{unpaired} or fasta_to_fastaRead($fastaFilesIn{read2}, $fastaFilesOut{read2}, $getUserSequence);
}
sub fasta_to_fastaRead{
    my ($fastaFileIn, $fastaFileOut, $getUserSequence) = @_;
    status("  $fastaFileIn...\n");
        if ($fastaFileIn eq $fastaFileOut){
            my $fastaTmp = $fastaFileIn;
            $fastaFileIn .= ".original";
            system("mv $fastaTmp $fastaFileIn");
        } 
        open my $inH, "<", $fastaFileIn or die "could not open $fastaFileIn";
        open my $outH, ">", $fastaFileOut or die "could not open $fastaFileOut";
        while (my $nameLine = <$inH>){
            $nameLine =~ m/^>/ or next;
            my $seq = <$inH>;
            $seq = &$getUserSequence($seq);
            $seq or next;
            print $outH $nameLine;
            print $outH $seq."\n"; 
        }
        close $inH;
        close $outH;
}

sub fastq_to_fastaSample{  
    my ($sampleIn, $sampleOut, $readType, $getUserSequence) = @_;
    getReadFiles($sampleIn, $readType, \my%fastqFilesIn);
    getReadFiles($sampleOut, 'fasta', \my%fastaFilesOut);
    fastq_to_fastaRead($fastqFilesIn{read1}, $fastaFilesOut{read1}, $getUserSequence);
    $param{unpaired} or fastq_to_fastaRead($fastqFilesIn{read2}, $fastaFilesOut{read2}, $getUserSequence);
}
sub fastq_to_fastaRead{
    my ($fastqFileIn, $fastaFileOut, $getUserSequence) = @_;
    status("  $fastqFileIn\n");
        open my $inH, "<", $fastqFileIn or die "could not open $fastqFileIn";
        open my $outH, ">", $fastaFileOut or die "could not open $fastaFileOut";
        while (my $nameLine = <$inH>){
            $nameLine =~ m/^\@/ or next;
            my $seq = <$inH>;            
            my $discard1 = <$inH>;
            my $discard2 = <$inH>;
            $seq = &$getUserSequence($seq);
            $seq or next;
            print $outH ">$nameLine";  
            print $outH $seq."\n"; 
        } 
        close $inH;
        close $outH;
}

our %qseqColumns = (machine=>0,run=>1,lane=>2,
                    tile=>3,x=>4,y=>5,
                    index=>6,read=>7,
                    sequence=>8,quality=>9,filter=>10);
sub qseq_to_fastaSample{  
    my ($sampleIn, $sampleOut, $getUserSequence) = @_;
    getDirectories($sampleIn, \my%dirsIn);
    getReadFiles($sampleOut, 'fasta', \my%fastaFilesOut);    
    open my $outH1, ">", $fastaFilesOut{read1} or die "could not open $fastaFilesOut{read1}";
    open my $outH2, ">", $fastaFilesOut{read2} or die "could not open $fastaFilesOut{read2}";
    my $nPassed;    
    foreach my $tileFile1 (<$dirsIn{qseq1}/*qseq.txt>){ #read files from qseq subdir
        if ($tileFile1 =~ m/.*\/(s_\d)_\d_(\d\d\d\d)_qseq.txt$/){
            my $lane = $1;  #learn lane and tile name from file name
            my $tile = $2;
            my $tileFile2 = "$dirsIn{qseq2}/$lane\_2_$tile\_qseq.txt";
            my %reads;
            qseqGetReads($tileFile1, \%reads, 1, $getUserSequence);  #fill quality filtered reads into hash
            $param{unpaired} or qseqGetReads($tileFile2, \%reads, 2, $getUserSequence);
            foreach my $readID (keys %reads){  #create correlated fasta files for pairs with two passed reads
                if ( defined $reads{$readID}{1} and ($param{unpaired} or defined $reads{$readID}{2}) ){
                    $nPassed++;
                    print $outH1 ">$readID\n";
                    print $outH1 "$reads{$readID}{1}\n";
                    unless ($param{unpaired}){
                        print $outH2 ">$readID\n";
                        print $outH2 "$reads{$readID}{2}\n";                    
                    }
                }
            }            
        }
    }
    close $outH1;
    close $outH2;
    status("    $nPassed pairs passed the quality filter on both reads and were kept\n");

}
sub qseqGetReads{
    my ($qseqFile, $readsRef, $read, $getUserSequence) = @_;
    status("  reading $qseqFile\n");
        open my $inH, "<", $qseqFile or die "could not open $qseqFile";
        while (<$inH>){
            chomp $_;
            my @in = split(/\t/, $_);
            if ($in[$qseqColumns{filter}]){ #at present just use Bustards filter for quality passing
                my $readID = join(":", @in[$qseqColumns{machine}..$qseqColumns{y}]);
                $$readsRef{$readID}{$read} = &$getUserSequence($in[$qseqColumns{sequence}]);  
                
            }
        }
        close $inH;
}

sub getUserSequence { #only take the bases asked for by the user
    my ($seq) = @_;
    chomp $seq;            
    return substr($seq, $param{baseOffset}, $param{readLength});  
}

sub getSampleReadLength { #overrides readLength to the entirety of sample input reads
    my ($inSample) = @_;
    status("discovering $inSample readLength\n");
    if ($param{readType} eq 'solexa_fastq' or  
        $param{readType} eq 'fastq' or
        $param{readType} eq 'fasta'){
        getReadFiles($inSample, $param{readType}, \my%readFilesIn);   
        open my $inH, "<", $readFilesIn{read1} or die "could not open $readFilesIn{read1}: $!\n";
        my $nameLine = <$inH>;
        my $seq = <$inH>;
        chomp $seq;
        $param{readLength} = length($seq);
        close $inH;
    } elsif ($param{readType} eq 'qseq'){
        getDirectories($inSample, \my%dirsIn);
        foreach my $tileFile1 (<$dirsIn{qseq1}/*qseq.txt>){ #read files from qseq subdir
            if ($tileFile1 =~ m/.*\/(s_\d)_\d_(\d\d\d\d)_qseq.txt$/){
                open my $inH, "<", $tileFile1 or die "could not open $tileFile1 $!\n";
                my $line = <$inH>;
                chomp $line;
                my @in = split(/\t/, $line);
                my $seq = $in[$qseqColumns{sequence}];
                $param{readLength} = length($seq);
                close $inH; 
                last;
            }
        }   
    } else {
        die "'$param{readType}' is not a recognized readType\n";
    }    
    status("  $inSample readLength = $param{readLength}\n");
}

my %mapBinFileHs; #larger datasets require using disk to prevent memory overruns during dup purging
sub purgeDuplicates{  #reject identical read pairs, within limitation that must be 100% identical
    my ($sample) = @_;
    getReadFiles($sample, 'fasta', \my%fastaFiles);
    getReadFiles($sample, 'purgedFasta', \my%purgedFiles);
    
    $param{keepDups} and return skipDupPurge(\%fastaFiles, \%purgedFiles);
    
    my ($inCounter, $outCounter, $countsRef);

    if($param{unpaired}){
        ($inCounter, $outCounter, $countsRef) = purgeDuplicatesUnpaired($sample, \%fastaFiles, \%purgedFiles);
    } else {
        ($inCounter, $outCounter, $countsRef) = purgeDuplicatesPairedReads($sample, \%fastaFiles, \%purgedFiles);
    }
    status("\n  $inCounter\tread pairs processed, yielding...\n");
    status("  $outCounter\tunique read pairs\n");
    
    getDirectories($sample, \my%dirs);
    my $countFile = "$dirs{sample}/pairCounts.xls";
    
    status("  writing pair count frequencies to\n    $countFile\n");
        open my $countFileH, ">", $countFile or die "could not open $countFile";
        foreach my $count (sort {$a <=> $b} keys %$countsRef){print $countFileH "$count\t$$countsRef{$count}\n"}
        close $countFileH;           
}

sub skipDupPurge {
    my ($fastaFilesRef, $purgedFilesRef) = @_;
    status("  skipping duplicate purge, copying and renaming reads...\n");
    skipDupPurgeRead($fastaFilesRef, $purgedFilesRef, 'read1');
    $param{unpaired} or skipDupPurgeRead($fastaFilesRef, $purgedFilesRef, 'read2');
    return 1;
}

sub skipDupPurgeRead {
    my ($fastaFilesRef, $purgedFilesRef, $read) = @_;
    status("  writing new purged fasta file\n    $$purgedFilesRef{$read}\n");
    open my $fastaH, "<", $$fastaFilesRef{$read} or die "could not open $$fastaFilesRef{$read}";
    open my $outH, ">", $$purgedFilesRef{$read} or die "could not open $$purgedFilesRef{$read}\n";
    my $outCounter = 0;
    while (my $line = <$fastaH>){
        $line =~ m/^>/ and next;
        length($line) > 1 or next; #only non-blank, non-name, i.e. sequence lines acted on
        $outCounter++;
        print $outH ">$outCounter\n";  #lose old pair names, just use an incrementing numerical pairID
        print $outH $line;  #newline is propogated from original file
    }
    close $outH; 
    close $fastaH;
    status("    $outCounter reads copied\n");
}

sub purgeDuplicatesUnpaired{  
    my ($sample, $fastaFilesRef, $purgedFilesRef) = @_;
    
    #TODO: purgeDuplicatesUnpaired not yet updated to use disk  
    #binning for memory efficiency during duplicate purging
    
    status("  collecting unique reads...\n");
        open my $fasta1H, "<", $$fastaFilesRef{read1} or die "could not open $$fastaFilesRef{read1}";
        my ($read, $inCounter, %reads); 
        while (my $line1 = <$fasta1H>){
            if ($line1 =~ m/^>/){
                $read = 1;
                $inCounter++;
            } elsif ($read) {
                $reads{$line1}++; 
                $read = undef;
            }
        }
        close $fasta1H;
    
    status("  writing new purged fasta file\n    $$purgedFilesRef{read1}\n ");
        open my $out1H, ">", $$purgedFilesRef{read1} or die "could not open $$purgedFilesRef{read1}\n";
        my ($outCounter, %counts) = 0;
        foreach my $read (keys %reads){
            $outCounter++;
            print $out1H ">$outCounter\n";  #lose old pair names, just use an incrementing numerical pairID
            print $out1H $read;  #newline is propogated from original file
            $counts{$reads{$read}}++;
        }
        close $out1H;

    return ($inCounter, $outCounter, \%counts);
}

sub purgeDuplicatesPairedReads{  #reject identical read pairs, within limitation that must be 100% identical  
    my ($sample, $fastaFilesRef, $purgedFilesRef) = @_;
    status("  collecting unique read pairs...\n");
        open my $fasta1H, "<", $$fastaFilesRef{read1} or die "could not open $$fastaFilesRef{read1}";
        open my $fasta2H, "<", $$fastaFilesRef{read2} or die "could not open $$fastaFilesRef{read2}";
        my $read; 
        my $inCounter = 0;
        while (<$fasta1H>){
            my $line1 = $_;
            my $line2 = <$fasta2H>;
            chomp $line1;
            chomp $line2;
            if ($line1 =~ m/^>/){
                my ($read1) = split("/", $line1);
                my ($read2) = split("/", $line2);
                ($read1 eq $read2) or die "corresponding read order in fasta files violated!";
                $read = $read1;
                $inCounter++;
            } elsif ($read) {
                my ($seq1, $seq2) = sort {$a cmp $b} ($line1, $line2);  #order sequences to account for reversed duplicates
                my $catSeq = "$seq1:$seq2";
                my $binSeq = substr($catSeq, 0, $param{dupBinSize});
                my $binFileH = getMapBinFileH($sample, $binSeq);
                print $binFileH "$catSeq\n"; 
                $read = 0;
            }
        }
        close $fasta1H;
        close $fasta2H;
        closeMapBinFileHs();

    status("  writing new purged fasta files\n    $$purgedFilesRef{read1}\n    $$purgedFilesRef{read2}");   
        open my $out1H, ">", $$purgedFilesRef{read1} or die "could not open $$purgedFilesRef{read1}\n";
        open my $out2H, ">", $$purgedFilesRef{read2} or die "could not open $$purgedFilesRef{read2}\n";
        my $outCounter = 0;
        my %counts;    
        foreach my $binSeq (keys %mapBinFileHs) { #purge duplicate pairs within each binned temporary read file
            my %pairs;
            my $binFile = getMapBinFile($sample, $binSeq);
            open my $binFileH, "<", $binFile or die "could not open $binFile\n";
            while(my $catSeq = <$binFileH>){
                chomp $catSeq;
                $pairs{$catSeq}++; #hold unique read pairs in hash and count the number of occurrences
            }
            close $binFileH;
            unlink $binFile;
            foreach my $catSeq(keys %pairs){
                $outCounter++;
                my ($seq1, $seq2) = split (":", $catSeq); 
                print $out1H ">$outCounter\n";  #lose old pair names, just use an incrementing numerical pairID
                print $out1H "$seq1\n";  
                print $out2H ">$outCounter\n";
                print $out2H "$seq2\n";
                $counts{$pairs{$catSeq}}++;            
            }
        }    
        close $out1H;
        close $out2H;    
        
    return ($inCounter, $outCounter, \%counts);      
}


sub getMapBinFileH {
    my ($sample, $binSeq) = @_;
    unless($mapBinFileHs{$binSeq}){
        my $binFile = getMapBinFile($sample, $binSeq);
        open my $binFileH, ">", $binFile or die "could not open $binFile\n";
        $mapBinFileHs{$binSeq} = $binFileH;
    }
    return $mapBinFileHs{$binSeq};
}
sub getMapBinFile {
    my ($sample, $binSeq) = @_;
    return "$sample\_$binSeq.txt";
}
sub closeMapBinFileHs {
    foreach my $binSeq (keys %mapBinFileHs) {
        my $binFileH =  $mapBinFileHs{$binSeq};
        close $binFileH;
    }
}

sub mapReads{
    my ($inputSample, $outputSample, $read) = @_;
    my $paramFile = "$outputSample\_map\_parameters.txt";     
    printParameters($paramFile, 'mapReads', $inputSample, $outputSample, $read);     
    if ($param{mapType} eq 'pass'){
        mapReads_PASS($inputSample, $outputSample, $read);
    } elsif ($param{mapType} eq 'bowtie') {
        mapReads_Bowtie($inputSample, $outputSample, $read);
    } elsif ($param{mapType} eq 'bwa') {
        mapReads_BWA($inputSample, $outputSample, $read);
    } else {
        die "$param{mapType} is not a recognized mapType";
    }
}

sub mapReads_PASS{
    my ($inputSample, $outputSample, $read) = @_;
    open my $statusFileH, ">", "progress.txt";
    reportStatus($statusFileH, "mapReads begun...");
    
    #pass parameters which are held constant
    getReadFiles($inputSample, 'purgedFasta', \my%purgedFiles);
    getMapFiles($outputSample, \my%gffFiles);
    my $inputFile = $purgedFiles{$read};  #reads as purged fasta
    (-e $inputFile) or die "$inputFile does not exist";    
    my $gffFile = $gffFiles{$read};  #final catenated data destination
    my $gffFileTMP = "$gffFile.TMP"; #holds data from each interation of pass
    my $oGFF = "-o $gffFileTMP";
    my $output  = "-gff -info_gff -aligned -not_aligned"; #output in gff format, including mismatch info
    my $refSeqFile = "$param{refSeqPath}/$param{refSeq}/$param{refSeq}.fa";
    my $refdb = "$refSeqFile.$param{longWord}.db";
    my $R = "-R $refdb";  #reference sequence derived from inputs
    my $g = "-g 1";   #allow alignment gaps
    my $maxHits = $param{maxHits} + 1; #+1 allows keeping of all up to and including maxHits
    my $maxBestHits = "-max_best_hits $maxHits";   
    my $pstPath = "$param{passPath}PST/W$param{shortWord}M1m0G0X0.pst";
    my $aligned = "aligned.fa";
    my $notAligned = "not_aligned.fa";
    my $tooManyHits = "too_many_hits.fa";
    my %kept;  #counter for successful read mapping
    unlink ($gffFile, $notAligned, $tooManyHits);  #clear any previous mapping output
    
    #initialize first run of pass against entire input file
    my $iQuery = "-i $inputFile";
    
    #loop through possible discrepancy thresholds in DESCENDING order
    for (my $NDisc = $param{maxDisc}; $NDisc>=0; $NDisc--){
        
        reportStatus($statusFileH, "  nDisc = $NDisc");
        
        #calculate pass parameters dependent on discrepancy threshold and run pass
        print "-" x 60; print "\n$NDisc mismatches\n"; print "-" x 60;
        my $matchThreshold = int((($param{readLength}-$NDisc)/$param{readLength})*100);
        my $fid = "-fid $matchThreshold";  #allowed discrepancies communicated as % match
        my $pstThreshhold = (2*$param{shortWord})-$NDisc;  #same as table provided by pass
        my $pst = "-pst $pstPath $pstThreshhold";
        unlink ($gffFileTMP, $aligned); #clear temp results of previous pass iteration
        runPass(($iQuery, $R, $fid, $g, $maxBestHits, $pst, $output, $oGFF));
        
        reportStatus($statusFileH, "    PASS mapping done");
               
        #score the aligned output against the discrepancy threshold
        #(it would be easier if pass did this for you, but it goes pretty fast)
        open my $gffFileTMPH, "<", $gffFileTMP;
        open my $gffFileH, ">>", $gffFile; #APPEND to catenated gff output file
        my %retry = ();
        while (<$gffFileTMPH>){
            if ($_ =~ m/Hits=$maxHits;/){ #read _exceeded_ maxHits, retry these at lower NDisc
                if ($_ =~ m/Name=(.+?);/){$retry{$1}=1}
            } else {  #read within discrepancy tolerance limits
                print $gffFileH $_;  #keep all mappings for the read
                if ($_ =~ m/Name=(.+?);/){$kept{$1}++}  #count all kept reads
            } 
        }
        close $gffFileH;
        close $gffFileTMPH;
        
        #retrieve just those aligned reads which exceeded the discrepancy threshold
        #(it would be easier if pass did this for you, but it goes pretty fast)
        open my $alignedH, "<", $aligned;
        open my $tooManyHitsH, ">", $tooManyHits; #write NEW retry file
        while (<$alignedH>){
            if ($_ =~ m/^>(.+)\n/){ #a name line
                if (defined $retry{$1}){ #name is on retry list, dump to retry file
                    print $tooManyHitsH $_;
                    my $seqLine = <$alignedH>;
                    print $tooManyHitsH $seqLine;
                }
            } 
        }
        close $tooManyHitsH;
        close $alignedH;
        
        reportStatus($statusFileH, "    PASS output parsed");

        #reset input files to the retry file
        $iQuery = "-i $tooManyHits";
    }
    unlink ($gffFileTMP, $aligned); #final output = catenated alignments as gff, not_aligned.fa and too_many_hits.fa
    my $kept = scalar(keys %kept);
    status("\n\n$kept reads kept over all iterations\n");
    reportStatus($statusFileH, "mapReads done\n$kept reads kept over all iterations\n");
    close $statusFileH;
    return $kept;
}

sub mapReads_Bowtie{
    my ($inputSample, $outputSample, $read) = @_;
    
    open my $statusFileH, ">", "progress.txt";
    reportStatus($statusFileH, "mapReads_Bowtie begun...");    
    
    if ($param{maxDisc} > 3){
        status("Bowtie allows up to 3 mismatches, changing maxDisc to 3...\n");
        $param{maxDisc} = 3;
    }

    #parameters which are held constant
    getReadFiles($inputSample, 'purgedFasta', \my%purgedFiles);
    my $inputFile = $purgedFiles{$read};  #reads as purged fasta    
    (-e $inputFile) or die "$inputFile does not exist";     
    getMapFiles($outputSample, \my%mapFiles);   
    my $outputFile = $mapFiles{$read};  #final catenated data destinations
    my $notAligned = "$outputFile.not_aligned";
    my $tooManyHits = "$outputFile.too_many_hits";
    my $indeterminate = "$outputFile.indeterminate";
    unlink ($outputFile, $notAligned, $tooManyHits, $indeterminate);  #clear any previous mapping output

    #loop through possible discrepancy thresholds in DESCENDING order
    my (@outputFiles, @tooManyHitsFiles, $tooManyHits_nDisc, @indeterminateFiles);
    for (my $nDisc = $param{maxDisc}; $nDisc>=0; $nDisc--){
        
        reportStatus($statusFileH, "  nDisc = $nDisc");
        
        my $indeterminate_nDisc = "$indeterminate.$nDisc";
        $tooManyHits_nDisc = "$tooManyHits.$nDisc";        
        my $output_nDisc = "$outputFile.$nDisc";
        $nDisc == $param{maxDisc} or push @indeterminateFiles, $indeterminate_nDisc;
        push @tooManyHitsFiles, $tooManyHits_nDisc;        
        push @outputFiles, $output_nDisc;

        #calculate pass parameters dependent on discrepancy threshold and run pass
        print "-" x 60; print "\n$nDisc mismatches\n"; print "-" x 60;
        my $bowtieCommand = "$param{bowtiePath}bowtie "; 
        $bowtieCommand .= "-f -B 1 ";  #fasta, report refSeq 1-referenced
        my $unFile = $indeterminate_nDisc;
        $nDisc == $param{maxDisc} and $unFile = $notAligned; 
        $bowtieCommand .= "--un $unFile ";  
        $bowtieCommand .= "--max $tooManyHits_nDisc ";
        $bowtieCommand .= " -v $nDisc ";
        $bowtieCommand .= " -k $param{maxHits} ";        
        $bowtieCommand .= " -m $param{maxHits} ";
        $bowtieCommand .= " $param{refSeqPath}/$param{refSeq}/$param{refSeq} ";
        $bowtieCommand .= " $inputFile ";
        $bowtieCommand .= " $output_nDisc "; 
        runBowtie($bowtieCommand);        
        
        reportStatus($statusFileH, "    Bowtie done");
        
        $inputFile = $tooManyHits_nDisc;
    }
    
    system("cat @indeterminateFiles > $indeterminate");
    unlink (@indeterminateFiles);  
    system("cat $tooManyHits_nDisc > $tooManyHits");
    unlink (@tooManyHitsFiles);    
    system("cat @outputFiles > $outputFile");
    unlink (@outputFiles);      
    
    reportStatus($statusFileH, "mapReads done\n");
    close $statusFileH;
}

sub mapReads_BWA{
    my ($inputSample, $outputSample, $read) = @_;
    $param{readType} eq 'solexa_fastq' or die "mapType 'bwa' requires readType 'solexa_fastq'\n";
    getReadFiles($inputSample, 'solexa_fastq', \my%fastqFiles);
    getMapFiles($outputSample, \my%mapFiles);   
    my $fastqFile = $fastqFiles{$read};
    my $saiFile = $mapFiles{$read};
    my $bwaCommand = "$param{bwaPath}bwa aln -n $param{maxDisc} -R $param{maxHits} -I ";
    $bwaCommand .= "$param{refSeqPath}/$param{refSeq}/$param{refSeq}.fa $fastqFile > $saiFile";
    runBWA($bwaCommand);    
}

sub pairBWAMaps{
    my ($sample) = @_;
    my ($inputSample, $outputSample) = splitSampleName($sample);
    getReadFiles($inputSample, 'solexa_fastq', \my%fastqFiles);
    getMapFiles($outputSample, \my%mapFiles);   
    my $fastqFile1 = $fastqFiles{read1};
    my $fastqFile2 = $fastqFiles{read2};
    my $saiFile1 = $mapFiles{read1};
    my $saiFile2 = $mapFiles{read2};
    my $bamFile = $mapFiles{sample};
    my $bwaCommand = "$param{bwaPath}bwa sampe -P -r '\@RG\tID:$outputSample' ";
    $bwaCommand .= "$param{refSeqPath}/$param{refSeq}/$param{refSeq}.fa $saiFile1 $saiFile2 $fastqFile1 $fastqFile2";
    $bwaCommand .= " | samtools view -Shb -o $bamFile -";
    runBWA($bwaCommand);    
}

sub reportStatus{
    my ($statusFileH, $message) = @_;
    my ($sec, $min, $hr, $day, $month, $year) = localtime(time);
    $year = $year + 1900;
    $month++;
    print $statusFileH "$message, $month/$day/$year at $hr:$min:$sec\n";
}

sub runPass{
    my (@passParameters) = @_;
    status("\nrunning PASS with parameters:\n  ");
    status(join("\n  ", @passParameters)."\n\n");
    system("$param{passPath}bin/pass @passParameters");
}

sub runBowtie{
    my ($commandLine) = @_;
    status("\nrunning Bowtie as:\n$commandLine\n\n  ");
    system($commandLine);
}

sub runBWA{
    my ($commandLine) = @_;
    status("\nrunning bwa as:\n$commandLine\n\n  ");
    system($commandLine);
}

1;

