use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir $rLibDir
            $TRUE $FALSE);
require "$libDir/Illumina.pl";
require "$libDir/formats.pl";
our ($method, $tmpDir, $maxMem,
     $arrayDir, $arrayFmt, $nameCol, $arrayType, $genomeFasta,
     $binSize, 
     $maxTLen, $minQual, $reference, $maxCN,
     $chroms,
     $outputDir, $modelName);
our ($inputDir);
my @globs;
our (%chrCN, @chroms, %chroms);
my ($refChrom, $refStart, $refEnd, $refPloidy,
    $refMean, $refStdev, @refCounts);
my ($cntH, $mdlH, $prbH, $cntFile, $mdlFile, $binFile);
my $chrom;           # the chromosome being processed
my $binE;            # last base of the first possible bin
my $thsBinC = 0;     # coverage of the current bin
my $nxtBinC = 0;     # coverage of next bin (for pairs that cross bin boundaries)
my $prec = 1e-2;     # sufficient significant digits for bins
my ($binWidth, $minBinCount, $maxBinCount, $maxBinI);
my ($nPrbTotal, $probeBedFile) = (0);

# manage options
sub setOptions_train {
    setOptionValue(\$method,     'method');
    setOptionValue(\$tmpDir,     'tmp-dir',         '/tmp');
    setOptionValue(\$maxMem,     'max-mem',         1000000000);
    #-------------------------------------
    setOptionValue(\$arrayDir,   'array-dir');
    setOptionValue(\$arrayFmt,   'array-format', 'Illumina');
    setOptionValue(\$nameCol,    'name-column',  'DNA_ID');
    setOptionValue(\$arrayType,  'array-type');
    setOptionValue(\$genomeFasta,'genome-fasta');
    #-------------------------------------
    setOptionValue(\$binSize,    'bin-size',        5000);
    setOptionValue(\$maxTLen,    'max-TLen',        1000);
    setOptionValue(\$minQual,    'min-qual',        5);
    setOptionValue(\$reference,  'reference');
    setOptionValue(\$maxCN,      'max-copy-number', 4);
    #-------------------------------------
    setOptionValue(\$chroms,     'chromosomes');
    #-------------------------------------
    setOptionValue(\$outputDir,  'output-dir');
    setOptionValue(\$modelName,  'model-name');
    #-------------------------------------
    -d $arrayDir or die "$error: $arrayDir does not exist or is not a directory\n";
    -f $genomeFasta or die "$error: genome fasta not found: $genomeFasta\n";
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem =~ m|\D| and die "$error: max-mem must be an integer number of bytes\n";    
    $binSize =~ m|\D| and die "$error: bin-size must be an integer number of base pairs\n";
    $maxCN =~ m|\D| and die "$error: max-mem must be an integer number\n";    
    if($method eq 'bam'){
        $reference or die "$error: --reference required for method 'bam'\n";
        $reference =~ m|(\w+):(\d+)-(\d+)/(\d+)| or die "malformed reference, expected format 'chrom:start-end/copy_number'\n";
        ($refChrom, $refStart, $refEnd, $refPloidy) = ($1, $2, $3, $4);
        $maxTLen < $binSize or die "$error: --max-TLen must be less than --bin-size\n";
    } elsif($method eq 'bed6'){
        $globs[0] or die "missing bed file as last command argument";
    } else {
        ($method eq 'M' or $method eq 'F') or die "$error: unknown method '$method', expected (bam|M|F)\n";
        $chrCN{chrX} = $method eq 'M' ? 1 : 2;
        $chrCN{chrY} = $method eq 'M' ? 1 : 0;   
    }
    @chroms = sort {$a cmp $b} split(",", $chroms);
    %chroms = map { $_ => 1 } @chroms;    
    -d $outputDir or die "$error: directory not found: $outputDir\n";
}
# main execution block
sub msvtools_train {
    (@globs) = @_;
    
    # initialize
    setOptions_train();
    print STDERR "$utility train: " . getTime(), "\n";
    
    # # open outputs
    $mdlFile = getOutFile('model', $modelName, 'bed.bgz');
    $binFile = getOutFile('bins',  $modelName, 'bed.bgz');
    $probeBedFile = getOutFile('probes', $arrayType, 'bed.gz');

    # train the bin model based on bam counts
    if($method eq 'bam'){
        print STDERR "training model from bam reads\n";
        $cntFile = getOutFile('bin_counts', $modelName, 'gz');
        openOutputStream($cntFile, \$cntH, $TRUE);
        getBinCoverage_train();
        closeHandles($cntH);
        $refMean  = mean (@refCounts);
        $refStdev = stdev(@refCounts);
        #$refMean  = 221.148738243302;
        #$refStdev = 36.8249415102796;
        print STDERR "$reference avg +/- sd = $refMean +/- $refStdev\n\n";
        solveCopyNumberHMM();

    # use an externally trained CN model, e.g., from svWGS
    } elsif($method eq 'bed6'){
        print STDERR "training model from external bed file CN assessment\n";
        openOutputStream($mdlFile, \$mdlH, $FALSE, $FALSE, "bgzip -c");
        openInputStream($globs[0], \my $inH, undef, undef, $globs[0] =~ m/\.gz$/);
        while (my $line = <$inH>){
            chomp $line;
            my @line = split("\t", $line);
            print $mdlH join("\t", @line[0..5]), "\n";
        }
        close($inH);
        
    # train the bin model based simply on sex of mammalian sample
    } else {
        print STDERR "training model from simple sex-based chromosome ploidy\n";
        openOutputStream($mdlFile, \$mdlH, $FALSE, $FALSE, "bgzip -c");
        foreach my $chrom(@chroms){
            my $CN = defined $chrCN{$chrom} ? $chrCN{$chrom} :  2; # autosome = CN2
            $CN or next; # model only has target chromosomes        
            print $mdlH join("\t", $chrom, 0, 275e6, '.', $CN, '+'), "\n"; # bigger than any chromosome
        }  
    }
    
    # index the bin model file
    print STDERR "$utility train: indexing model file(s)\n";
    closeHandles($mdlH);
    system("tabix -p bed $mdlFile");
    -e $binFile and system("tabix -p bed -f $binFile");

    # calculate the probe (not bin) GC content for normalization
    print STDERR "extracting probe positions\n";
    $inputDir = $arrayDir;
    openOutputStream($probeBedFile, \$prbH, $TRUE, $FALSE, "sort -k1,1 -k2,2n -S 2G");
    parseIlluminaProbes(\&writeProbeSpanBinFile);
    closeHandles($prbH);
    print STDERR join("\t", commify($nPrbTotal), "probes recovered from index array"), "\n";
    print STDERR "calculating probe GC values\n";
    my $gcFile = getOutFile('probes', $arrayType, 'gc.bed.gz');
    system("bedtools nuc -fi $genomeFasta -bed $probeBedFile | cut -f 4,6 | gzip -c > $gcFile");
}

# calculate fixed-width bin coverage
sub getBinCoverage_train {
    my @bamFiles;
    foreach my $glob(@globs){ push @bamFiles, glob($glob) }
    my $bamFiles = join(" ", @bamFiles);
    
    # proceed one requested chromosome at a time
    foreach my $chrom_(@chroms){
        $chrom = $chrom_;
        $binE = $binSize;
        $thsBinC = 0;
        $nxtBinC = 0;
        print STDERR "\t", $chrom, "\n";
        
        # thread through all well-mapped pairs
        my $f = 32; # 2 (data set when I wrote this didn't set bit 2, proper pair)
        my $F = 4 + 8 + 16 + 256 + 512 + 1024 + 2048;
        my $command1 = @bamFiles == 1 ?
                       "samtools view -u $bamFiles $chrom" :
                       "samtools merge -R $chrom -u - $bamFiles";
        open my $bamH, "-|", "$command1 | samtools view -q $minQual -f $f -F $F -"
            or die "$error: could not open bam merge stream, chrom $chrom\n";
        while (my $aln = <$bamH>) {
            my @bam = split("\t", $aln, 10);
            $bam[6] eq '=' or next;         
            $bam[8] <= 0 and next;
            $bam[8] > $maxTLen and next; 
            my $runS = $bam[3];
            my $runE = $runS + $bam[8] - 1;
            while($runS > $binE){ commitBin() }
            if($runE > $binE){
                my $runC = 1 / $bam[8]; 
                $thsBinC += $runC * ($binE - $runS + 1);
                $nxtBinC += $runC * ($runE - $binE);
            } else {
                $thsBinC += 1;
            }            
        }  
        closeHandles($bamH); 
        commitBin();
    }
}
sub commitBin {
    my $binS = $binE - $binSize + 1;
    if($thsBinC){
        print $cntH join("\t",
            $chrom, $binS, $binE, int($thsBinC/$prec + 0.5) * $prec
        ), "\n";
        if($chrom eq $refChrom and
           $binS >= $refStart and
           $binE <= $refEnd){
            push @refCounts, $thsBinC;
        }
    }
    $binE += $binSize;
    $thsBinC = $nxtBinC;
    $nxtBinC = 0;
}

# solve for the modal copy number across all samples
sub solveCopyNumberHMM {
    my ($mean1, $stdev1) = adjustSizeStats(1);
    $ENV{LIB_DIR}     = $libDir;
    $ENV{R_LIB_DIR}   = $rLibDir;
    $ENV{MODEL_NAME}  = $modelName;
    $ENV{DATAFILE}    = $cntFile;
    $ENV{MEAN_1}      = $mean1;
    $ENV{STDEV_1}     = $stdev1; 
    $ENV{MAX_CN}      = $maxCN;
    $ENV{EPROB_FILE}  = "$tmpDir/msvtools.$modelName.emiss.prob.".int(rand(1e6)).".txt";
    $ENV{TRAIN_RDATA} = getOutFile('model', $modelName, 'RData');
    $ENV{HMM_FILE}    = $mdlFile;
    $ENV{BIN_FILE}    = $binFile;
    $ENV{PLOT_DIR}    = "$outputDir/plots";
    $ENV{BIN_SIZE}    = $binSize;
    -d $ENV{PLOT_DIR} or mkdir $ENV{PLOT_DIR};
    system("Rscript $libDir/train_CN_HMM.R") == 0 or die "error solving CN HMM\n";  
}
sub adjustSizeStats {
    my ($CN) = @_;
    return ($refMean  / $refPloidy * $CN, 
            $refStdev / $refPloidy * $CN); 
}

# ues bedtools nuc to get the percent GC for DNA surrounding every probe position
sub writeProbeSpanBinFile { 
    my ($f, $hdr) = @_;
    $nPrbTotal++;
    my $start = $$f[$$hdr{POS}] - 101;
    print $prbH join("\t", 
        $$f[$$hdr{CHROM}],
        $start < 0 ? 0 : $start,
        $$f[$$hdr{POS}] + 100,
        $$f[$$hdr{PRB_NAME}]
    ), "\n";
}

1;
    
#=========================================================================
# SAM/BAM
#-------------------------------------------------------------------------
#  1	QNAME	Query template/pair NAME
#  2	FLAG	bitwise FLAG   NOTE: always 0 or 16 when using tophat
#  3	RNAME	Reference sequence NAME
#  4	POS	    1-based leftmost POSition/coordinate of clipped sequence
#  5	MAPQ	MAPping Quality (Phred-scaled)
#  6	CIGAR	extended CIGAR string
#  7	MRNM	Mate Reference sequence NaMe (?=? if same as RNAME)
#  8	MPOS	1-based Mate POSistion
#  9	TLEN	inferred Template LENgth (insert size)
#  10	SEQ	    query SEQuence on the same strand as the reference
#  11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
#  12+	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE
#-------------------------------------------------------------------------
# flag bits
#-------------------------------------------------------------------------
#  1    0x1    template having multiple segments in sequencing
#  2    0x2    each segment properly aligned according to the aligner
#  4    0x4    segment unmapped
#  8    0x8    next segment in the template unmapped
#  16   0x10   SEQ being reverse complemented
#  32   0x20   SEQ of the next segment in the template being reverse complemented
#  64   0x40   the first segment in the template
#  128  0x80   the last segment in the template
#  256  0x100  secondary alignment
#  512  0x200  not passing filters, such as platform/vendor quality controls
#  1024 0x400  PCR or optical duplicate
#  2048 0x800  supplementary alignment
#-------------------------------------------------------------------------
# CIGAR operations
#-------------------------------------------------------------------------
#  M 0 alignment match (can be a sequence match or mismatch)
#  I 1 insertion to the reference
#  D 2 deletion from the reference
#  N 3 skipped region from the reference
#  S 4 soft clipping (clipped sequences present in SEQ)
#  H 5 hard clipping (clipped sequences NOT present in SEQ)
#  P 6 padding (silent deletion from padded reference)
#  = 7 sequence match
#  X 8 sequence mismatch
#=========================================================================


# $ bedtools nuc

# Tool:    bedtools nuc (aka nucBed)
# Version: v2.30.0
# Summary: Profiles the nucleotide content of intervals in a fasta file.

# Usage:   bedtools nuc [OPTIONS] -fi <fasta> -bed <bed/gff/vcf>

# Options: 
#         -fi     Input FASTA file

#         -bed    BED/GFF/VCF file of ranges to extract from -fi

#         -s      Profile the sequence according to strand.

#         -seq    Print the extracted sequence

#         -pattern        Report the number of times a user-defined sequence
#                         is observed (case-sensitive).

#         -C      Ignore case when matching -pattern. By defaulty, case matters.

#         -fullHeader     Use full fasta header.
#                 - By default, only the word before the first space or tab is used.

# Output format: 
#         The following information will be reported after each BED entry:
#             1) %AT content
#             2) %GC content
#             3) Number of As observed
#             4) Number of Cs observed
#             5) Number of Gs observed
#             6) Number of Ts observed
#             7) Number of Ns observed
#             8) Number of other bases observed
#             9) The length of the explored sequence/interval.
#             10) The seq. extracted from the FASTA file. (opt., if -seq is used)
#             11) The number of times a user's pattern was observed.
#                 (opt., if -pattern is used.)


# initialize the copy number HMM for genome based on reference region
#sub initializeCopyNumberHMM {
#    
#    # set the boundaries for the bin-count HMM
#    my $nSD   = 3.5;  # helps determine HMM boundary limits
#    my $nBins = 100;  # number of divisions, i.e. observation states, in the HMM (actually get one more)
#    $maxBinI = $nBins + 2;
#    my ($mean1, $stdev1) = adjustSizeStats(1);       # CN=1, smallest bin count
#    my ($meanX, $stdevX) = adjustSizeStats($maxCN);  # largest bin count
#    my $minCount = $mean1 - $nSD * $stdev1; # largest and smallest bin counts to allow in HMM
#    my $maxCount = $meanX + $nSD * $stdevX;    
#    $binWidth = int(($maxCount - $minCount) / $nBins + 0.5); # width of each bin
#    $maxBinCount = getBinCount($maxCount); # rounded bin values at the limit bins in unit varBin length
#    $minBinCount = getBinCount($minCount);
#    
#    # generate a set of copy number emission probabilities for the genome based on reference
#    $ENV{eProbFile_CN} = "$tmpDir/msvtools.$modelName.emiss.prob.".int(rand(1e6)).".txt";
#    $ENV{refMean}     = $refMean;
#    $ENV{refStdev}    = $refStdev;
#    $ENV{refPloidy}   = $refPloidy;    
#    $ENV{binWidth}    = $binWidth;
#    $ENV{maxBinCount} = $maxBinCount;
#    $ENV{minBinCount} = $minBinCount;
#    $ENV{maxCN}       = $maxCN;
#    system("Rscript $libDir/set_CN_emissions.R") == 0
#        or die "error setting CN emission probabilities\n";      
#}
#
#sub getBinCount {
#    my ($binCount) = @_;
#    return int($binCount / $binWidth + 0.5) * $binWidth; 
#}

    ## convert bin data to observation indices, correlated to bin indices
    #my $tmpFile = "$tmpDir/msvtools.$modelName.CN.segment.".int(rand(1e6)).".data";  
    #open(my $tmpH, ">", $tmpFile) or die "could not open $tmpFile for writing: $!\n";
    #while(my $line = <$cntH>){
    #    chomp $line;
    #    my ($chrom, $start, $end, $count) = split("\t", $line);
    #    print $tmpH join("\t", getBinI($count), $chrom, $start, $end, $count), "\n"; 
    #}
    #close $tmpH;
    
    
    ## set null hypotheses for different copy number states
    #my @nullH;
    #foreach my $CN(1..$maxCN){
    #    push @nullH, [adjustSizeStats($CN)]; # mean, stdev
    #}
    #
    ## convert resulting CN state indices to copy number
    #open(my $hmmH, "-|", "cat $tmpFile | ".
    #                     "segment -e $ENV{eProbFile_CN} -z 0.1 -p 0.99 | ".
    #                     "groupBy -g 1,2 -c 3,4,5,5,5 -o min,max,count,mean,stdev") 
    #    or die "could not open CN segment stream: $!\n";  
    #while (my $line = <$hmmH>) {
    #    chomp $line;
    #    my ($CN, $chrom, $start, $end, $nBins, $meanCount, $sdCount) = split("\t", $line);
    #    my @ts;
    #    foreach my $nullH(@nullH){
    #        my $Z = ($meanCount - $$nullH[0]) / ($$nullH[1] / sqrt($nBins));
    #        my $s = $sdCount / $$nullH[1];
    #        push @ts, $Z / $s;
    #    }
    #    $CN > $maxCN and $CN = $maxCN;
    #    print $mdlH join("\t", $chrom, $start, $end, '.', $CN, '+',
    #                           $nBins, $meanCount, $sdCount,
    #                           join(",", @ts) ), "\n";
    #}
    #close $hmmH;
    #unlink $tmpFile;
    
    
    #sub getBinI {
#    my ($binCount) = @_;
#    my $countBin = getBinCount($binCount);
#    if ($countBin < $minBinCount ) {
#        return 0;
#    } elsif($countBin > $maxBinCount){
#        return $maxBinI;
#    } else {
#        return ($countBin - $minBinCount) / $binWidth + 1;
#    }
#}