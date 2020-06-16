use strict;
use warnings;

#========================================================================
# 'common.pl' has miscellaneous subs common to many commands and utilities
#========================================================================

#========================================================================
# define variables
#------------------------------------------------------------------------
use vars qw($utility $command $libDir @args $error
            $slurp $genomeFasta
            %chroms %prbCol %svCol);
our ($TRUE, $FALSE) = (1, undef);
our ($inH, $tmpH, $outH, $rptH);
our ($inputFile, $region, $ouputFile,
	 $inputDir, $tmpDir, $outputDir, $plotDir);
#========================================================================

#========================================================================
# array manipulations
#------------------------------------------------------------------------
sub filterArrayLine {
    my ($line, $hdr, $sub) = @_;
    chomp $line; 
    $line =~ s/\r//g;
    my @f = split("\t", $line);
    $f[$$hdr{CHROM}] =~ m/^chr/ or $f[$$hdr{CHROM}] = "chr$f[$$hdr{CHROM}]";    
    $chroms{$f[$$hdr{CHROM}]} or return;               
    $f[$$hdr{POS}] =~ m/^\d+$/ or return;    
    if($f[$$hdr{LRR}] eq 'NaN' or $f[$$hdr{BAF}] eq 'NaN'){
        $f[$$hdr{LRR}] = $f[$$hdr{BAF}]= 'NA'; # keep failed probes in an R-friendly manner
    }
    &$sub(\@f, $hdr); 
}
sub segmentArray {
    my ($name, $datFile, $plotsDir, $rData,
        $extZyg, $lrrLim, $badPrbFreq, $transProb, $preservation) = @_;
    # make sure refine has segmentation defaults; segment provides user-exposed options
    $extZyg       or $extZyg       = 0.85; 
    $lrrLim       or $lrrLim       = "-1.5,1.0";
    $badPrbFreq   or $badPrbFreq   = 1e-3;
    $transProb    or $transProb    = 1e-4;
    $preservation or $preservation = "0.99,0.95,0.9,0.8";
    # set parameters
    my $hmmFile = getOutFile('HMM',      $name, 'RData');
    my $prbFile = getOutFile('probes',   $name, 'bed.bgz');
    my $segFile = getOutFile('segments', $name, 'bed.bgz');
    my $svsFile = getOutFile('SVs',      $name, 'bed');
    $ENV{LIB_DIR}      = $libDir;
    $ENV{MODEL_NAME}   = $name;
    $ENV{DATAFILE}     = $datFile;
    $ENV{PLOT_DIR}     = "$plotsDir/$name";
    $ENV{HMM_FILE}     = $hmmFile;
    $ENV{PROBES_FILE}  = $prbFile;
    $ENV{SEGMENTS_FILE}= $segFile;
    $ENV{SVS_FILE}     = $svsFile;
    $ENV{TRAIN_RDATA}  = $rData;
    $ENV{PROBE_COL}    = join(",", sort {$prbCol{$a} <=> $prbCol{$b}} keys %prbCol);
    $ENV{SV_COL}       = join(",", sort {$svCol{$a} <=> $svCol{$b}} keys %svCol);
    $ENV{extreme_zyg}     = $extZyg;
    $ENV{lrr_lim}         = $lrrLim;
    $ENV{BAD_PROBE_FREQ}  = $badPrbFreq;
    $ENV{PERSISTENCE}     = 1 - $transProb;
    $ENV{PRESERVATION}    = $preservation;
    mkdir $plotsDir;      
    mkdir $ENV{PLOT_DIR};     
    system("Rscript $libDir/segment_HMM.R") and die "$error: segment_HMM.R failed\n";
    #unlink $datFile;

    # index the model file
    print STDERR "$utility: indexing model files\n";
    -e $prbFile and system("tabix -f -p bed $prbFile");
    -e $segFile and system("tabix -f -p bed $segFile");
    if(-e $svsFile){
        system("cat $svsFile | bgzip -c > $svsFile.bgz");
        system("tabix -f -p bed $svsFile.bgz");
    }
}
#========================================================================

#========================================================================
# sequence manipulation
#------------------------------------------------------------------------
sub revComp {
    my $seq = reverse(uc($_[0]));
    $seq =~ tr|ACGT|TGCA|;
    return $seq;
}
sub getRegSeqs {
	@_ or return [];
	open(my $faH, "-|", "samtools faidx $genomeFasta ".join(" ", @_))
		or die "could not open samtools faidx stream: $!\n";
	my $i = -1;
	my @regSeqs;
	while (my $line = <$faH>) {
		chomp $line;
		if ($line =~ m|^>|) { $i++ } else { $regSeqs[$i] .= uc($line) }
	}
	close $faH;
	return \@regSeqs;
}
#========================================================================

#========================================================================
# miscellaneous
#------------------------------------------------------------------------
sub reportCount {
    my ($n, $desc) = @_;
    print STDERR join("\t", sprintf("$utility $command: %13s", commify($n)), $desc), "\n";
}
sub commify {
    local $_  = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}
sub min {
    my ($v1, $v2) = @_;
    $v1 <= $v2 ? $v1 : $v2;
}
sub max {
    my ($v1, $v2) = @_;
    $v1 >= $v2 ? $v1 : $v2;
}
sub median {
    my (@data) = sort {$a <=> $b} @_;
    my $i = @data / 2;
    my $upper = $data[$i];
    @data % 2 and return $upper;
    my $lower = $data[$i - 1];    
    return($lower + ($upper - $lower) / 2);
}
sub mean{
    my (@data) = @_;
    @data == 0 and die "no values sent to mean\n";
    my $sum = 0;
    foreach (@data) { $sum += $_ }
    return $sum / @data;
}
sub stdev{
    my (@data) = @_;
    @data <= 1 and return (@data, 0);
    my $mean = mean(@data);
    my $sqSum = 0;
    foreach(@data) { $sqSum += ($mean-$_) ** 2 }
    return ($mean, ($sqSum / (@data-1)) ** 0.5);
}
sub stdev_center{ # exclude high and low values and return mean/stdev of remainder
    my (@data) = @_;
    @data >= 3 or return (0, 0);    
    @data = sort {$a <=> $b} @data;
    return stdev(@data[1..$#data-1]); # or could apply more complex outlier exclusion
}
sub roundCount {
    my ($val, $scalar) = @_;
    $scalar or $scalar = 1000;
    if($val >= 0){
        int($val * $scalar + 0.5) / $scalar;
    } else {
        -(int(-$val * $scalar + 0.5) / $scalar);
    }
}
#========================================================================

1;
