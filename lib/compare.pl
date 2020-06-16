use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            %prbCol %svCol @RLL_COL @RLL_SMP_COL %cmpCol);
require "$libDir/formats.pl";
our ($inputDir, $samples, $outputDir, $sizeThreshold, $groupName, $tmpDir);
our (@samples, %svs, $collFile);
my ($collH);
my $regionN = 0;
use constant {
    LARGE => 'large',
    SMALL => 'small'
};

# manage options
sub setOptions_compare {
    setOptionValue(\$inputDir,   'input-dir');
    setOptionValue(\$samples,    'samples');   
    setOptionValue(\$outputDir,  'output-dir');
    setOptionValue(\$sizeThreshold,  'size-threshold', 2000000);       
    setOptionValue(\$groupName,  'group-name');
    setOptionValue(\$tmpDir,  'tmp-dir', '/tmp');
    #-------------------------------------
    -d $inputDir  or die "$error: directory not found: $inputDir\n";   
    -d $outputDir or die "$error: directory not found: $outputDir\n";
    @samples = sort {$a cmp $b} split(",", $samples);
    -d $tmpDir or die "$error: directory not found: $tmpDir\n";
}
# main execution block
sub msvtools_compare {

    # initialize
    setOptions_compare();    
    print STDERR "$utility compare: " . getTime(), "\n";
    
    # write the samples file (used by the web interface)
    my $cmpFile = getOutFile('samples', $groupName, 'txt');
    open my $cmpH, ">", $cmpFile or die "could not open $cmpFile\n";
    print $cmpH join(" ", @samples), "\n";
    print $cmpH join(" ", generateBlindNames(@samples)), "\n";
    close $cmpH;
    
    # load input SV BED file for each sample
    map { loadBedFile($_) } @samples;

    # collapse the SV spans over all samples into query regions
    # attempt to avoid overly zealous fusions with size threshold
    foreach my $size(SMALL, LARGE){
        keys %{$svs{$size}} or next;
        print STDERR "\n$size CNVs\n\n";    
        $collFile = getTmpFile('collapsed', 'regions.'.$size, 'tmp.gz');
        collapseBED($size);
            
        # use R to calculate the log likelihoods for all SVs for all sample pairs
        my $gz = "$groupName.$size";
        my $cmpFile = getOutFile('SVs', $gz, 'bed');
        $ENV{LIB_DIR}      = $libDir;
        $ENV{DATAFILE}     = $collFile;
        $ENV{PLOT_DIR}     = "$outputDir/plots";
        $ENV{SIZE_THRESHOLD} = $sizeThreshold;
        $ENV{SIZE}         = $size;
        $ENV{GROUP_NAME}   = $gz;
        $ENV{SAMPLES_DIR}  = $inputDir;
        $ENV{SAMPLES}      = $samples;
        $ENV{RLL_COL}      = join(",", @RLL_COL);
        $ENV{RLL_SMP_COL}  = join(",", @RLL_SMP_COL);
        $ENV{CMP_COL}      = join(",", sort {$cmpCol{$a} <=> $cmpCol{$b}} keys %cmpCol);
        $ENV{CMP_FILE}     = $cmpFile;
        mkdir $ENV{PLOT_DIR};
        system("Rscript $libDir/compare_SVs.R") and die "$error: compare_SVs.R failed\n";
        unlink $collFile;        
    }
}

sub generateBlindNames {
    my (@smp) = @_;
    my @alpha = 'A'..'Z';
    my @nmric = 0..9;
    my @blind = map {
        join("",
             $alpha[int(rand(26))],
             $alpha[int(rand(26))],
             '-',
             $nmric[int(rand(10))],
             $nmric[int(rand(10))]
        );
    } @smp; 
}

sub loadBedFile {
    my ($sample) = @_;
    my $bedFile = getInFile('segment', 'SVs', $sample, 'bed');
    openInputStream($bedFile, \my $bedH);
    my $cI  = $svCol{CHROM};
    my $sI  = $svCol{START}; # so that regions will include flanking probes
    my $eI  = $svCol{END};   # used to determine SV uniqueness
    my $fsI = $svCol{FLANK_START}; # so that regions will include flanking probes
    my $feI = $svCol{FLANK_END};   # used to determine SV uniqueness
    while (my $line = <$bedH>){
        $line =~ m|^\s*#| and next;  # ignore comment lines
        chomp $line;
        $line =~ s/\r//g;
        my @sv = split("\t", $line);
        my $size = ($sv[$eI] - $sv[$sI] > $sizeThreshold) ? LARGE : SMALL;
        push @{$svs{$size}{$sv[$cI]}{$sv[$fsI]}{$sv[$feI]}}, \@sv;
    }
    closeHandles($bedH);
}

sub collapseBED {
    my ($size) = @_;
    openOutputStream($collFile, \$collH, $TRUE);
    print $collH join("\t",
        qw(CHROM START END N_SVS N_SAMPLES SAMPLES REGION_ID)), "\n";
    foreach my $chrom(sort {$a cmp $b} keys %{$svs{$size}}){
        my ($regionStart, @svs);
        my $regionEnd = 0;        
        foreach my $start(sort {$a <=> $b} keys %{$svs{$size}{$chrom}}){
            if(!$regionEnd){
                $regionStart = $start;
            } elsif($start > $regionEnd){
                printSVRegion($chrom, $regionStart, $regionEnd, \@svs);
                ($regionStart, $regionEnd) = ($start, 0);
                @svs = ();
            }
            foreach my $end(sort {$a <=> $b} keys %{$svs{$size}{$chrom}{$start}}){
                $regionEnd >= $end or $regionEnd = $end;
                foreach my $sv(@{$svs{$size}{$chrom}{$start}{$end}}){
                    push @svs, $sv;
                }
            }
        }
        printSVRegion($chrom, $regionStart, $regionEnd, \@svs);
    }
    closeHandles($collH);
}
sub printSVRegion {
    my ($chrom, $start, $end, $svs) = @_;
    my (%svSamples);
    foreach my $sv(@$svs){ $svSamples{$$sv[$svCol{SAMPLE}]} = 1 }
    my $nSVs = @$svs;    
    my $nSamples = keys %svSamples;
    my $samples = join(",", sort {$a cmp $b} keys %svSamples);
    $regionN++;
    print $collH join("\t",
        $chrom, $start, $end, $nSVs, $nSamples, $samples, $regionN), "\n";   
}

1;
