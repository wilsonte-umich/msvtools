use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE);
require "$libDir/Illumina.pl";
our ($inputDir, $arrayFmt, $project, $nameCol,
     $sample,
     $outputDir, $plotDir, $maxMem);
my ($prbH, $lrrH, $bafH);
my (%LRR, %BAF);
my %lrr = (scalar=>20,  min=>-1, max=>1);
my %baf = (scalar=>40,  min=>0,  max=>1);

# manage options
sub setOptions_index {
    setOptionValue(\$inputDir,     'input-dir');
    setOptionValue(\$arrayFmt,     'array-format', 'Illumina');    
    setOptionValue(\$project,      'project');
    setOptionValue(\$nameCol,      'name-column',  'DNA_ID');
    setOptionValue(\$sample,       'sample');
    setOptionValue(\$outputDir,    'output-dir');
    #setOptionValue(\$plotDir,      'plot-dir');
    setOptionValue(\$maxMem,       'max-mem',       1000000000);
    #-------------------------------------
    -d $inputDir or die "$error: directory not found: $inputDir\n";
    $maxMem =~ m|\D| and die "$error: max-mem must be an integer number of bytes\n";
}

# main execution block
sub msvtools_index {

    # initialize
    setOptions_index();    
    print STDERR "$utility index: " . getTime(), "\n";
    
    # open outputs
    my $prbFile = getOutFile('probes', $sample, 'bgz');
    my $lrrFile = getOutFile('lrr',    $sample, 'txt');    
    my $bafFile = getOutFile('baf',    $sample, 'txt'); 
    openOutputStream($prbFile, \$prbH, $FALSE, $FALSE,
                     "sort -S $maxMem"."b -k1,1 -k2,2n | bgzip -c");
    openOutputStream($lrrFile, \$lrrH);    
    openOutputStream($bafFile, \$bafH); 
 
    # do the work
    if($arrayFmt eq 'Illumina'){
        indexIllumina();
    } else {
        die "unknown array format: $arrayFmt\n";
    }
    
    # print the histograms
    print STDERR "$utility index: print histogram data\n";
    printIndexHistogram($lrrH, \%lrr, \%LRR);
    printIndexHistogram($bafH, \%baf, \%BAF);

    # index the array file
    print STDERR "$utility index: indexing probes file\n";
    closeHandles($prbH, $lrrH, $bafH);
    system("tabix -s 1 -b 2 -e 2 $prbFile");
}

sub printIndexedLine {
    my ($chrom, $pos,
        $qual,
        $x, $y, $r, $theta,
        $lrr, $baf) = @_;
    print $prbH join("\t",
        $chrom, $pos,
        $qual,
        $x, $y, $r, $theta,
        $lrr, $baf
    ), "\n";
    $LRR{binArrayData($lrr, $lrr{scalar})}++;
    $BAF{binArrayData($baf, $baf{scalar})}++;
}
sub binArrayData {
    my ($val, $scalar) = @_;
    $scalar or $scalar = 1000;
    if($val >= 0){
        int($val * $scalar + 0.5);
    } else {
        -(int(-$val * $scalar + 0.5));
    }
}

sub printIndexHistogram {
    my ($outH, $prm, $dat) = @_;
    my $nPrb = 0;
    my $min = binArrayData($$prm{min}, $$prm{scalar});
    my $max = binArrayData($$prm{max}, $$prm{scalar});
    for (my $bin = $min; $bin <= $max; $bin++){
        $$dat{$bin} = ($$dat{$bin} || 0);
        $nPrb += $$dat{$bin};
    }
    for (my $bin = $min; $bin <= $max; $bin++){
        print $outH join("\t", $bin/$$prm{scalar}, roundCount($$dat{$bin} / $nPrb)), "\n";
    } 
}

1;
