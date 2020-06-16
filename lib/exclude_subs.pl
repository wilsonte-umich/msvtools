use strict;
use warnings;

#========================================================================
# 'exclude.pl' supplies the exlude function, adapted from bed_util crossing
#========================================================================

#========================================================================
# define variables
#------------------------------------------------------------------------
use vars qw($utility $command);
my $error;
my (%exclRegions, $exclBinSize);
my ($collRefBases, $maxCollRefLength, $maxCollRefEnd, $nCollRefFeatures) = (0, 0, 0, 0);
#========================================================================


#========================================================================
# collect and process the excluded, i.e. reference, features
#------------------------------------------------------------------------
sub initializeExclude {
    my ($excludefile) = @_;
    $excludefile or return;
    $error = "$utility $command error";
    (%exclRegions, $exclBinSize) = ();
    ($collRefBases, $maxCollRefLength, $maxCollRefEnd, $nCollRefFeatures) = (0, 0, 0, 0);
    loadBedFile($excludefile, \%exclRegions);
    analzyeExcludeChroms();   
    my $aveCollRefLength = int($collRefBases / $nCollRefFeatures);    
    $exclBinSize = int($maxCollRefEnd / 100);  # nominally split largest chromosome into 100 bins for faster searching (but higher memory use)
    my $minRefBinSize = $aveCollRefLength * 2;
    $exclBinSize = $exclBinSize > $minRefBinSize ?  $exclBinSize : $minRefBinSize;   # but use larger bins if reference features are especially long
    binRefFeatures();  
}
#========================================================================


#========================================================================
# check whether a feature should be excluded
#------------------------------------------------------------------------
sub checkExclude {  
    my ($chrom, $testStart, $testEnd) = @_;
    $exclRegions{$chrom} or return;
    #($testEnd - $testStart) > 0 or die "$error: malformed request to checkExclude: $chrom, $testStart, $testEnd\n";
    if ($testEnd - 1 < $testStart) {
        ($testStart, $testEnd) = ($testEnd - 1, $testStart + 1);
    }
    my ($startBin, $endBin) = getCrossedBins($testStart, $testEnd);
    for (my $bin = $startBin; $bin <= $endBin; $bin += $exclBinSize){  # check all bins crossed by query feature
        $exclRegions{$chrom}{$bin} or next;
        foreach my $refRegion(@{$exclRegions{$chrom}{$bin}}){
            $$refRegion[0] <= $testEnd and $$refRegion[1] >= $testStart and return $$refRegion[2];         
        }
    }
    return undef;
}
#========================================================================


#========================================================================
# worker subs
#------------------------------------------------------------------------
sub loadBedFile {
    my ($bedFile, $featuresHash) = @_;
    my $inH;
    if($bedFile =~ m/\.gz$/ or $bedFile =~ m/\.bgz$/){
        open $inH, "-|", "gunzip -c $bedFile" or die "$error: could not open $bedFile: $!\n";
    } else {
        open $inH, "<", $bedFile or die "$error: could not open $bedFile: $!\n";
    }
    while (my $line = <$inH>){
        $line =~ m|^\s*#| and next;  # ignore comment lines
        chomp $line;
        $line =~ s/\r//g;
        my ($chrom, $start, $end, $name, $score) = split("\t", $line, 6);        
        $chrom or next;  # ignore blank lines                
        parseInt(\$start);  # start and end must be present and integer numbers
        parseInt(\$end);               
        push @{$$featuresHash{$chrom}}, [$start, $end, $score];
    }  
    close $inH;
}
sub parseInt {  # validate and uncommify integer values
    my ($int) = @_;  # int passed as reference
    defined $$int or die "$error: missing start or end in exclude-file\n";
    $$int =~ s/,//g;
    $$int =~ m|\D| and die "$error: invalid integer in exclude-file\n";
}
#------------------------------------------------------------------------
sub analzyeExcludeChroms {
    foreach my $chrom(keys %exclRegions){
        foreach my $refRegion(@{$exclRegions{$chrom}}){  # collect reference feature stats
            my ($start, $end) = @$refRegion;
            my $featureLength = $end - $start;
            $collRefBases += $featureLength;
            $maxCollRefLength >= $featureLength or $maxCollRefLength = $featureLength;
            $maxCollRefEnd >= $end or $maxCollRefEnd = $end;
            $nCollRefFeatures++;
        }     
    }
}
#------------------------------------------------------------------------
sub binRefFeatures {  # split reference features into bins on each chromosome strand to speed up crossing search
    my %binned;
    foreach my $chrom(keys %exclRegions){
        foreach my $refRegion(@{$exclRegions{$chrom}}){
            my ($startBin, $endBin) = getCrossedBins(@$refRegion);
            for (my $bin = $startBin; $bin <= $endBin; $bin += $exclBinSize){
                push @{$binned{$chrom}{$bin}}, $refRegion;
            }
        }
    } 
    %exclRegions = %binned;
}
sub getCrossedBins {  # the lowest and highest bins crossed by a region
    my ($start, $end) = @_;
    return (int($start     / $exclBinSize) * $exclBinSize,
            int(($end - 1) / $exclBinSize) * $exclBinSize);  # end bins are converted to 0-referenced start coordinates
}
#========================================================================

1;
