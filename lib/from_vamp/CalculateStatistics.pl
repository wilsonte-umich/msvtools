#!/usr/bin/perl -w
use strict;
use warnings;

#use Algorithm::CurveFit;

#########################################################################
#CalculateStatistics.pl find the main Convergent Pairs peak,
#describes that peak by curve-fitting a Gaussian through non-linear regression,
#and deposits the result Normal statistics into the Stats table.
#########################################################################

use vars(qw(%param %types %fields %refSeqs));

sub calcStats { #provides direct call to calcStats separate from extractPairs
    my ($sample) = @_;
    my $pairsTable = getTableName('Pairs', $sample);
    my $histogramTable = getTableName('h', $pairsTable);
    my $statsTable = newTable('Stats', $sample);   
    calculateNormalStats($histogramTable, $statsTable);
}

sub calculateNormalStats{
    my ($histogramTable, $statsTable) = @_;
    #get initial stat estimates from fairly broad range of fragment sizes
    my ($initMin, $initMax) = (100, 1000);
    $param{isCircles} and ($initMin, $initMax) = (1000, 10000);
    my ($mean, $stDev, $amp, $mode) = getNormalStats($histogramTable, $initMin, $initMax);
    
    # my $SD3 = 3 * $stDev;
    # #get final values by recalculating based on 1st round estimates
    # ($mean, $stDev, $amp, $mode) = getNormalStats($histogramTable, $mean - $SD3, $mean + $SD3, $amp, $mean, $stDev);
    # #round to integer values
    # $mean = int($mean + 0.5); 
    # $stDev = int($stDev + 0.5);
    # $amp = int($amp + 0.5);
    # $mode = int($mode + 0.5);
    
    #calculate derivative Normal peak values
    my $SD3 = 3 * $stDev;
    my $minNormal = $mean - $SD3;
    my $maxNormal = $mean + $SD3;
    my $minInsertion = $mode - $minNormal;
    my $minDeletion = $maxNormal - $mode;
    #save and report results
    updateStat($statsTable, 'meanNormal', $mean);
    updateStat($statsTable, 'stDevNormal', $stDev);
    updateStat($statsTable, 'minNormal', $minNormal);
    updateStat($statsTable, 'maxNormal', $maxNormal);
    updateStat($statsTable, 'ampNormal', $amp);
    updateStat($statsTable, 'modeNormal', $mode);
    updateStat($statsTable, 'minFragSize', 250); #default minFragSize set here
    $param{isCircles} and updateStat($statsTable, 'maxReverseNormal', 750); #default maxReverseNormal set here
    updateStat($statsTable, 'minInsertion', $minInsertion);
    updateStat($statsTable, 'minDeletion', $minDeletion);
    status("  Normal peak statistics:\n");
    status("      mean, stDev = $mean +/- $stDev bp\n");
    status("      peak tip at $mode bp\n");
    status("      peak ranges from $minNormal to $maxNormal bp\n");
    status("      minimum size of detectable Insertion = $minInsertion bp\n");
    status("      minimum size of detectable Deletion = $minDeletion bp\n");
    status("      NOTE!  use Graphs.pl to visualize the accuracy of and adjust these values!\n");
}

sub getNormalStats{
    my ($histogramTable, $lowerLimit, $upperLimit, $aGuess, $mGuess, $sGuess) = @_;
    my $sql = "SELECT X, Y
                FROM $histogramTable
                WHERE X >= $lowerLimit AND X <= $upperLimit
                  AND SERIES = $types{Pairs}{Convergent}
                ORDER BY X ";
    return fitGaussian($sql, $aGuess, $mGuess, $sGuess, $histogramTable);
}

sub fitGaussian{
    my ($sql, $aGuess, $mGuess, $sGuess, $histogramTable) = @_;
    #find the tip of the peak, assign to maxN and mode
    # runSQL("SELECT Max(Y) FROM ($sql)");
    # my ($maxY) = fetchRowArray();
    # runSQL("SELECT Avg(X)
            # FROM ($sql)
            # WHERE Y > (0.90 * $maxY)");
    # my ($mode) = fetchRowArray();
    # #establish Gaussian value guesses if none provided
    # defined $aGuess or $aGuess = $maxY;
    # defined $mGuess or $mGuess = $mode;
    # defined $sGuess or $sGuess = $mGuess/10;
  
    #use R instead of Algorithm::CurveFit, which is unreliable
    runSQL($sql, \my($xValue, $yValue));
    my $outFile = "$histogramTable.csv";
    open my $outH, ">", $outFile or die "could not open $outFile\n";
    print $outH "SIZE,COUNT\n";
    while (fetchRow()){ print $outH "$xValue,$yValue\n" }
    close $outH;
    my $rCommand = "Rscript $param{vampPath}/bin/CalculateStatistics.R $outFile";
    my $fitResults = qx/$rCommand/;
    my ($meanExpected, $stdDevExpected, $minExpected, $maxExpected, $amplitude) = split(" ", $fitResults);
    return ($meanExpected, $stdDevExpected, $amplitude, $meanExpected); 

#    #parse histogram data into required parallel x and y arrays
#    runSQL($sql, \my($xValue, $yValue));
#    my @xData;
#    my @yData;
#    while (fetchRow()){
#        push @xData, $xValue;
#        push @yData, $yValue;
#    }
#    #run curve_fit
#    my $formula = 'a*(exp(-((x-m)^2)/(2*(s^2))))';  #Gaussian equation
#    my $variable = 'x';   
#    my @parameters = (
#        # Name    Guess   Accuracy
#        ['m',     $mGuess,     1],   #mean (mu)
#        ['s',     $sGuess,     1],    #standard deviation (sigma)
#        ['a',     $aGuess,    100] );   #amplitude
#    my $max_iter = 10; #prevent endless loops
#    my $square_residual = Algorithm::CurveFit->curve_fit(
#        formula            => $formula, 
#        params             => \@parameters,
#        variable           => $variable,
#        xdata              => \@xData,  
#        ydata              => \@yData,
#        maximum_iterations => $max_iter);
#    #return mean, stdev, amplitude, mode
#    return ($parameters[0][1], $parameters[1][1], $parameters[2][1], $mode); 
}

sub updateStat{
    my ($statsTable, $statName, $statValue) = @_;
    runSQL("DELETE FROM $statsTable WHERE STATNAME = '$statName'");
    runSQL("INSERT INTO $statsTable (STATNAME, STATVALUE) VALUES ('$statName', $statValue)");
}

sub getStatistics{
    my ($statsTable, $statsRef) = @_;
    runSQL("SELECT STATNAME, STATVALUE FROM $statsTable", \my($statName, $statValue));
    while(fetchRow()){$$statsRef{$statName} = $statValue}
}

sub calculateCoverageStats{
    my ($histogramTable, $mapType, $statsTable) = @_; #mapType = 'FMap' or 'RMap'
    #get max count from coverage histogram
    runSQL("SELECT Max(Y) FROM $histogramTable");
    my ($maxY) = fetchRowArray();
    #revise coverage histogram to the peak above any noise or extended tails
    runSQL("SELECT X, Y FROM $histogramTable WHERE Y > (0.05 * $maxY)", \my($coverage, $N));
    #calculate mean and stdev by standard statistics
    #NOTE:  will NOT be very accurate for low coverages subject to Poisson skewing near 0 !!
    my @coverages;
    while (fetchRow()){ for my $i(1..$N){ push @coverages, $coverage } }
    my $count = scalar @coverages;
    my $sum;
    foreach my $coverage(@coverages){$sum += $coverage}
    my $coverageMean = $sum / $count;
    my $sumDevSq;
    foreach my $coverage(@coverages){ $sumDevSq += abs($coverage - $coverageMean) ** 2 }
    my $coverageStDev = ($sumDevSq / $count) ** 0.5;
    #calculate coverage limits for Set and Disc finding
    my $SD3 = 3 * $coverageStDev;
    my $minCoverage = $coverageMean - $SD3;
    my $maxCoverage = $coverageMean + $SD3;
    #round to nearest integer
    $coverageMean = int($coverageMean + 0.5);
    $coverageStDev = int($coverageStDev + 0.5);
    $minCoverage = int($minCoverage + 0.5);
    $maxCoverage = int($maxCoverage + 0.5);
    #insist on coverage = 2 at minimum
    $minCoverage >= 2 or $minCoverage = 2;
    #save and report values
    updateStat($statsTable, "meanCoverage$mapType", $coverageMean);
    updateStat($statsTable, "stDevCoverage$mapType", $coverageStDev);
    updateStat($statsTable, "minCoverage$mapType", $minCoverage);
    updateStat($statsTable, "maxCoverage$mapType", $maxCoverage);
    status("  coverage statistics as determined from $mapType\n");
    status("      meanCoverage$mapType = $coverageMean\n");
    status("      stDevCoverage$mapType = $coverageStDev\n");
    status("      minCoverage$mapType = $minCoverage\n");
    status("      maxCoverage$mapType = $maxCoverage\n");
}

#coerceStats helps deal with the high failure rate of Algorithm::CurveFit
#it expects that Graphs.pl was used to set minNormal, maxNormal and modeNormal
#these boundary values help better define the main peak for Gaussain fitting
#coerceStats should be used when extractPairs fails, before parseFragments
sub coerceStats{
    my ($sample) = @_;
    status("calculating statistics using minNormal, modeNormal and maxNormal...\n");
    my $pairsTable = getTableName('Pairs', $sample);
    my $histogramTable = getTableName('h', $pairsTable);
    my $statsTable = getTableName('Stats', $sample);   
    getStatistics($statsTable, \my%stats);   
    my $lowerLimit = $stats{minNormal};
    my $upperLimit = $stats{maxNormal};
    my $aGuess = undef;
    my $mGuess = $stats{modeNormal};
    my $sGuess = int(($upperLimit - $lowerLimit)/6);
    getNormalStats($histogramTable, $lowerLimit, $upperLimit, $aGuess, $mGuess, $sGuess);
}

1;

