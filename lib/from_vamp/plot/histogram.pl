#!/usr/bin/perl
use strict;
use warnings;

#wholly generic functions for creating sample histogram plots
#generates val by freq plot, with optional log binning,
#then prints and passes the csv file to R for plotting

use vars(qw(%param %command));

sub plotHistogram {
    my ($sourceSql, $where, $fileRoot,
        $xSql, $xName, $xlab, $xmin, $xmax, $xlog, $binSize,
        $vertical) = @_;
    
    #adjust passed parameters
    $fileRoot or $fileRoot = $sourceSql;
    my $main = $fileRoot;
    $xName or $xName = $xSql;
    $xlab or $xlab = $xName;
    $xmin or $xmin = "min(data\$$xName, na.rm=TRUE)";
    $xmax or $xmax = "max(data\$$xName, na.rm=TRUE)";
    
    #retrieve the data
    my $histSql = getHistogramSql($sourceSql, $where, $xSql, $xName, $xlog, $binSize);
    $fileRoot =~ m/\/$/ or $fileRoot .= '.';
    my $file = printPlotData($histSql, "$fileRoot$xName.histogram", $xName, 'freq');

    #pass to R for plotting
    my $r = intializeRPlot($file, undef, $main,
                           $xlab, $xmin, $xmax, $xlog,
                           'Frequency', 0, "max(data\$freq, na.rm=TRUE)", undef);
    foreach my $x (@$vertical){ addRVerticalRule($r, $x) }
    addRLines($r, $xName, 'freq', 'blue');
    createRPlot($r, $file);
    
    return $file;
}

sub plotPairedHistograms {
    my ($sourceSql1, $where1, $xSql1, $xName1, 
        $sourceSql2, $where2, $xSql2, $xName2, 
        $xName, $xlab, $xmin, $xmax, $xlog, $binSize,
        $fileRoot, $vertical) = @_;
        
    #adjust passed parameters
    my $main = $fileRoot;
    $xName1 or $xName1 = $xSql1;
    $xName2 or $xName2 = $xSql2;
    $xlab or $xlab = $xName;
    $xmin or $xmin = "min(data\$$xName, na.rm=TRUE)";
    $xmax or $xmax = "max(data\$$xName, na.rm=TRUE)";
    
    #retrieve the data
    my $hist1 = getHistogramSql($sourceSql1, $where1, $xSql1, $xName1, $xlog, $binSize);
    my $hist2 = getHistogramSql($sourceSql2, $where2, $xSql2, $xName2, $xlog, $binSize);
    my ($crosstabSql) = getTwoSampleCrosstabRankSql($xName1, $xName2, $hist1, $hist2, 'val', 'freq', $xName);
    $crosstabSql = "SELECT $xName, $xName1, $xName2 FROM ($crosstabSql) ORDER BY $xName";
    $fileRoot =~ m/\/$/ or $fileRoot .= '.';
    my $file = printPlotData($crosstabSql, "$fileRoot$xName1.$xName2.histogram", $xName, $xName1, $xName2);

    #pass to R for plotting
    my $r = intializeRPlot($file, undef, $main,
                           $xlab, $xmin, $xmax, $xlog,
                           'Frequency', 0, "max(data\$$xName1, data\$$xName2, na.rm=TRUE)", undef);
    foreach my $x (@$vertical){ addRVerticalRule($r, $x) }
    addRLines($r, $xName, $xName1, 'blue');
    addRLines($r, $xName, $xName2, 'red');
    addRLegend($r, [$xName1, $xName2], ['blue', 'red']);
    createRPlot($r, $file);
        
    return $file;
}

sub getHistogramSql {
    my ($sourceSql, $where, $xSql, $xName, $xlog, $binSize) = @_;
    
    #adjust passed parameters
    if($where){ "\U$where" =~ /^WHERE/ or $where = "WHERE $where" } else { $where = "" }  
    
    #bin the values
    my $xVal = "round(($xSql)/$binSize)*$binSize";
    my $posX = "case when $xSql > 0 then $xSql else null end"; #prevent log of 0 or negative numbers
    $xlog and $xVal = "power(10, round(log(10, $posX)/$binSize)*$binSize)";
    my $valSql =   "SELECT $xVal $xName FROM ($sourceSql) $where";
    
    #calculates binned value frequencies
    my $countSql = "SELECT $xName val, count(*) N FROM ($valSql) GROUP BY $xName";
    my $sumSql =   "SELECT sum(N) FROM ($countSql)";
    my $histSql =  "SELECT val, N/($sumSql) freq FROM ($countSql) ORDER BY val";
    
    return $histSql;
}

1;

