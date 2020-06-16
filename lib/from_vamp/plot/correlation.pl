#!/usr/bin/perl
use strict;
use warnings;

#wholly generic functions for creating value correlation plots
#prints and passes the csv file to R for plotting

use vars(qw(%param %command));

sub plotCorrelation {
    my ($sourceSql, $where, $include, $fileRoot,
        $xSql, $xName, $xlab, $xmin, $xmax, $xlog,
        $ySql, $yName, $ylab, $ymin, $ymax, $ylog,
        $diagonal, $horizontal, $vertical, $noPlot) = @_;
    
    #adjust passed parameters
    if($where){ "\U$where" =~ /^WHERE/ or $where = "WHERE $where" } else { $where = "" }
    my $includeList = "";
    scalar(@$include) and $includeList = join(",", @$include).",";
    $fileRoot or $fileRoot = $sourceSql;
    my $main = $fileRoot;
    $xName or $xName = $xSql;
    $yName or $yName = $ySql;
    $xlab or $xlab = $xName;
    $ylab or $ylab = $yName;
    $xmin or $xmin = "min(data\$$xName, na.rm=TRUE)";
    $xmax or $xmax = "max(data\$$xName, na.rm=TRUE)";
    $ymin or $ymin = "min(data\$$yName, na.rm=TRUE)";
    $ymax or $ymax = "max(data\$$yName, na.rm=TRUE)";

    #extract the values
    my $valSql =   "SELECT $includeList $xSql $xName, $ySql $yName FROM ($sourceSql) $where";
    runSQL($valSql);
    #my $rows = fetchAllHashRef();

    #print the values to a csv file (for passing and also storage)
    $fileRoot =~ m/\/$/ or $fileRoot .= '.';
    my $file = "$param{inputPath}/$fileRoot$xName.$yName.correlation.csv";
    open my $fileH, ">", $file or die "could not open $file for writing: $!\n";
    print $fileH "$includeList$xName,$yName\n";
    while (my @vals = fetchRowArray()){
        fixOracleNulls(\@vals);
        print $fileH join(",", @vals)."\n";
    }
    close $fileH;
    
    $noPlot and return $file;
    
    #pass to R for plotting
    my $r = intializeRPlot($file, undef, $main,
                           $xlab, $xmin, $xmax, $xlog,
                           $ylab, $ymin, $ymax, $ylog);
    $diagonal and addRDiagonalRule($r);
    foreach my $y (@$horizontal){ addRHorizonatalRule($r, $y) }
    foreach my $x (@$vertical){ addRVerticalRule($r, $x) }
    addRPoints($r, $xName, $yName, 'blue', undef, 0.1);
    createRPlot($r, $file,);

    return $file;    
}


sub getCorrelationRegression_LeastSquares {
    my ($sourceSql, $where, $xSql, $ySql) = @_;
    if($where){ "\U$where" =~ /^WHERE/ or $where = "WHERE $where" } else { $where = "" }
    my $regrSql = "SELECT REGR_SLOPE($ySql, $xSql) slope, 
                          REGR_INTERCEPT($ySql, $xSql) intercept
                   FROM ($sourceSql) $where";
    runSQL($regrSql, \my($slope,$intercept));
    fetchRow();
    return ($slope, $intercept);      
}

sub getCorrelationRegression_TheilIncomplete {
    my ($sourceSql, $where, $xSql, $ySql) = @_;
    
    #created data table with filtering
    if($where){ "\U$where" =~ /^WHERE/ or $where = "WHERE $where" } else { $where = "" }
    my $filterSql = "SELECT $xSql x, $ySql y FROM ($sourceSql) $where";     
    my $filterTable = "TheilSen_filtered";
    dropTable($filterTable);
    runSQL("CREATE TABLE $filterTable AS $filterSql");
     
    #split data into lower and upper half and join points by rank within each half
    runSQL("SELECT median(x) FROM $filterTable", \my($xMedian));
    fetchRow();  
    my $lower = getThielHalfSql($filterTable, $xMedian, '<');
    my $upper = getThielHalfSql($filterTable, $xMedian, '>');
    my $joinSql = "SELECT lower.x x1, lower.y y1, upper.x x2, upper.y y2
                   FROM ($lower) lower, ($upper) upper
                   WHERE lower.pid = upper.pid"; 

    #calculate regression slope (m) as the median of the slopes of all ranked point pairs
    my $mSql = "SELECT median((y2 - y1) / (x2 - x1)) m FROM ($joinSql)";
    runSQL($mSql, \my($m));
    fetchRow();  
    
    #calculate regression intercept (b) as the median of all possible point intercepts, given m
    my $bSql = "SELECT median(y - $m*x) b FROM $filterTable";
    runSQL($bSql, \my($b));
    fetchRow();  

    return ($m, $b);
}
sub getThielHalfSql {
    my ($filterTable, $xMedian, $operator) = @_;
   return "SELECT sum(1) OVER (ORDER BY x, y ROWS UNBOUNDED PRECEDING) pid, x, y 
           FROM $filterTable 
           WHERE x $operator $xMedian";
}

sub getCorrelationRegression_TheilSen {
    my ($sourceSql, $where, $xSql, $ySql) = @_;
    
    #created data table with filtering and point IDs
    if($where){ "\U$where" =~ /^WHERE/ or $where = "WHERE $where" } else { $where = "" }
    my $filterSql = "SELECT $xSql x, $ySql y FROM ($sourceSql) $where";    
    $filterSql = "SELECT x, y, sum(1) OVER (ORDER BY x, y ROWS UNBOUNDED PRECEDING) pid
                  FROM ($filterSql)
                  ORDER BY x, y";    
    my $filterTable = "TheilSen_filtered";
    dropTable($filterTable);
    runSQL("CREATE TABLE $filterTable AS $filterSql");

    #calculate regression slope (m) as the median of the slopes of all possible point pairs
    my $joinSql = "SELECT p1.x x1, p1.y y1, p2.x x2, p2.y y2
                   FROM $filterTable p1, $filterTable p2
                   WHERE p1.pid > p2.pid
                     AND p2.x != p1.x"; 
    my $mSql = "SELECT median((y2 - y1) / (x2 - x1)) m FROM ($joinSql)";
    runSQL($mSql, \my($m));
    fetchRow();  
    
    #calculate regression intercept (b) as the median of all possible point intercepts, given m
    my $bSql = "SELECT median(y - $m*x) b FROM $filterTable";
    runSQL($bSql, \my($b));
    fetchRow();  

    return ($m, $b);
}
                                        
1;

