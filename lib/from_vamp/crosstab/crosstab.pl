#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command));

#crosstab.pl provides wholly generic internal functions
#for generating Oracle crosstab SQL

sub getTwoSampleCrosstabRankSql {
    #generates crosstab of just two inputs, allowing for calculation of
    #both sample ratios and sample fractions
    #assumes input data are normalized to same scale
    my ($sample1, $sample2, $inputSql1, $inputSql2, $groupBySql, $sumSql, $groupByName, $sumName) = @_;
    my $crosstabSql;
    ($crosstabSql, $groupByName) = getSampleCrosstabSql($groupBySql, $sumSql, $groupByName, $sumName,
                                                        [$sample1, $inputSql1], [$sample2, $inputSql2]);         
    my ($ratioSql, $ratio1Name, $ratio2Name, $fractionName) = getTwoSampleRatioSql($crosstabSql, $groupByName, $sample1, $sample2);                                                            
    my $rankSql = getCrosstabRankSql($ratioSql, $groupByName, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName);
    return ($rankSql, $ratio1Name, $ratio2Name, $fractionName);
}
sub getNormalizedCrosstabSql {
    my ($test, $ref) = @_;
    my ($testSql, $testRatio1Name, $testRatio2Name, $testFractionName) = getTwoSampleCrosstabRankSql(@$test);
    my ($refSql,  $refRatio1Name,  $refRatio2Name,  $refFractionName) =  getTwoSampleCrosstabRankSql(@$ref);
    my ($sample1, $sample2, $inputSql1, $inputSql2, $groupBySql, $sumSql, $groupByName, $sumName) = @$test;    
    $groupByName or $groupByName = $groupBySql;
    my $normSql = "SELECT t.$groupByName, t.$sample1, t.$sample2, 
                          t.$testRatio1Name / nullif(r.$refRatio1Name,0) $testRatio1Name,
                          t.$testRatio2Name / nullif(r.$refRatio2Name,0) $testRatio2Name,
                          t.$testFractionName
                   FROM ($testSql) t, ($refSql) r
                   WHERE t.$groupByName = r.$groupByName";
    my $rankSql = getCrosstabRankSql($normSql, $groupByName, $sample1, $sample2, $testRatio1Name, $testRatio2Name, $testFractionName);                    
    return ($rankSql, $testRatio1Name, $testRatio2Name, $testFractionName);
}

sub getSampleCrosstabSql {
    my ($groupBySql, $sumSql, $groupByName, $sumName, @samples) = @_;
    $groupByName or $groupByName = $groupBySql;
    $sumName or $sumName = $sumSql;   
    my (@unionSql, @sampleSql);
    foreach my $sample(@samples){
        my ($sampleName, $inputSql) = @$sample;
        push @unionSql, getCrosstabInputSql($sampleName, $inputSql, $groupBySql, $sumSql, $groupByName, $sumName);
        push @sampleSql, "sum(decode(sample,'$sampleName',$sumName,0)) $sampleName"
    }
    my $unionSql = join(" UNION ALL ", @unionSql);
    my $sampleSql = join(", ", @sampleSql);
    my $crosstabSql = "SELECT $groupByName, $sampleSql
                       FROM  ( $unionSql )
                       GROUP BY $groupByName";
    return ($crosstabSql, $groupByName);
}
sub getCrosstabInputSql {
    my ($sampleName, $inputSql, $groupBySql, $sumSql, $groupByName, $sumName) = @_;
    return "SELECT '$sampleName' sample,
                    $groupBySql $groupByName,
                    $sumSql $sumName
            FROM ($inputSql)";
}

sub getSampleMultiCrosstabSql {
    my ($groupBys, $sums, $samples, $where) = @_;
    my ($groupBySqls, $groupByNames) = parseCrosstabMultiVals($groupBys);
    my ($sumSqls, $sumNames, @sumNames) = parseCrosstabMultiVals($sums);
    my (@unionSql, @sampleSql, @sampleSumNames);
    foreach my $sample(@$samples){
        my ($sampleName, $inputSql) = @$sample;
        push @unionSql, "SELECT '$sampleName' sample, $groupBySqls, $sumSqls FROM ($inputSql)";
        foreach my $sumName(@sumNames){
            my $sampleSumName = "$sampleName\_$sumName";
            push @sampleSql, "sum(decode(sample,'$sampleName',$sumName,0)) $sampleSumName";
            push @sampleSumNames, $sampleSumName;
        }
    }
    my $unionSql = join(" UNION ALL ", @unionSql);
    my $sampleSql = join(",", @sampleSql);
    my $sampleSumNames = join(",", @sampleSumNames);
    my $crosstabSql = "SELECT $groupByNames, $sampleSql
                       FROM  ( $unionSql )
                       GROUP BY $groupByNames";    
    if ($where){
        "\U$where" =~ m/^WHERE/ or $where = "WHERE $where";
        $crosstabSql = "SELECT * FROM ($crosstabSql) $where";
    }
    return ($crosstabSql, $groupByNames, $sampleSumNames);
}
sub parseCrosstabMultiVals {
    my ($vals) = @_;
    my (@sqls, @names);
    foreach my $val(@$vals){
        my ($sql, $name) = @$val;
        $name or $name = $sql;
        push @sqls, "$sql $name";
        push @names, $name;
    }
    my $sqls = join(",", @sqls);
    my $names = join(",", @names);
    return ($sqls, $names, @names);
}  

sub getTwoSampleRatioSql {
    my ($crosstabSql, $groupByName, $sample1, $sample2) = @_;
    my $ratio1Name = "$sample1\_$sample2";
    my $ratio2Name = "$sample2\_$sample1";
    my $fractionName = "fraction_$sample1"; 
    my $zeroSql = "SELECT $groupByName, nvl($sample1,0) $sample1, nvl($sample2,0) $sample2
                   FROM ($crosstabSql)";
    my $ratioSql = "SELECT $groupByName, $sample1, $sample2,
                           $sample1/nullif($sample2,0) $ratio1Name,
                           $sample2/nullif($sample1,0) $ratio2Name,
                           $sample1/nullif($sample1 + $sample2,0) $fractionName
                    FROM ($zeroSql)";
    return ($ratioSql, $ratio1Name, $ratio2Name, $fractionName);
}
sub getCrosstabRankSql {
    my ($sql, $groupByName, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName) = @_;
    my $rankSql = "SELECT sum(1) OVER (ORDER BY $ratio1Name ROWS UNBOUNDED PRECEDING) rank,
                          $groupByName, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName
                   FROM ($sql)
                   GROUP BY $groupByName, $sample1, $sample2, $ratio1Name, $ratio2Name, $fractionName
                   ORDER BY $ratio1Name";
    return $rankSql;
}

sub updateCrosstabTable {
    my ($crosstabTable, $selectSql) = @_;
    my $tableTMP = "$crosstabTable\_tmp";
    $tableTMP =~ s/_//g;
    dropTable($tableTMP);
    runSQL("CREATE TABLE $tableTMP AS $selectSql");    
    dropTable($crosstabTable);
    runSQL("RENAME $tableTMP TO $crosstabTable");
}

1;

