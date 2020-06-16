#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

requireFolders_("$param{vampPath}/bin/crosstab");
requireFolders_("$param{vampPath}/bin/plot");
    
#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{xx} or $command{xx} = ['multiThread', '24:00:00', 5000, 0];
defined $param{FDR} or $param{FDR} = 0.1; #false discovery rate
defined $command{runDESeq} or $command{runDESeq} = ['singleThread', '2:00:00', 2000, 0];
defined $command{runDESeqNorm} or $command{runDESeqNorm} = ['singleThread', '2:00:00', 2000, 0];
defined $command{mergeDESeq} or $command{mergeDESeq} = ['singleThread', '2:00:00', 2000, 0];

#$testSamples and $refSamples are comma-delimited lists of samples to aggregate, such as:
#    nf0h1,nf0h2
#$testNormSamples and $refNormSamples are comma-delimited lists of samples to normalize to, such that:
#    $xxSample[i] = normalizeDECounts($xxSample[i], $xxNormSample[i])
#$xxNormSamples can be 0, in which case $testSample[i] is used directly, 
#otherwise each $xxNormSample[i] should be provided

#NOTE: normSamples was used in an attempt to allow stability assessments
#when for example two cell lines were being compared that have different 0h references
#This is NOT a sound approach as count information gets normalized away.  Scaling the 
#counts back up is invalid - genes with low expression (and therefore high variance)
#can end up all over the place after normalization.  The approach must instead somehow
#allow for propagation of count errors into the 6h/0h normalization, or something like that.
#Or, just skip the DESeq approach and just use single sample comparisons and only pay attention 
#to the extreme outliers.

sub runDESeq {
    my ($testName, $testSamples, $refName,  $refSamples) = @_;  
    runDESeqNorm($testName, $testSamples, 0, $refName,  $refSamples, 0);
}
sub runDESeqNorm {
    my ($testName, $testSamples, $testNormSamples, 
        $refName,  $refSamples,  $refNormSamples) = @_; 
    my @testSamples = split(",", $testSamples);
    my @refSamples =  split(",", $refSamples); 
    my @testNormSamples = split(",", $testNormSamples);
    my @refNormSamples =  split(",", $refNormSamples);
    my @samples = (@testSamples, @refSamples);   
    my @normSamples = (@testNormSamples, @refNormSamples);   
    my @testConds = ($testName) x scalar(@testSamples);
    my @refConds = ($refName) x scalar(@refSamples);    
    my $conds = join(",", @testConds, @refConds);   
    my $nSamples = scalar(@samples); 
    checkDESeqDirs();
    my @gmaps = getDESeqArray(\@samples, \@normSamples, $nSamples);
    my ($crosstabSql, $groupByName) = getSampleCrosstabSql("name2", "coverage", undef, undef, @gmaps);
    my $fileRoot = "$testName\_$refName";
    my $dataFile = printPlotData($crosstabSql, "plots/$fileRoot.raw_counts", "name2", @samples);
    my $rScript = "$param{vampPath}/bin/mRnaSeq/DESeq.R";
    system("Rscript $rScript $dataFile $conds $testName $refName $nSamples $param{FDR} $param{inputPath} $fileRoot");
}
sub checkDESeqDirs {
    my @dirs = qw(  plots
                    plots/MA
                    tables   
                    tables/all
                    tables/outliers);
    foreach my $dir(@dirs) { 
        my $path = "$param{inputPath}/$dir";
        -d $path or mkdir $path;
    }
}
sub getDESeqArray {
    my ($samples, $normSamples, $nSamples) = @_;
    my @gmaps;
    foreach my $i(0..($nSamples-1)){
        my $sample = $$samples[$i];
        my $normSample = $$normSamples[$i];
        my $gMap = getTableName('GMap', "$sample\_$param{geneType}");
        $normSample and $gMap = normalizeDECounts($normSample, $gMap);
        push @gmaps, [$sample, $gMap];
    }
    return @gmaps;
}
sub normalizeDECounts {
    my ($normSample, $gMap) = @_;
    #see note above, this is not a good approach...
    my $normGMap = getTableName('GMap', "$normSample\_$param{geneType}");
    my $normSql = "SELECT g.name2, nvl(g.coverage/decode(n.coverage,0,1,nvl(n.coverage,1)),0) normCov
                   FROM $gMap g, $normGMap n
                   WHERE g.name2 = n.name2";
    my $sampleSumSql = "SELECT sum(coverage) FROM $gMap";                   
    my $normSumSql = "SELECT sum(normCov) FROM ($normSql)";
    my $scalarSql = "($sampleSumSql)/($normSumSql)";   
    return "SELECT name2, round(normCov * $scalarSql) coverage FROM ($normSql)"; 
}

sub mergeDESeq {
    my ($relInputPath1, $testName1, $refName1,
        $relInputPath2, $testName2, $refName2) = @_; 
    my $key1 = "$testName1\_$refName1";      
    my $key2 = "$testName2\_$refName2";  
    my $inFile1 = "$param{inputPath}/$relInputPath1/tables/all/$key1.all.csv";
    my $inFile2 = "$param{inputPath}/$relInputPath2/tables/all/$key2.all.csv";
    my %genes;
    readDESeqOutput($inFile1, $key1, \%genes, $testName1, $refName1);
    readDESeqOutput($inFile2, $key2, \%genes, $testName2, $refName2);
    my $outputPath = "$param{inputPath}/merged";
    -d $outputPath or mkdir $outputPath;
    my $outFile = "$outputPath/$key1\_$key2.all.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile\n";
    print $outFileH ",";
    print $outFileH "$key1,,,,,,";
    print $outFileH "$key2,,,,,\n";
    print $outFileH "GENE,";
    print $outFileH "normMeanCount,$testName1,$refName1,foldChange,padj,FDR_$param{FDR},";
    print $outFileH "normMeanCount,$testName2,$refName2,foldChange,padj,FDR_$param{FDR}\n";
    foreach my $name2(sort {$a cmp $b} keys %genes){
        print $outFileH "$name2,";
        if($genes{$name2}{$key1}){
            print $outFileH join(",", @{$genes{$name2}{$key1}}).",";
        } else {
            print $outFileH ",,,,,,";
        }
        if($genes{$name2}{$key2}){
            print $outFileH join(",", @{$genes{$name2}{$key2}})."\n";
        } else {
            print $outFileH ",,,,,\n";
        } 
    }
    close $outFileH;
}
sub readDESeqOutput {
    my ($file, $key, $hash, $testName, $refName) = @_;
    open my $fileH, "<", $file or die "could not open $file\n";
    my $header = <$fileH>;
    chomp $header;
    my @header = split(",", $header);
    my %header;
    foreach my $i(0..(scalar(@header)-1)){ $header{$header[$i]} = $i } 
    while (my $line = <$fileH>){
        chomp $line;
        my @line = split(",", $line);
        $$hash{$line[$header{id}]}{$key} = 
             [ $line[$header{normMeanCount}],
               $line[$header{$refName}],
               $line[$header{$testName}],
               $line[$header{foldChange}],
               $line[$header{padj}],
               $line[$header{"FDR_$param{FDR}"}] ];
    }
    close $fileH;
}

1;
    
