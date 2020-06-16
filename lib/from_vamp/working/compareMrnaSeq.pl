#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{xx} or $command{xx} = ['multiThread', '24:00:00', 5000, 0];
defined $command{compareMrnaSeq} or $command{compareMrnaSeq} = ['singleThread', '1:00:00', 2000, 0];

sub compareMrnaSeq {
    my ($sample1, $sample2) = @_;
    my ($crosstabTable, $crosstabFile) = crosstabMrnaSeq($sample1, $sample2);  
    
    my $sample1Limit = 0.3;
    createMrnaCrosstabPlots($crosstabFile, $sample1, $sample2, $sample1Limit, $crosstabTable);
    printMrnaCrosstabLimitList($sample1, $sample2, $sample1Limit, $crosstabTable);
}

sub crosstabMrnaSeq {
    my ($sample1, $sample2) = @_;
    status("comparing hit densities per gene for samples $sample1 and $sample2...\n");    
    my $gMapTable1 = getTableName('GMap', "$sample1\_$param{geneType}");
    my $gMapTable2 = getTableName('GMap', "$sample2\_$param{geneType}");
    my $file = "$param{inputPath}/compareMrnaSeq_$sample1\_$sample2.csv";
    open my $fileH, ">", $file or die "could not open $file\n";   
    my $rankSQL = getTwoSampleCrosstab($sample1, $sample2, $gMapTable1, $gMapTable2, "name2", "NORMALIZEDDENSITY");                     
    my $refGeneGenesTable = getRefGeneTableName();                               
    my $geneSQL = "SELECT g.chromosome, g.corrstart_, g.end_, g.strand, r.*
                   FROM ($rankSQL) r, $refGeneGenesTable g
                   WHERE r.name2 = g.name2
                   ORDER BY r.rank";       
    my $crosstabTable = "gmap_$sample1\_$sample2\_$param{geneType}_ct";
    dropTable($crosstabTable);
    runSQL("CREATE TABLE $crosstabTable AS $geneSQL");                
    print $fileH "CHROMOSOME,START,END,STRAND,RANK,GENE,$sample1,$sample2,$sample1\/$sample2,$sample2\/$sample1\n";
    runSQL($geneSQL);
    my $genes = fetchAllHashRef();
    foreach my $gene(@$genes){
        print $fileH "$$gene{CHROMOSOME},$$gene{CORRSTART_},$$gene{END_},$$gene{STRAND},$$gene{RANK},$$gene{NAME2},";
        my ($s1, $s2, $s12, $s21) = 
        ($$gene{"\U$sample1"},           $$gene{"\U$sample2"},
         $$gene{"\U$sample1\_$sample2"}, $$gene{"\U$sample2\_$sample1"});
        defined $s1 or $s1 = "";
        defined $s2 or $s2 = "";
        defined $s12 or $s12 = "";
        defined $s21 or $s21 = "";
        print $fileH "$s1,$s2,$s12,$s21\n";
    }
    close $fileH;
    return ($crosstabTable, $file);
}
sub getTwoSampleCrosstab {
    my ($sample1, $sample2, $table1, $table2, $groupBy, $sum) = @_;
    my $sampleSQL1 = "SELECT '$sample1' sample, t1.* FROM $table1 t1";
    my $sampleSQL2 = "SELECT '$sample2' sample, t2.* FROM $table2 t2";
    my $crosstabSQL = "SELECT $groupBy,
                              sum(decode(sample,'$sample1',$sum,0)) $sample1,
                              sum(decode(sample,'$sample2',$sum,0)) $sample2
                       FROM  ( $sampleSQL1 UNION ALL $sampleSQL2 )
                       GROUP BY $groupBy";            
    my $ratioSQL = "SELECT $groupBy, $sample1, $sample2,
                           $sample1/nullif($sample2,0) $sample1\_$sample2,
                           $sample2/nullif($sample1,0) $sample2\_$sample1
                    FROM ($crosstabSQL)";
    my $rankSQL = "SELECT sum(1) OVER (ORDER BY $sample2\_$sample1 ROWS UNBOUNDED PRECEDING) rank,
                          $groupBy, $sample1, $sample2, $sample1\_$sample2, $sample2\_$sample1
                   FROM ($ratioSQL)
                   GROUP BY $groupBy, $sample1, $sample2, $sample1\_$sample2, $sample2\_$sample1
                   ORDER BY $sample2\_$sample1";                                     
    return $rankSQL;          
}

sub createMrnaCrosstabPlots {
    my ($crosstabFile, $sample1, $sample2, $sample1Limit, $crosstabTable) = @_;
    status("creating correlation plot...\n"); 
    system("Rscript $param{vampPath}/bin/mRnaSeq/correlation.R $crosstabFile $sample1 $sample2");
    status("creating two sample rank plot...\n"); 
    system("Rscript $param{vampPath}/bin/mRnaSeq/twoSampleRank.R $crosstabFile $sample1 $sample2"); 
    status("creating normalized density histogram plot...\n"); 
    my $gMapTable1 = getTableName('GMap', "$sample1\_$param{geneType}");
    my $gMapTable2 = getTableName('GMap', "$sample2\_$param{geneType}");
    my $hMapTable1 = getTableName('h', $gMapTable1);
    my $hMapTable2 = getTableName('h', $gMapTable2);
    my $hFile1 = createHistogramPlotFile($hMapTable1);
    my $hFile2 = createHistogramPlotFile($hMapTable2);
    system("Rscript $param{vampPath}/bin/mRnaSeq/normDensHistogram.R $hFile1 $hFile2 $sample1 $sample2 $param{inputPath}"); 
    status("creating sample ratio rank plot...\n"); 
    my $srFile = printSampleRatioRankFile($sample1, $sample2, $sample1Limit, $crosstabTable);
    system("Rscript $param{vampPath}/bin/mRnaSeq/sampleRatioRank.R $srFile $sample1 $sample2"); 
    status("creating crosstab ratio histogram plots...\n"); 
    plotCrosstabRatioHistogram($crosstabTable, "$sample1\_$sample2");
    plotCrosstabRatioHistogram($crosstabTable, "$sample2\_$sample1");
    unlink($hFile1, $hFile2, $srFile);        
}
sub createHistogramPlotFile {
    my ($hTable, $fileRoot) = @_; #$hTable can also be SQL, must have fields X and Y
    $fileRoot or $fileRoot = $hTable;
    my $countSQL = "SELECT sum(Y) FROM ($hTable)";
    runSQL("SELECT X, Y / ($countSQL) Y FROM ($hTable) ORDER BY X",\my($x,$y));
    my $hFile = "$param{inputPath}/$fileRoot.csv";
    open my $fileH, ">", $hFile;
    print $fileH "X,Y\n";
    while(fetchRow()){ print $fileH "$x,$y\n" }
    close $fileH;
    return $hFile;
}
sub printSampleRatioRankFile {
    my ($sample1, $sample2, $sample1Limit, $crosstabTable) = @_;
    runSQL("SELECT sum(1) over (order by $sample2\_$sample1 rows unbounded preceding) X, $sample2\_$sample1 Y 
            FROM $crosstabTable 
            WHERE $sample1 >= $sample1Limit
            ORDER BY $sample2\_$sample1",\my($x,$y));
    my $srFile = "$param{inputPath}/$crosstabTable.toPlot.csv";
    open my $fileH, ">", $srFile;
    print $fileH "X,Y\n";
    while(fetchRow()){ print $fileH "$x,$y\n" }
    close $fileH;
    return $srFile;
}
sub printMrnaCrosstabLimitList {
    my($sample1, $sample2, $sample1Limit, $crosstabTable) = @_;
    runSQL("SELECT sum(1) over (order by $sample2\_$sample1 rows unbounded preceding) RANK, 
                   name2,chromosome,corrstart_,end_,strand,$sample1,$sample2,$sample2\_$sample1
            FROM $crosstabTable 
            WHERE $sample1 >= $sample1Limit
            ORDER BY $sample2\_$sample1",
            \my($rank,$gene,$chrom,$start,$end,$strand,$s1,$s2,$s21));
    my $outFile = "$param{inputPath}/$crosstabTable.filtered.csv";
    open my $fileH, ">", $outFile;
    print $fileH "RANK,GENE,CHROMOSOME,START,END,STRAND,$sample1,$sample2,$sample2 \/ $sample1\n";
    while(fetchRow()){ 
        print $fileH "$rank,$gene,$chrom,$start,$end,$strand,";
        $s1 or $s1 = 0;
        $s2 or $s2 = 0; 
        $s21 or $s21 = 0; 
        print $fileH "$s1,$s2,$s21\n";
    }
    close $fileH;                 
}

sub plotCrosstabRatioHistogram {
    my ($crosstabTable, $ratioField) = @_;
    my $hSQL = "SELECT $ratioField X, count(*) Y
                FROM $crosstabTable
                GROUP BY $ratioField";
    my $fileRoot = "$crosstabTable.$ratioField";
    createHistogramPlotFile($hSQL, $fileRoot);            
    

    
           
}






1;
