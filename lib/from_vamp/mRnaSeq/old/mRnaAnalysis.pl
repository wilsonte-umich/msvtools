#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
#defined $command{addMrnaMaps} or $command{addMrnaMaps} = ['multiThread', '24:00:00', 5000, 0];
defined $param{geneType} or $param{geneType} = "Unq"; 
defined $param{isMrna} or $param{isMrna} = 0; 
defined $command{parseGeneHits} or $command{parseGeneHits} = ['multiThread', '12:00:00', 2000, 0];
defined $command{stabilome} or $command{stabilome} = ['singleThread', '1:00:00', 2000, 0];

defined $command{compareMrnaSeq} or $command{compareMrnaSeq} = ['singleThread', '4:00:00', 5000, 0];

defined $command{testRunSql} or $command{testRunSql} = ['singleThread', '2:00:00', 500, 0];

sub parseGeneHits {
    my ($sample) = @_;
    getIsMrna(\my%isMrna);
    my $isMrna = $isMrna{$sample};
    #my $gMapTable = newTable('GMap', "$sample\_$param{geneType}");  
    #calculateGeneDensity($sample, $gMapTable, $isMrna);
    
    my $gMapTable = getTableName('GMap', "$sample\_$param{geneType}");  
    calculateNormalizedGeneDensity($sample, $gMapTable, $isMrna);
    
    #fix table names
}

sub calculateGeneDensity {
    my ($sample, $gMapTable, $isMrna) = @_;
    status("calculating hit density per gene...\n    Chr: ");
    foreach my $chrom (1..nChrom()){
        status("$chrom ");
        my $densitySql;
        if ($isMrna) {
            $densitySql = getMrnaDensitySQL($sample, $chrom);
        } else {
            $densitySql = getGeneDensitySQL($sample, $chrom);
        } 
        my $gMapFile = "$gMapTable.csv";
        open my $gMapFileH, ">", $gMapFile;
        runSQL($densitySql, \my($name2,$coverage,$density));   
        while(fetchRow()){ print $gMapFileH "$name2,$coverage,$density\n"} 
        close $gMapFileH;
        loadData($gMapFile, $gMapTable, ",", "NAME2, COVERAGE, DENSITY"); 
    } 
    status("\n");  
}

sub getIsMrna {
    my ($isMrna) = @_;
    $param{isMrna} or return;
    my @isMrna = split(",", $param{isMrna});
    foreach my $sample(@isMrna){ $$isMrna{$sample} = 1 }
}

sub calculateNormalizedGeneDensity {
    my ($sample, $gMapTable, $isMrna) = @_;
    status("calculating normalized hit density per gene...\n");
    my $genomeHitCount = getGenomeGeneHitCount($sample);
    my $genomeSize = getGenomeGeneSize($isMrna);
    my $genomeDensity = $genomeHitCount / $genomeSize;
    status("    $sample has $genomeDensity average gene hit density\n");
    runSQL("UPDATE $gMapTable SET NORMALIZEDDENSITY = DENSITY / $genomeDensity");
}

#sub getGenesByChrom {
#    my ($chrom, $genes) = @_; 
#    my $refGeneGenesTable = getRefGeneTableName();
#    runSQL("SELECT name2 FROM $refGeneGenesTable WHERE chromosome = $chrom",\my($name2));
#    while(fetchRow()){$$genes{$name2}++}
#}


sub stabilome {
    my ($earlySample, $lateSample) = @_;
    status("comparing hit densities per gene for early sample = $earlySample, late sample = $lateSample...\n");    
    my $earlyTable = getTableName('GMap', "$earlySample\_$param{geneType}");
    my $lateTable = getTableName('GMap', "$lateSample\_$param{geneType}");
    my $file = "$param{inputPath}/stabilome_$earlySample\_$lateSample.csv";
    open my $fileH, ">", $file or die "could not open $file\n";   
    my $rankSQL = getTwoSampleCrosstab($earlySample, $lateSample, $earlyTable, $lateTable, 1, 1,
                                       "name2", "NORMALIZEDDENSITY");
    my $refGeneGenesTable = getRefGeneTableName();                               
    my $geneSQL = "SELECT g.chromosome, g.corrstart_, g.end_, g.strand, r.*
                   FROM ($rankSQL) r, $refGeneGenesTable g
                   WHERE r.name2 = g.name2
                   ORDER BY r.rank";  
    print $fileH "CHROMOSOME,START,END,STRAND,RANK,GENE,$earlySample,$lateSample,$earlySample\/$lateSample,$lateSample\/$earlySample\n";
    runSQL($geneSQL);
    my $genes = fetchAllHashRef();
    foreach my $gene(@$genes){
        print $fileH "$$gene{CHROMOSOME},$$gene{CORRSTART_},$$gene{END_},$$gene{STRAND},$$gene{RANK},$$gene{NAME2},";
        my ($s1, $s2, $s12, $s21) = 
        ($$gene{"\U$earlySample"},             $$gene{"\U$lateSample"},
         $$gene{"\U$earlySample\_$lateSample"},$$gene{"\U$lateSample\_$earlySample"});
        defined $s1 or $s1 = "";
        defined $s2 or $s2 = "";
        defined $s12 or $s12 = "";
        defined $s21 or $s21 = "";
        print $fileH "$s1,$s2,$s12,$s21\n";
    }
    close $fileH;
    status("        see output file $file...\n");
}



sub getGeneScalars { #normalize samples onto the same scale for inter-sample comparisons
    my ($sample1, $sample2) = @_;
    $param{scaled} or return (1,1);
    status("setting samples to common number scale...\n");
    my $count1 = getGenomeGeneHitCount($sample1);
    my $count2 = getGenomeGeneHitCount($sample2);   
    my $maxCount = $count1; #scale to sample with _greatest_ hit count
    $maxCount >= $count2 or $maxCount = $count2;  
    my $scalar1 = $maxCount / $count1;
    my $scalar2 = $maxCount / $count2;
    my $scaleRef = $sample1;
    $maxCount == $count2 and $scaleRef = $sample2;
    status("    samples scaled to $scaleRef, which has the greatest hit count\n"); 
    return ($scalar1, $scalar2);
} 

sub getGenomeGeneHitCount {
    my ($sample) = @_;  
    my $gmapTable = getTableName('GMap', "$sample\_$param{geneType}");
    runSQL("SELECT sum(coverage) N from $gmapTable", \my$count);
    fetchRow();
    status("    $sample has $count gene hits\n");
    return $count;    
}

sub getGenomeGeneSize {
    my ($isMrna) = @_;
    my $refGeneGenesTable = getRefGeneTableName();
    my $sizeField = "geneSize";
    $isMrna and $sizeField = "mRnaSize";
    runSQL("SELECT sum($sizeField) N from $refGeneGenesTable", \my$size);
    fetchRow();
    status("    $param{refSeq} has $size bp in genes\n");
    return $size; 
}


# sub compareMrnaSeq {
    # my ($sample1, $sample2, $isMrna1, $isMrna2) = @_;
    # my ($scalar1, $scalar2) = getSampleScalars($sample1, $sample2);
    # status("comparing hit densities per gene for $sample1 and $sample2...\n");
    # my $file = "$param{inputPath}/geneCompare_$sample1\_$sample2.csv";
    # open my $fileH, ">", $file or die "could not open $file\n";             
    # print $fileH "GENE,CHROMOSOME,STRAND,START,END,";
    # print $fileH "$sample1\_nHits,$sample2\_nHits,$sample1\_hitDensity,$sample2\_hitDensity,$sample1\_rank,$sample2\_rank\n";
    # my $sql = getCrosstabDensitySQL($sample1, $sample2, $isMrna1, $isMrna2, $scalar1, $scalar2);
    # runSQL($sql);
    # my $genes = fetchAllHashRef();
    # foreach my $gene(@$genes){
        # print $fileH "$$gene{NAME2},$$gene{CHROMOSOME},$$gene{STRAND},$$gene{CORRSTART_},$$gene{END_},";
        # my ($s1nh, $s2nh, $s1hd, $s2hd, $s1r, $s2r) = 
           # ($$gene{"\U$sample1\_nHits"},     $$gene{"\U$sample2\_nHits"},
            # $$gene{"\U$sample1\_hitDensity"},$$gene{"\U$sample2\_hitDensity"},
            # $$gene{"\U$sample1\_rank"},      $$gene{"\U$sample2\_rank"});
        # defined $s1nh or $s1nh = 0;
        # defined $s2nh or $s2nh = 0;
        # defined $s1hd or $s1hd = 0;
        # defined $s2hd or $s2hd = 0;
        # print $fileH "$s1nh,$s2nh,$s1hd,$s2hd,$s1r,$s2r\n";
    # }
    # close $fileH;
    # status("        see output file $file...\n");
# }


#sub testRunSql {
#    
#    my $sql = getInGeneHitsSQL("bru_30m");
#    
#    runSQL($sql);
#    my $rows = fetchAllHashRef();
#    my $counter = 0;
#    foreach my $row(@$rows){
#        foreach my $key (sort {$a cmp $b} keys %$row){
#            print "$key\t$$row{$key}\n";
#        }
#        print "\n";
#        $counter++;
#        $counter > 1 and last;      
#    }
#}




1;
