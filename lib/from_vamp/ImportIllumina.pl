#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#microarray analysis, imports Illumina SNP arrays into VAMP
#########################################################################

use vars(qw(%param %types %fields %fieldNames %refSeqs));  #common parameters used by most analysis scripts
my ($runDir, @samples, $fileH, %dataFields);

sub importIllumina{
    my ($run) = @_;
    status("loading Illumina array run $run...\n");
    $runDir = "$param{inputPath}/Array/$run";
    getIlluminaNames();
    loadIlluminaSampleData();
}

sub importIlluminaCNVs{
    #this sub expects that importArray and analyzeArray have already been run
    my ($run) = @_;
    status("importing Illumina called CNVs for run $run...\n");
    $runDir = "$param{inputPath}/Array/$run";
    #getIlluminaNames();
    loadIlluminaCNVData();
}

sub getIlluminaNames{
    status("  reading sample names from DNAReport file\n");
    while (<$runDir/*DNAReport.csv>){
        open $fileH, "<", $_;
        my ($line, @fields);
        do { $line = <$fileH>; @fields = split(",", $line) } until ($fields[0] eq 'Row'); 
        my $DNA_Name_i;
        foreach my $i(0..(scalar(@fields) - 1)){if($fields[$i] eq 'DNA_Name'){$DNA_Name_i = $i; last}}
        while($line = <$fileH>){
            @fields = split(",", $line);
            $fields[$DNA_Name_i] or last;
            $fields[$DNA_Name_i] =~ s/\s/\_/g;
            push @samples, $fields[$DNA_Name_i]; 
        }
        close $fileH;
        last;
    }
}

sub loadIlluminaSampleData{
    status("  loading array data into database for samples:\n");
    foreach my $sample(@samples){
        status("    $sample\n");
        my @files = glob "$runDir/*FinalReport\_$sample\.txt";
        FRFILE: foreach my $file(@files){   
            open $fileH, "<", $file; 
            scrollToIlluminaHeader();            
            getArrayFields($fileH, \%dataFields);  
            my $arrayTable = newTable('Array', arrayTableName($sample, 'Illumina', $param{refSeqBase}));    
            my $arrayFile = writeIlluminaData($arrayTable);
            close $fileH;
            loadData($arrayFile, $arrayTable, ',', $fieldNames{Array});
            calculateR_BF($arrayTable);     
        } 
    }
}

sub scrollToIlluminaHeader{
    my $identifier;
    while(!$identifier or !($identifier =~ m/\[Data\]/)){
        my $line = <$fileH>; 
        ($identifier) = split("\t", $line);
    }
}

sub getArrayFields{
    my ($fileH, $fieldsRef) = @_;
    my $line = <$fileH>;
    chomp $line;
    $line =~ s/\r//g;
    my @labels = split("\t", $line);
    for my $i(0..((scalar @labels) - 1)){$$fieldsRef{$labels[$i]} = $i}
}

sub writeIlluminaData{
    my ($arrayTable) = @_;
    my $outFile = "$arrayTable.csv";
    open my $outH, ">", $outFile;
    while (my $line = <$fileH>){
        chomp $line; 
        $line =~ s/\r//g;
        my @fields = split("\t", $line);
        my $pos = $fields[$dataFields{Position}];                  
        $pos =~ m/^\d+$/ or next;   
        my $chr = $fields[$dataFields{Chr}];   
        $chr or next;  #for some reason, there are chromosome 0 showing up
        $chr eq 'XY' and $chr = 'X'; #????
        $chr eq 'MT' and $chr = 'M';  
        my $chrom = $refSeqs{$param{refSeqBase}}{"chr$chr"};  
        my $ratio = 2 ** $fields[$dataFields{'Log R Ratio'}];        
        my $bFreq = $fields[$dataFields{'B Allele Freq'}];                 
        my $zygosity = abs(0.5 - $bFreq) + 0.5;  
        my $cnvValue = -1;
        print $outH join(",", $chrom, $pos, 
                              $ratio, $ratio, $bFreq, $zygosity, $zygosity, 
                              $cnvValue, $fields[$dataFields{'GC Score'}], $fields[$dataFields{'X'}], $fields[$dataFields{'Y'}], $fields[$dataFields{'Theta'}], 
                              $fields[$dataFields{'R'}], $bFreq)."\n";
#                             CHROMOSOME, POSITION, 
#                             RATIO, NORMALIZEDRATIO, BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, 
#                             CNVVALUE, GCSCORE, X, Y, THETA, 
#                             R, BF
    }
    close $outH;    
    return $outFile;
}

#now that I understand better how Illumina calculates Rexp, this R correction isnt' really necessary for CNV finding
#also, side by side comparison of VAMP corrected Rtest/Rref is no better or worse than Illumina RATIOtest/RATIOref
sub calculateR_BF{
    my ($arrayTable) = @_;
    my $aaRMean = getArrayMean($arrayTable, 'R', 'THETA < 0.2');
    my $abRMean = getArrayMean($arrayTable, 'R', 'THETA > 0.2 AND THETA < 0.8');     
    my $bbRMean = getArrayMean($arrayTable, 'R', 'THETA > 0.8');
    my $corrRSQL = "SELECT CHROMOSOME, POSITION, 
                            RATIO, NORMALIZEDRATIO, BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, 
                            CNVVALUE, GCSCORE, 
                            CASE WHEN THETA < 0.2 THEN X / $aaRMean
                                 WHEN THETA < 0.8 THEN X / $abRMean
                                 ELSE X / $bbRMean
                                 END X,
                            CASE WHEN THETA < 0.2 THEN Y / $aaRMean
                                 WHEN THETA < 0.8 THEN Y / $abRMean
                                 ELSE Y / $bbRMean
                                 END Y,
                            THETA, 
                            CASE WHEN THETA < 0.2 THEN R / $aaRMean
                                 WHEN THETA < 0.8 THEN R / $abRMean
                                 ELSE R / $bbRMean
                                 END R
                     FROM $arrayTable";   
    my $calcBFSQL = "SELECT CHROMOSOME, POSITION, 
                          RATIO, NORMALIZEDRATIO, BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, 
                          CNVVALUE, GCSCORE, X, Y, THETA, R,
                          nvl((Y / nullif((Y + X), 0)),0.5) BF
                   FROM ($corrRSQL) ";    #NOTE; BF is NOT affected by the R correction above, it cancels out!
    updateTable('Array', $arrayTable, $calcBFSQL);    
}

sub getArrayMean{
    my ($arrayTable, $field, $thetaFilter) = @_;
    my $pointsSQL = "SELECT $field FROM $arrayTable WHERE GCSCORE >= 0.15 AND $thetaFilter";
    my $roundSQL = "SELECT Trunc($field/0.01)*0.01 $field FROM ($pointsSQL)";
    my $histSQL = "SELECT $field, Count($field) N FROM ($roundSQL) GROUP BY $field";
    runSQL("SELECT Max(N) FROM ($histSQL)", \my($maxN));    
    fetchRow();       
    runSQL("SELECT $field FROM ($histSQL) WHERE N = $maxN", \my($fieldAtMaxN));
    fetchRow();    
    my $minField = $fieldAtMaxN * 0.5;
    my $maxField = $fieldAtMaxN * 1.5;    
    runSQL("SELECT Avg($field) FROM ($pointsSQL) WHERE $field >= $minField AND $field <= $maxField", \my($mean));
    fetchRow(); 
    return $mean;
}

#Illumina FinalReport fields
#note, this list has been reordered!
#============================================================================
#snp and sample identifiers
#----------------------------------------------------------------------------
#SNP index
#SNP Name
#Sample Name
#Sample ID
#============================================================================
#chromsome and position
#----------------------------------------------------------------------------
#Chr
#Position
#============================================================================
#the base values of the snp
#----------------------------------------------------------------------------
#Allele1 - Top
#Allele2 - Top
#SNP
#============================================================================
#quality scores
#of the first three, GenCall trumps all, i.e. filtering against low GC also removes lowest GT and CS
#----------------------------------------------------------------------------
#GC Score = GenCall Score is a quality metric that indicates the reliability of _each genotype call_
            #range 0 to 1, low scores bad
            #Illumina recommends that you use a GenCall Score cutoff of 0.15 for Infinium products (Duo and Quad)
#GT Score = GenTrain Score, Measure of the cluster quality for the SNP
            #A number between 0 and 1 indicating how well the samples clustered for this locus
#Cluster Sep = Measure of the cluster separation for the SNP that ranges between 0 and 1
#Top Genomic Sequence = ??????????????
#============================================================================
#the raw image analysis numbers, generally not useful since must be normalized which is non-trivial
#----------------------------------------------------------------------------
#X Raw
#Y Raw
#============================================================================
#normalized polar and cartesian coordinates of allele1 x allele2
#----------------------------------------------------------------------------
#Theta
#R
#X
#Y
#============================================================================
#the calculated probe output values, derived from normalized coordinates
#----------------------------------------------------------------------------
#B Allele Freq
#Log R Ratio
#============================================================================
#CNV calls from Illumina's algorithm, CNVValue = copy number
#----------------------------------------------------------------------------
#CNV Value
#CNV Confidence
#============================================================================

sub loadIlluminaCNVData{
    status("  extracting CNVs:\n");
    #find CNVs file
    my @files = glob "$runDir/*CNV*\.csv"; 
    my $nFiles = scalar @files;   
    $nFiles > 1 and die "found more than one possible file with expected name $runDir/*CNV*\.csv";
    $nFiles < 1 and die "could not find a file with expected name $runDir/*CNV*\.csv";
    my $file = $files[0];
    #read header line
    open $fileH, "<", $file; 
    my $line = <$fileH>;      
    chomp $line;
    $line =~ s/\r//g;
    my @labels = split(",", $line);
    for my $i(0..((scalar @labels) - 1)){$dataFields{$labels[$i]} = $i}      
    #extract and store CNVs
    my %cnvs;
    while (my $line = <$fileH>){
        chomp $line; 
        $line =~ s/\r//g;
        my @fields = split(",", $line);
        my $sampleID = $fields[$dataFields{'SampleID'}];
        $sampleID or next;          
        $sampleID =~ m/^(.*) \[.*\]/;
        $sampleID = $1;
        my $chr = $fields[$dataFields{Chr}];      
        $chr eq 'XY' and $chr = 'X'; #????
        $chr eq 'MT' and $chr = 'M';  
        my $chrom = $refSeqs{$param{refSeqBase}}{"chr$chr"};          
        my $minPos = $fields[$dataFields{'Start'}];
        my $maxPos = $fields[$dataFields{'End'}];   
        my $cnvValue = $fields[$dataFields{'Value'}];   
        push @{$cnvs{$sampleID}}, [$chrom, $minPos, $maxPos, $cnvValue];
    }
    close $fileH;
    #put CNVs into NMA table if it exists
    $fieldNames{NMA} =~ s/, DESCRIPTION//;    
    foreach my $sampleID(keys %cnvs){
        my $arrayTable = getTableName('Array', arrayTableName($sampleID, 'Illumina', $param{refSeqBase}));
        tableExists($arrayTable) or next;       
        my $nmaTable = getTableName('NMA', $sampleID);  
        tableExists($nmaTable) or next;   
        status("    sample $sampleID\n");
        my $nmaFile = "$nmaTable.csv"; 
        open my $nmaFileH, ">", $nmaFile; 
        runSQL("DELETE FROM $nmaTable WHERE NMATYPE = 'External'");         
        runSQL("SELECT Max(NMAID) FROM $nmaTable", \my($nmaID));   
        fetchRow();
        foreach my $cnvRef(@{$cnvs{$sampleID}}){
            my ($chrom, $minPos, $maxPos, $cnvValue) = @$cnvRef;
            $chrom or next;
            $nmaID++;  
            my $size = $maxPos - $minPos + 1;            
            my $ratioSQL = "SELECT (CASE RATIO WHEN 0 THEN 0.001 ELSE RATIO END) RATIO
                            FROM $arrayTable
                            WHERE CHROMOSOME = $chrom AND POSITION >= $minPos and POSITION <= $maxPos";         
            runSQL("SELECT Count(RATIO), Avg(RATIO), StdDev(RATIO) FROM ($ratioSQL)", \my($nProbesInSet, $testMean, $testStdDev));
            fetchRow();    
            $nProbesInSet or (print "found no probes on chromosome $chrom, $minPos to $maxPos\n" and next);               
            my $normRatio = logBase2($testMean);       
            runSQL("SELECT nvl(Log(2, Avg(Abs(BFREQUENCY - 0.5) + 0.5)/0.5), 0)
                    FROM $arrayTable
                    WHERE CHROMOSOME = $chrom AND POSITION >= $minPos and POSITION <= $maxPos
                      AND RATIO > 0.333 AND RATIO < 3.0
                      AND BFREQUENCY > 0.25 AND BFREQUENCY < 0.75 ", \my($normZyg)); 
            fetchRow();                                  
            runSQL("SELECT TESTMEAN NVALMEAN, TESTSTDEV NVALSTDEV
                      FROM $nmaTable
                      WHERE NMATYPE = 'RATIO_STATS_$chrom'", \my($nvalMean, $nvalStDev));      
            fetchRow();    
            $nvalStDev or (print "could not find RATIO_STATS in table $nmaTable for chromosome $chrom" and next); 
            my $Z = abs($normRatio - $nvalMean) * sqrt($nProbesInSet) / $nvalStDev;  
            print $nmaFileH join(",", $nmaID, 'External', -1, -1,
                                      $chrom, $minPos, $maxPos, $size, $nProbesInSet, 
                                      $testMean, $testStdDev, 0, 0,                                       
                                      $normRatio, $normZyg, $Z,
                                      $cnvValue)."\n";            
        }
        close $nmaFileH;
        loadData($nmaFile, $nmaTable, ",", $fieldNames{NMA});
    }
}

#============================================================================
#Illumina CNV output fields
#these are manually dumped into CSV format file
#the single file contains called CNVs over the entire run, i.e. mutiple samples
#============================================================================
#SampleID = the same as DNA_Name from DNA Report file PLUS [sample#], e.g. 8488[1]
#           therefore must strip off [#] to get to true DNA name/sample ID
#----------------------------------------------------------------------------
#BookmarkType
#----------------------------------------------------------------------------
#Chr
#Start
#End
#Size
#----------------------------------------------------------------------------
#Author
#CreatedDate
#----------------------------------------------------------------------------
#Value = the called copy number, i.e. 0, 1, 3, 4, etc.
#----------------------------------------------------------------------------
#Comment
#============================================================================


1;


