use strict;
use warnings;

# subs specific to the handling of Illumina array data formats

use vars qw($utility $error
            $inputDir $arrayFmt $project $nameCol
            $sample $outputDir $maxMem
            $FALSE
            @chroms %chroms);

#-------------------------------------------------------------------------------
# variables
#-------------------------------------------------------------------------------
my %colNames = (
    CHROM   =>  'Chr',
    POS     =>  'Position',
    #QUAL    =>  'GC Score',
    #X       =>  'X',
    #Y       =>  'Y',
    #R       =>  'R',
    #Theta   =>  'Theta',
    LRR     =>  'Log R Ratio',
    BAF     =>  'B Allele Freq',
    #SNP_ID  =>  'SNP Index',
    PRB_NAME=>  'SNP Name'
);
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# extract one or all sample _names_ from the DNAReport file and print
#-------------------------------------------------------------------------------
sub getIlluminaNames{
    my ($drH, $hdr) = openIlluminaDNAReport();
    while(my $line = <$drH>){
        my @f = split(",", $line);
        my $name = $f[$$hdr{$nameCol}] or last;
        $name =~ s/\s/\_/g;
        print "$name\n"; 
    }
    close $drH; 
}
sub getIlluminaName{
    my ($returnCol) = @_;
    $returnCol or $returnCol = 'Row';    
    my ($drH, $hdr) = openIlluminaDNAReport();
    while(my $line = <$drH>){
        my @f = split(",", $line);
        my $name = $f[$$hdr{$nameCol}] or next;
        $name =~ s/\s/\_/g;
        $name eq $sample and return $f[$$hdr{$returnCol}];
    }
    close $drH; 
}
sub openIlluminaDNAReport {
    my $path = "$inputDir/Reports/$project"."*DNAReport.csv";
    my ($inH, $line) = openIlluminaReport($path, 'Row,');
    my $i = 0;
    my %hdr = map { $_ => $i++ } split(",", $line);
    defined $hdr{$nameCol} or die "$error: could not find column $nameCol in file $path\n";
    return ($inH, \%hdr);        

}
sub openIlluminaReport {
    my ($path, $leader) = @_;
    my @files = glob($path);
    if(@files == 0){
        die "$error: could not find $path\n";
    } elsif(@files > 1){
        die "$error: $path matched more than one file, expected only one\n";
    } else {
        open my $inH, "<", $files[0];
        my $line;
        do { $line = <$inH> } until ($line =~ m/^$leader/);
        return ($inH, $line);      
    } 
}
sub checkFileGlobExists {
    my ($path) = @_;
    my @files = glob($path);
    return @files > 0;
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# extract and index data for one array by sample name
#-------------------------------------------------------------------------------
sub indexIllumina{
    my ($sub) = @_;
    my $sampleRow = getIlluminaName();
    print STDERR "$utility index: loading Illumina array $sample = #$sampleRow\n";
    my $path = "$inputDir/Reports/FinalReports_StandardFormat/$project*FinalReport$sampleRow.txt";
    if(! checkFileGlobExists($path)){ # Prjt_167_Glover_14Aug2013b_FinalReport_27691.txt
        my $sampleId = getIlluminaName('DNA_Name');
        $path = "$inputDir/Reports/FinalReports_StandardFormat/$project*FinalReport_$sampleId.txt";
    }
    if(! checkFileGlobExists($path)){ 
        my $sampleId = getIlluminaName('DNA_ID');
        $path = "$inputDir/Reports/FinalReports_StandardFormat/$project*FinalReport_$sampleId.txt";
    }
    my ($inH, $line) = openIlluminaReport($path, '\[Data\]');
    $line = <$inH>;
    my $i = 0;
    my %hdr = map { $_ => $i++ } split("\t", $line);
    foreach my $key(keys %colNames){
        $hdr{$key} = $hdr{$colNames{$key}};
    }
    while (my $line = <$inH>){
        filterArrayLine($line, \%hdr, $sub);
    }
    closeHandles($inH);
}
sub parseIlluminaProbes{
    my ($sub) = @_;
    print STDERR "$utility index: parsing Illumina prototype array\n";
    my $reportDir = "$inputDir/Reports/FinalReports_StandardFormat";
    my @files = glob("$reportDir/*FinalReport*txt");
    @files or die "no file found: $reportDir\n";
    print STDERR  $files[0], "\n";
    my ($inH, $line) = openIlluminaReport($files[0], '\[Data\]');
    $line = <$inH>;
    my $i = 0;
    my %hdr = map { $_ => $i++ } split("\t", $line);
    foreach my $key(keys %colNames){
        $hdr{$key} = $hdr{$colNames{$key}};
    }
    while (my $line = <$inH>){
        filterArrayLine($line, \%hdr, $sub);
    }
    closeHandles($inH);
}
#-------------------------------------------------------------------------------





#sub importIlluminaCNVs{
#    #this sub expects that importArray and analyzeArray have already been run
#    my ($run) = @_;
#    print STDERR "importing Illumina called CNVs for run $run\n";
#    $runDir = "$param{inputPath}/Array/$run";
#    #getIlluminaNames();
#    loadIlluminaCNVData();
#}


#sub getArrayMean{
#    my ($arrayTable, $field, $thetaFilter) = @_;
#    my $pointsSQL = "SELECT $field FROM $arrayTable WHERE GCSCORE >= 0.15 AND $thetaFilter";
#    my $roundSQL = "SELECT Trunc($field/0.01)*0.01 $field FROM ($pointsSQL)";
#    my $histSQL = "SELECT $field, Count($field) N FROM ($roundSQL) GROUP BY $field";
#    runSQL("SELECT Max(N) FROM ($histSQL)", \my($maxN));    
#    fetchRow();       
#    runSQL("SELECT $field FROM ($histSQL) WHERE N = $maxN", \my($fieldAtMaxN));
#    fetchRow();    
#    my $minField = $fieldAtMaxN * 0.5;
#    my $maxField = $fieldAtMaxN * 1.5;    
#    runSQL("SELECT Avg($field) FROM ($pointsSQL) WHERE $field >= $minField AND $field <= $maxField", \my($mean));
#    fetchRow(); 
#    return $mean;
#}

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


#Sample ID
#SNP Index
#Y
#SNP Name
#Top Genomic Sequence
#Allele1 - Design
#B Allele Freq
#Allele1 - Plus
#Plus/Minus Strand
#Sample Name
#Position
#ILMN Strand
#GT Score
#Allele1 - Top
#Sample Group
#X Raw
#GC Score
#Log R Ratio
#Y Raw
#CNV Value
#SNP
#CNV Confidence
#
#Customer Strand
#Allele1 - Forward
#Allele2 - Plus
#Allele2 - Top
#Sample Index
#Allele2 - Forward
#Allele1 - AB
#X
#Allele2 - AB
#SNP Aux
#Theta
#Chr
#Cluster Sep
#R
#Allele2 - Design


#sub loadIlluminaCNVData{
#    status("  extracting CNVs:\n");
#    #find CNVs file
#    my @files = glob "$runDir/*CNV*\.csv"; 
#    my $nFiles = scalar @files;   
#    $nFiles > 1 and die "found more than one possible file with expected name $runDir/*CNV*\.csv";
#    $nFiles < 1 and die "could not find a file with expected name $runDir/*CNV*\.csv";
#    my $file = $files[0];
#    #read header line
#    open $fileH, "<", $file; 
#    my $line = <$fileH>;      
#    chomp $line;
#    $line =~ s/\r//g;
#    my @labels = split(",", $line);
#    for my $i(0..((scalar @labels) - 1)){$dataFields{$labels[$i]} = $i}      
#    #extract and store CNVs
#    my %cnvs;
#    while (my $line = <$fileH>){
#        chomp $line; 
#        $line =~ s/\r//g;
#        my @fields = split(",", $line);
#        my $sampleID = $fields[$dataFields{'SampleID'}];
#        $sampleID or next;          
#        $sampleID =~ m/^(.*) \[.*\]/;
#        $sampleID = $1;
#        my $chr = $fields[$dataFields{Chr}];      
#        $chr eq 'XY' and $chr = 'X'; #????
#        $chr eq 'MT' and $chr = 'M';  
#        my $chrom = $refSeqs{$param{refSeqBase}}{"chr$chr"};          
#        my $minPos = $fields[$dataFields{'Start'}];
#        my $maxPos = $fields[$dataFields{'End'}];   
#        my $cnvValue = $fields[$dataFields{'Value'}];   
#        push @{$cnvs{$sampleID}}, [$chrom, $minPos, $maxPos, $cnvValue];
#    }
#    close $fileH;
#    #put CNVs into NMA table if it exists
#    $fieldNames{NMA} =~ s/, DESCRIPTION//;    
#    foreach my $sampleID(keys %cnvs){
#        my $arrayTable = getTableName('Array', arrayTableName($sampleID, 'Illumina', $param{refSeqBase}));
#        tableExists($arrayTable) or next;       
#        my $nmaTable = getTableName('NMA', $sampleID);  
#        tableExists($nmaTable) or next;   
#        status("    sample $sampleID\n");
#        my $nmaFile = "$nmaTable.csv"; 
#        open my $nmaFileH, ">", $nmaFile; 
#        runSQL("DELETE FROM $nmaTable WHERE NMATYPE = 'External'");         
#        runSQL("SELECT Max(NMAID) FROM $nmaTable", \my($nmaID));   
#        fetchRow();
#        foreach my $cnvRef(@{$cnvs{$sampleID}}){
#            my ($chrom, $minPos, $maxPos, $cnvValue) = @$cnvRef;
#            $chrom or next;
#            $nmaID++;  
#            my $size = $maxPos - $minPos + 1;            
#            my $ratioSQL = "SELECT (CASE RATIO WHEN 0 THEN 0.001 ELSE RATIO END) RATIO
#                            FROM $arrayTable
#                            WHERE CHROMOSOME = $chrom AND POSITION >= $minPos and POSITION <= $maxPos";         
#            runSQL("SELECT Count(RATIO), Avg(RATIO), StdDev(RATIO) FROM ($ratioSQL)", \my($nProbesInSet, $testMean, $testStdDev));
#            fetchRow();    
#            $nProbesInSet or (print "found no probes on chromosome $chrom, $minPos to $maxPos\n" and next);               
#            my $normRatio = logBase2($testMean);       
#            runSQL("SELECT nvl(Log(2, Avg(Abs(BFREQUENCY - 0.5) + 0.5)/0.5), 0)
#                    FROM $arrayTable
#                    WHERE CHROMOSOME = $chrom AND POSITION >= $minPos and POSITION <= $maxPos
#                      AND RATIO > 0.333 AND RATIO < 3.0
#                      AND BFREQUENCY > 0.25 AND BFREQUENCY < 0.75 ", \my($normZyg)); 
#            fetchRow();                                  
#            runSQL("SELECT TESTMEAN NVALMEAN, TESTSTDEV NVALSTDEV
#                      FROM $nmaTable
#                      WHERE NMATYPE = 'RATIO_STATS_$chrom'", \my($nvalMean, $nvalStDev));      
#            fetchRow();    
#            $nvalStDev or (print "could not find RATIO_STATS in table $nmaTable for chromosome $chrom" and next); 
#            my $Z = abs($normRatio - $nvalMean) * sqrt($nProbesInSet) / $nvalStDev;  
#            print $nmaFileH join(",", $nmaID, 'External', -1, -1,
#                                      $chrom, $minPos, $maxPos, $size, $nProbesInSet, 
#                                      $testMean, $testStdDev, 0, 0,                                       
#                                      $normRatio, $normZyg, $Z,
#                                      $cnvValue)."\n";            
#        }
#        close $nmaFileH;
#        loadData($nmaFile, $nmaTable, ",", $fieldNames{NMA});
#    }
#}
#
##============================================================================
##Illumina CNV output fields
##these are manually dumped into CSV format file
##the single file contains called CNVs over the entire run, i.e. mutiple samples
##============================================================================
##SampleID = the same as DNA_Name from DNA Report file PLUS [sample#], e.g. 8488[1]
##           therefore must strip off [#] to get to true DNA name/sample ID
##----------------------------------------------------------------------------
##BookmarkType
##----------------------------------------------------------------------------
##Chr
##Start
##End
##Size
##----------------------------------------------------------------------------
##Author
##CreatedDate
##----------------------------------------------------------------------------
##Value = the called copy number, i.e. 0, 1, 3, 4, etc.
##----------------------------------------------------------------------------
##Comment
##============================================================================


1;


