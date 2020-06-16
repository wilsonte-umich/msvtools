use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE);
require "$libDir/Illumina.pl";
our ($inputDir, $arrayFmt, $project, $nameCol);
my ($inH);

# manage options
sub setOptions_list {
    setOptionValue(\$inputDir,     'input-dir');
    setOptionValue(\$arrayFmt,     'array-format', 'Illumina');    
    setOptionValue(\$project,      'project');
    setOptionValue(\$nameCol,      'name-column',  'DNA_ID');
    #-------------------------------------
    -d $inputDir or die "$error: directory not found: $inputDir\n";    
}

# main execution block
sub msvtools_list {
    
    # initialize
    setOptions_list();    
    print STDERR "$utility list: " . getTime(), "\n";
    
    if($arrayFmt eq 'Illumina'){
        getIlluminaNames();
    } else {
        die "unknown array format: $arrayFmt\n";
    }
}

1;
