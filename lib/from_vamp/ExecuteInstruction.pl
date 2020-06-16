#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#executeInstruction.pl is called by vamp.pl via the qsub job file.
#executeInstruction.pl is responsible for actually calling the  
#work scripts needed by the queued job.
#########################################################################

my ($paramPass, $command, @parameters) = @ARGV;
our $paramCat = $paramPass;

#initialize status reporting
my $pJoin = join("_", @parameters);
my $statusFile = "$command\_$pJoin.status.txt";
my ($nTries, $statusFileH) = (0);
while($nTries < 2 and !(open $statusFileH, ">>", $statusFile)){
    my ($sec, $min, $hr, $day, $month, $year) = localtime(time);
    $statusFile = "$command\_$day\_$month\_$year\_$hr\_$min.status.txt";
    $nTries++;
}
close $statusFileH;
unlink $statusFile;
reportTime("\nVAMP started: ");

#recover vamp's operating parameters and the command to execute
our %param = ();
our %command = (); #not used, prevents user required scripts from balking
my @paramJoins = split("##", $paramPass);
foreach my $paramJoin (@paramJoins) {
    my ($param, $value) = split("#", $paramJoin);
    $param{$param} = $value
}
$param{passPath} or $param{passPath} = "";
$param{bowtiePath} or $param{bowtiePath} = "";
$param{blastPath} or $param{blastPath} = "";
$param{samPath} or $param{samPath} = "";
$param{tophatPath} or $param{tophatPath} = "";

#check Sun Grid job dependencies (since SGE job hold does not check job success)
if($param{qType} eq 'sunGrid' and !($param{noQueue})){
    sleep 60; #make sure q actions have had a chance to finish
    my $dependencies = pop @parameters;
    if ($dependencies) {
        foreach my $jobID(split(",", $dependencies)){
            my $qacct = qx/qacct -j $jobID/;
            $qacct or die "failed dependency on job ID $jobID: qacct has no knowledge of job\n";
            my @qacctValues = split("\n", $qacct);
            my %qacctValues;
            foreach my $qacctValue(@qacctValues){
                $qacctValue =~ s/ {1,}/\t/; #NOT greedy!
                my ($key, $value) = split("\t", $qacctValue);
                $key or next;
                $value or next;
                $qacctValues{$key} = $value;
            }
            $qacctValues{exit_status} != 0 and die "failed dependency on job ID $jobID: exit_status = $qacctValues{exit_status}\n";
        }
    }
}

#initialize schema and nested requires
our (%types, %fields, %fieldNames, %partitions, %refSeqs, $uid, %archiveTables);
require "$param{vampPath}/bin/Schema.pl"; 

#provide feedback regarding command and parameters
status("\nexecuting VAMP command:\n  $command @parameters\n\n");
status("with parameters:\n");
foreach my $param (sort {$a cmp $b} keys %param){$param{$param} and status("  $param $param{$param}\n")}
status("\n");

#execute command
openOracle();
my $sub = \&{$command};
&$sub(@parameters);
closeDB();

#close status reporting
reportTime("\nVAMP done: ");

sub reportTime{
    my ($message) = @_;
    my ($sec, $min, $hr, $day, $month, $year) = localtime(time);
    $year = $year + 1900;
    $month++;
    status("$message $month/$day/$year at $hr:$min:$sec\n");
}

sub status{ #newline must come from status string
    my ($status) = @_;
    print $status;
    open my $statusFileH, ">>", $statusFile or die "could not open $statusFile";
    print $statusFileH $status;
    close $statusFileH;
}

sub printParameters{
    my ($paramFile, $command, @parameters) = @_;
    open my $paramFileH, ">", $paramFile;  
    print $paramFileH "Command Parameters:\n";
    print $paramFileH join(" ", $command, @parameters)."\n\n";
    print $paramFileH "General Parameters:\n";
    foreach my $param (sort {$a cmp $b} keys %param){$param{$param} and print $paramFileH "$param $param{$param}\n"}  
    close $paramFileH;    
}

1;

