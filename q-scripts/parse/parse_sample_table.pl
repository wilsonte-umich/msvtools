use strict;
use warnings;

# get parameters
my ($column, $control) = @ARGV;

# open table and collect and verify the header
open my $inH, "<", $ENV{SAMPLE_TABLE} or die "could not open $ENV{SAMPLE_TABLE}\n";
my $header = <$inH>;
chomp $header;
$header =~ s/\r//g;
my $i = 0;
my %h = map {$_ => $i++} split(",", $header);
foreach my $col(qw(Experiment Cell_Line Cell_Clone Control Sample_ID)){
    defined $h{$col} or die "column $col missing from file $ENV{SAMPLE_TABLE}\nif you're sure it is there, run dos2unix\n";
}

# identify the columns that are returned only one time, i.e. unique instances
my %unique = map {$_ => 1} qw(Experiment Cell_Line Cell_Clone);

# collect all requested sample data
my %enc;
while(my $line = <$inH>){
    chomp $line;
    $line =~ s/\r//g;
    my @f = split(",", $line);
    
    # apply the filters
    if($ENV{EXPERIMENT}){
        $f[$h{Experiment}]  eq $ENV{EXPERIMENT} or next;
    }
    if($ENV{CELL_LINE}){
        $f[$h{Cell_Line}]   eq $ENV{CELL_LINE} or next;
    }
    if($ENV{CELL_CLONE}){
        $f[$h{Cell_Clone}]  eq $ENV{CELL_CLONE} or next;
    }
    if($control){ # a filter to restrict to control samples when refining the cell clone's model
        uc($f[$h{Control}]) eq 'YES' or next;
    }    
    
    # enforce uniqueness and print
    my $val = $f[$h{$column}] or next; # skip blank lines
    if(!$unique{$column} or !$enc{$val}){
        print "$val\n";
    }    
    $enc{$val}++; 
}
close $inH;
