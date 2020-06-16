#!/usr/bin/perl
use strict;
use warnings;

########################################################
#FindIntersection is executes a comparsion of two sets
#of chromosomal segments to each other and finds and scores
#overlap matches between elements in each set.
########################################################

use vars(qw(%param));

our (%intInstrs, %mappedFields, %renamedFields, %revMappedFields);

require "$param{vampPath}/bin/HandleDB_DBI.pl";

sub findIntersection{
    my($instructionsRef) = @_;
    my ($inBlock, $tableType);
    foreach my $instrRef (@$instructionsRef){  
        my ($command, @parameters) = @$instrRef;   
        if ($command and $command =~ m/\s*(.+)/) { $command = $1 }
        $command or next;
        if($command eq 'findIntersection'){$inBlock = 1; next};
        $inBlock or next;
        my $parameters = join(" ", @parameters);
        if ($command and $command =~ m/\s*(.+)/) { #strip leading and trailing white space
            $command = $1;
            if ($command eq 'query' or $command eq 'target'){
                $tableType = $command;
                my ($inTable, $alias) = split(" ", $parameters);
                $intInstrs{$tableType}{table} = $inTable; 
                ($command eq 'target' and !$alias) and $alias = $inTable;
                $alias and $intInstrs{$tableType}{alias} = $alias;                   
            } elsif ($command eq 'output'){ 
                $intInstrs{$command} = "\U$parameters";
                executeFindIntersection(getFISQL());  
                %intInstrs = ();  %mappedFields = ();  %renamedFields = (); %revMappedFields = (); 
            } elsif ($command eq 'outputType'){ 
                $intInstrs{$command} = $parameters; 
            } elsif ($command eq 'fromCode'){ 
                $intInstrs{$command} = 1;           
            } elsif ($command eq 'exit'){
                goto allDone;
            } else {
                $intInstrs{$tableType}{$command} = $parameters;
            }
        }
    }
}

sub getFISQL{
    my $querySQL = getFIDataSQL('query');
    my $targetSQL = getFIDataSQL('target');  
    my ($qc, $qs, $qe, $ql);
    if($intInstrs{query}{alias}){
        ($qc, $qs, $qe, $ql) = ("CHROM_$intInstrs{query}{alias}",
                                "START_$intInstrs{query}{alias}",
                                "END_$intInstrs{query}{alias}",
                                "SIZE_$intInstrs{query}{alias}" );    
    } else {
        ($qc, $qs, $qe, $ql) = ("$intInstrs{query}{chrom}",
                                "$intInstrs{query}{start}",
                                "$intInstrs{query}{end}",
                                "SIZE_" ); 
    }
    my ($tc, $ts, $te, $tl, $isMatch, $score) = ("CHROM_$intInstrs{target}{alias}",
                                                 "START_$intInstrs{target}{alias}",
                                                 "END_$intInstrs{target}{alias}",
                                                 "SIZE_$intInstrs{target}{alias}",
                                                 "MATCH_$intInstrs{target}{alias}",
                                                 "SCORE_$intInstrs{target}{alias}" );                                                          
    my $joinSQL = "SELECT $renamedFields{query}, $renamedFields{target},
                            CASE WHEN $ts Is Null THEN 0
                                 WHEN ($ts <= $qs AND $te >= $qe) THEN 1
                                 WHEN ($ts > $qs AND $te < $qe) THEN 2
                                 WHEN ($ts <= $qe AND $te <= $qe) THEN 3.1
                                 ELSE 3.2 END OVERLAPTYPE
                     FROM ($querySQL), ($targetSQL)
                     WHERE $tc(+) = $qc AND $ts(+) <= $qe and $te(+) >= $qs ";  #could adjust this to allow adjacency rather than overlap
    my $overlapSQL =  "SELECT $renamedFields{query}, $renamedFields{target},
                            CASE OVERLAPTYPE WHEN 0 THEN 0 ELSE 1 END $isMatch,
                            CASE OVERLAPTYPE WHEN 0 THEN 0
                                             WHEN 1 THEN $ql
                                             WHEN 2 THEN $tl
                                             WHEN 3.1 THEN $te - $qs + 1
                                             ELSE $qe - $ts + 1 END OVERLAP    
                        FROM ($joinSQL) ";
    my $scoreSQL = "SELECT $renamedFields{query}, $renamedFields{target}, $isMatch, 
                            Round(Least(nvl(OVERLAP/$ql, 0), nvl(OVERLAP/$tl, 0)), 2) * 100 $score
                    FROM ($overlapSQL)"; 
    my $outSQL;
    if ($intInstrs{outputType} and $intInstrs{outputType} eq 'bestScore'){
        $outSQL = "SELECT $revMappedFields{query}, Max($isMatch) $isMatch, Max($score) $score
                    FROM ($scoreSQL)
                    GROUP BY $renamedFields{query}";
    } else {
        $outSQL = $scoreSQL;
    }
    return $outSQL;
}           

sub getFIDataSQL{
    my ($tableType, $shortTable) = @_;
    my $p = $intInstrs{$tableType};
    $$p{chrom} or $$p{chrom} = 'CHROMOSOME';
    $$p{start} or $$p{start} = 'START_';
    $$p{end} or $$p{end} = 'END_'; 
    if($$p{alias}){
        push @{$mappedFields{$tableType}}, "$$p{chrom} CHROM_$$p{alias}";
        push @{$mappedFields{$tableType}}, "$$p{start} START_$$p{alias}";
        push @{$mappedFields{$tableType}}, "$$p{end} END_$$p{alias}";
        push @{$mappedFields{$tableType}}, "$$p{end} - $$p{start} + 1 SIZE_$$p{alias}";
        push @{$renamedFields{$tableType}}, "CHROM_$$p{alias}";
        push @{$renamedFields{$tableType}}, "START_$$p{alias}";
        push @{$renamedFields{$tableType}}, "END_$$p{alias}"; 
        push @{$renamedFields{$tableType}}, "SIZE_$$p{alias}";
        push @{$revMappedFields{$tableType}}, "CHROM_$$p{alias} $$p{chrom}";
        push @{$revMappedFields{$tableType}}, "START_$$p{alias} $$p{start}";
        push @{$revMappedFields{$tableType}}, "END_$$p{alias} $$p{end}"; 
        push @{$revMappedFields{$tableType}}, "SIZE_$$p{alias} SIZE_";    
    } else {
        push @{$mappedFields{$tableType}}, "$$p{chrom}";
        push @{$mappedFields{$tableType}}, "$$p{start}";
        push @{$mappedFields{$tableType}}, "$$p{end}";
        push @{$mappedFields{$tableType}}, "$$p{end} - $$p{start} + 1 SIZE_";
        push @{$renamedFields{$tableType}}, "$$p{chrom}";
        push @{$renamedFields{$tableType}}, "$$p{start}";
        push @{$renamedFields{$tableType}}, "$$p{end}";
        push @{$renamedFields{$tableType}}, "SIZE_";
        push @{$revMappedFields{$tableType}}, "$$p{chrom}";
        push @{$revMappedFields{$tableType}}, "$$p{start}";
        push @{$revMappedFields{$tableType}}, "$$p{end}"; 
        push @{$revMappedFields{$tableType}}, "SIZE_";  
    }
    if ($$p{include}){
        $$p{include} =~ m/\s*(.*)\s*$/ and $$p{include} = $1; 
        $$p{include} =~ s/,/ /g;   
        $$p{include} =~ s/  / /g;
        $$p{include} = [split(" ", $$p{include})];
        foreach my $field(@{$$p{include}}){
            if($$p{alias}){
                push @{$mappedFields{$tableType}}, "$field $field\_$$p{alias}";
                push @{$renamedFields{$tableType}}, "$field\_$$p{alias}"; 
                push @{$revMappedFields{$tableType}}, "$field\_$$p{alias} $field";             
            } else {
                push @{$mappedFields{$tableType}}, "$field";
                push @{$renamedFields{$tableType}}, "$field"; 
                push @{$revMappedFields{$tableType}}, "$field"; 
            }
        }
    } 
    $mappedFields{$tableType} = join(", ", @{$mappedFields{$tableType}});
    $renamedFields{$tableType} = join(", ", @{$renamedFields{$tableType}});
    $revMappedFields{$tableType} = join(", ", @{$revMappedFields{$tableType}});
    my $where = '';
    $$p{where} and $where = " WHERE $$p{where} ";
    return "SELECT $mappedFields{$tableType} FROM $$p{table} $where";
}

sub executeFindIntersection{
    my ($outSQL) = @_;
    $intInstrs{fromCode} or openOracle();
    if(!$intInstrs{fromCode} and tableExists($intInstrs{output})){
        print "table $intInstrs{output} already exists.  OK to replace it? <Y to continue>";
        my $response = <stdin>;
        chomp $response;
        "\U$response" eq 'Y' or die;
    }
    dropTable($intInstrs{output});
    $intInstrs{fromCode} or print "creating table $intInstrs{output}...\n";
    runSQL("CREATE TABLE $intInstrs{output} AS $outSQL");  
    $intInstrs{fromCode} or closeDB();                          
}


1;
