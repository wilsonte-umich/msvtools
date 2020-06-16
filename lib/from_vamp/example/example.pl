#!/usr/bin/perl
use strict;
use warnings;

#=============================================================================================================
#this statement makes various vamp global variables available to this script
#   %param hash holds vamp parameters passed from the instructions file
#   e.g. $param{refSeq} might be equal to 'hg18'
use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs));
#=============================================================================================================


#=============================================================================================================
#this section defines callable parameters and commands in this script
#format:
#   defined $param{paramName} or $param{paramName} = <defaultValue>; 
#   defined $command{commandName} or $command{commandName} = [<threadType>, <wall time required>, <Mb RAM>, 0];
#<threadType> must be one of the following:
#   'singleThread' = command takes multiple parameters
#   'multiThread' = command takes only one parameter; parameter lists call many instances of command
#<wall time required> and <Mb RAM> are default values in same format as vamp 'job' parameter
defined $param{exampleParam1} or $param{exampleParam1} = 'someText'; 
defined $param{exampleParam2} or $param{exampleParam2} = 55; 
defined $command{exampleCommand1} or $command{exampleCommand1} = ['singleThread', '1:00:00', 2000, 0];
defined $command{exampleCommand2} or $command{exampleCommand2} = ['multiThread', '6:00:00', 5000, 0];
#=============================================================================================================


#=============================================================================================================
#this section contains the subs that define the callable commands in this script
#these would be called in a vamp instructions file as follows:
#   require example
#   exampleParam2 333
#   exampleCommand1 fred wilma
#   exampleCommand2 barney betty
#the result would be the following:
#   $param{exampleParam1} defaults to 'someText'
#   $param{exampleParam2} is set to 333
#   exampleCommand1 is called once with passed values 'fred' and 'wilma'
#   exampleCommand2 is called twice, once with passed value 'barney', once with 'betty'
#these commands simply pass parameters and input values on to Python and return the Python result
#the same can be done for R, or any other Linux program
#of course, you could also use Perl and VAMP functions to do work before calling Python
#--------------------------------------------------------------------------------------------------------------
sub exampleCommand1 {
    my ($passedValue1, $passedValue2) = @_; #singleThread commands can take multiple inputs
    
    my $pythonScript = "$param{vampPath}/bin/example/pythonScript.py"; #you must specify the full path, including vampPath
    my $pythonCommand = "$pythonScript $passedValue1 $passedValue2 $param{exampleParam1} $param{exampleParam2}";
    my $pythonResult = qx/$pythonCommand/; #calls Python and captures the Python output
    
    status("$pythonResult\n"); #print the result using vamp 'status' function
}
#--------------------------------------------------------------------------------------------------------------
sub exampleCommand2 { 
    my ($passedValue) = @_; #multiThread commands take only one input!

    my $pythonScript = "$param{vampPath}/bin/example/pythonScript.py"; #you must specify the full path, including vampPath
    my $pythonCommand = "$pythonScript $passedValue $param{exampleParam1} $param{exampleParam2}";
    my $pythonResult = qx/$pythonCommand/; #calls Python and captures the Python output
    
    status("$pythonResult\n"); #print the result using vamp 'status' function
}
#=============================================================================================================


#=============================================================================================================
#the following line ensure that the script returns true when required by Perl
#do not delete!
1;
#=============================================================================================================


