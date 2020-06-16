#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

sub wget{
    chdir "/home/wilsonte/arul/";
    system("wget --ftp-user=wilsonte --ftp-password=6myU7kM ftp://141.214.6.79/*");
}

1;

