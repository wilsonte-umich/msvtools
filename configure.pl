#!/usr/bin/perl
use strict;
use warnings;

our $utilityDescription = "applications for exploring microarray data for structural variants";
our @prerequisites = qw(
R
cat
awk
gunzip
gzip
sort
bgzip
tabix);
require "./lib/configure.pl";
