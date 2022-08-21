#!/usr/bin/perl
my $outline=$ARGV[0];
my $line=<STDIN>;
chomp($line);
my @cases=split(/\t/,$line);
shift(@cases);
shift(@cases);
foreach my $case(@cases){print "$case\n";}
