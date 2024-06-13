#!/usr/bin/perl
my $outline=$ARGV[0];
my $line=<STDIN>;
chomp($line);
$line=~s/\r//g;
my @cases=split(/\t/,$line);
shift(@cases);#chr
shift(@cases);#pos
shift(@cases);#ref
shift(@cases);#alt
foreach my $case(@cases){print "$case\n";}
