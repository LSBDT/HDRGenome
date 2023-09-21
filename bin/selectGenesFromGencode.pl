#!/usr/bin/perl
while(<STDIN>){
    if(/^#/){next;}
    chomp;
    s/\r//g;
    my $Line=$_;
    my ($seqname,$source,$feature,$start,$end,$score,$strand,$frame,$attribute)=split(/\t/);
    if($feature ne "gene"){next;}
    my @attrs=split(/\;\s*/,$attribute);
    my $geneId;
    my $geneName;
    my $geneType;
    foreach my $attr(@attrs){
        #gene_id "ENSG00000223972.5"
        if($attr=~/gene_id\s\"(.+)\"/){$geneId=$1;}
        if($attr=~/gene_name\s\"(.+)\"/){$geneName=$1;}
        if($attr=~/gene_type\s\"(.+)\"/){$geneType=$1;}
    }
    if(!defined($geneId)){print STDERR "geneId not found for $line\n";next;}
    print "$seqname\t$start\t$end\t$geneName\t$score\t$strand\n";
}
