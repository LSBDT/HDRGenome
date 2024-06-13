#!/usr/bin/perl
use strict 'vars';
use Cwd;
use File::Basename;
use File::Temp qw/tempfile tempdir/;
use FileHandle;
use Getopt::Std;
use IO::File;
use Time::localtime;
############################## HEADER ##############################
my ($program_name,$program_directory,$program_suffix)=fileparse($0);
$program_directory=Cwd::abs_path($program_directory);
my $program_path="$program_directory/$program_name";
my $program_version="2023/06/13";
############################## OPTIONS ##############################
use vars qw($opt_h $opt_n);
getopts('hn:');
############################## MAIN ##############################
#chromosome	position	NA18939_v2	NA18940_v2	NA18941_v2	NA18942_v2
my $line=<STDIN>;
chomp($line);
$line=~s/\r//g;
if($line=~/^#\s*(.+)$/){$line=$1;}
my @labels=split(/\t/,$line);
shift(@labels);#chr
shift(@labels);#pos
shift(@labels);#ref
shift(@labels);#alt
my $name=defined($opt_n)?$opt_n:"HDR";
print "track name=\"$name\" description=\"$name\" itemRgb=\"On\"\n";
while(<STDIN>){
    my ($chr,$start,$ref,$alt,@data)=split(/\t/);
    my $end=$start+1;
    my $size=scalar(@data);
    my $hetCount=0;
    my $homCount=0;
    my $nonCount=0;
    for(my $i=0;$i<$size;$i++){
        if(($data[$i]&1)>0){$hetCount++;}
        elsif(($data[$i]&2)>0){$homCount++;}
        else{$nonCount++;}
    }
    my $score=int(($homCount+0.5*$hetCount)/$size*1000);
    my $r=int(255*$homCount/$size);
    my $g=0;#int(255*$nonCount/$size);
    my $b=int(255*$hetCount/$size);
    print "$chr\t$start\t$end\t$homCount|$hetCount|$nonCount\t$score\t.\t$start\t$end\t$r,$g,$b\n";
}
############################## printTable ##############################
sub printTable{
	my @out=@_;
	my $return_type=$out[0];
	if(lc($return_type) eq "print"){$return_type=0;shift(@out);}
	elsif(lc($return_type) eq "array"){$return_type=1;shift(@out);}
	elsif(lc($return_type) eq "stderr"){$return_type=2;shift(@out);}
	else{$return_type= 2;}
	printTableSub($return_type,"",@out);
}
sub printTableSub{
	my @out=@_;
	my $return_type=shift(@out);
	my $string=shift(@out);
	my @output=();
	for(@out){
		if(ref($_)eq"ARRAY"){
			my @array=@{$_};
			my $size=scalar(@array);
			if($size==0){
				if($return_type==0){print $string."[]\n";}
				elsif($return_type==1){push(@output,$string."[]");}
				elsif($return_type==2){print STDERR $string."[]\n";}
			}else{
				for(my $i=0;$i<$size;$i++){push(@output,printTableSub($return_type,$string."[$i]=>\t",$array[$i]));}
			}
		} elsif(ref($_)eq"HASH"){
			my %hash=%{$_};
			my @keys=sort{$a cmp $b}keys(%hash);
			my $size=scalar(@keys);
			if($size==0){
				if($return_type==0){print $string."{}\n";}
				elsif($return_type==1){push( @output,$string."{}");}
				elsif($return_type==2){print STDERR $string."{}\n";}
			}else{
				foreach my $key(@keys){push(@output,printTableSub($return_type,$string."{$key}=>\t",$hash{$key}));}
			}
		}elsif($return_type==0){print "$string\"$_\"\n";}
		elsif($return_type==1){push( @output,"$string\"$_\"");}
		elsif($return_type==2){print STDERR "$string\"$_\"\n";}
	}
	return wantarray?@output:$output[0];
}