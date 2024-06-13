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
my $program_version="2024/06/13";
############################## OPTIONS ##############################
use vars qw($opt_h $opt_o $opt_q $opt_t);
getopts('ho:qt:');
############################## HELP ##############################
sub help{
  print "\n";
  print "Command: $program_name [option] TABLE > SCORE\n";
  print "Arguments:\n";
  print " TABLE  VCF table calculated from vcftable.pl\n";
  print " SCORE  VCF variant score in vcftable format\n";
  print "Options:\n";
  print "    -o  Output file (default='STDOUT')\n";
  print "    -t  Quality value threshold for low quality (default<50)\n";
  print "Flag:\n";
  print "     0  wild\n";
  print "     1  hetero\n";
  print "     2  homo\n";
  print "     4  deletion\n";
  print "     8  insertion\n";
  print "    16  multi allelic (column5.count(',')>0)\n";
  print "    32  low quality (column6<QV50)\n";
  print "\n";
  print "Combinations:\n";
  print " 1 = het(1)\n";
  print " 2 = hom(2)\n";
  print " 5 = het(1)+del(4)\n";
  print " 6 = hom(2)+del(4)\n";
  print " 9 = het(1)+ins(8)\n";
  print "10 = hom(2)+ins(8)\n";
  print "17 = het(1)       +multi(16)\n";
  print "18 = hom(2)       +multi(16)\n";
  print "21 = het(1)+del(4)+multi(16)\n";
  print "22 = hom(2)+del(4)+multi(16)\n";
  print "25 = het(1)+ins(8)+multi(16)\n";
  print "26 = hom(2)+ins(8)+multi(16)\n";
  print "33 = het(1)                 +lowqual(32)\n";
  print "34 = hom(2)                 +lowqual(32)\n";
  print "37 = het(1)+del(4)          +lowqual(32)\n";
  print "38 = hom(2)+del(4)          +lowqual(32)\n";
  print "41 = het(1)+ins(8)          +lowqual(32)\n";
  print "42 = hom(2)+ins(8)          +lowqual(32)\n";
  print "49 = het(1)       +multi(16)+lowqual(32)\n";
  print "50 = hom(2)       +multi(16)+lowqual(32)\n";
  print "53 = het(1)+del(4)+multi(16)+lowqual(32)\n";
  print "54 = hom(2)+del(4)+multi(16)+lowqual(32)\n";
  print "57 = het(1)+ins(8)+multi(16)+lowqual(32)\n";
  print "58 = hom(2)+ins(8)+multi(16)+lowqual(32)\n";
  print "\n";
  print "Condition of multi allelic:\n";
  print "0/1 het\n";
  print "1/1 hom\n";
  print "0/2 het multi-allelic\n";
  print "1/2 het multi-allelic\n";
  print "2/2 hom multi-allelic\n";
  print "x/y hom multi-allelic (y>1 is multi allelic) (y==x is hom)\n";
  print "\n";
  print "Author: akira.hasegawa\@riken.jp\n";
  print "Update: $program_version\n";
  print "\n";
}
############################## MAIN ##############################
if($ARGV[0]eq"sortsubs"){sortSubs();exit();}
elsif(defined($opt_h)||scalar(@ARGV)<1){help();exit(1);}
my $vcftable=$ARGV[0];
variantScore($vcftable);
############################## getDate ##############################
sub getDate{
	my $delim=shift();
	my $time=shift();
	if(!defined($delim)){$delim="";}
	if(!defined($time)||$time eq ""){$time=localtime();}
	else{$time=localtime($time);}
	my $year=$time->year+1900;
	my $month=$time->mon+1;
	if($month<10){$month="0".$month;}
	my $day=$time->mday;
	if($day<10){$day="0".$day;}
	return $year.$delim.$month.$delim.$day;
}
############################## openFile ##############################
sub openFile{
	my $path=shift();
	if($path=~/^(.+\@.+)\:(.+)$/){
		if($path=~/\.gz(ip)?$/){return IO::File->new("ssh $1 'gzip -cd $2'|");}
		elsif($path=~/\.bz(ip)?2$/){return IO::File->new("ssh $1 'bzip2 -cd $2'|");}
		elsif($path=~/\.bam$/){return IO::File->new("ssh $1 'samtools view $2'|");}
		else{return IO::File->new("ssh $1 'cat $2'|");}
	}else{
		if($path=~/\.gz(ip)?$/){return IO::File->new("gzip -cd $path|");}
		elsif($path=~/\.bz(ip)?2$/){return IO::File->new("bzip2 -cd $path|");}
		elsif($path=~/\.bam$/){return IO::File->new("samtools view $path|");}
		else{return IO::File->new($path);}
	}
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
		if(ref( $_ ) eq "ARRAY"){
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
############################## sortSubs ##############################
sub sortSubs{
	my $path="$program_directory/$program_name";
	my $reader=openFile($path);
	my @headers=();
	my $name;
	my $blocks={};
	my $block=[];
	my $date=getDate("/");
	my @orders=();
	while(<$reader>){
		chomp;s/\r//g;
		if(/^#{30}\s*(\S+)\s*#{30}$/){
			$name=$1;
			if($name!~/^[A-Z]+$/){push(@{$block},$_);last;}
		}elsif(/^my \$program_version=\"\S+\";/){$_="my \$program_version=\"$date\";";}
		push(@headers,$_);
	}
	while(<$reader>){
		chomp;s/\r//g;
		if(/^#{30}\s*(\S+)\s*#{30}$/){
			$blocks->{$name}=$block;
			push(@orders,$name);
			$name=$1;
			$block=[];
		}
		push(@{$block},$_);
	}
	close($reader);
	if(defined($name)){$blocks->{$name}=$block;push(@orders,$name);}
	my ($writer,$file)=tempfile(DIR=>"/tmp",SUFFIX=>".pl");
	foreach my $line(@headers){print $writer "$line\n";}
	foreach my $key(sort{$a cmp $b}@orders){foreach my $line(@{$blocks->{$key}}){print $writer "$line\n";}}
	close($writer);
	return system("mv $file $path");
}
############################## variantScore ##############################
sub variantScore{
  my $vcftable=shift();
  my $reader=openFile($vcftable);
  my $label=<$reader>;chomp($label);
  print "$label\n";
  while(<$reader>){
    chomp;
    my ($chromosome,$position,$ref,$alt,@values)=split(/\t/);
    my $line="$chromosome\t$position\t$ref\t$alt";
    foreach my $value(@values){
      if(($value&1)>0){$value=1;}
      elsif(($value&2)>0){$value=2;}
      else{$value=0;}
    }
    my $size=scalar(@values);
    my $total=2*($size-1);
    for(my $i=0;$i<$size;$i++){
      my $score=0;
      for(my $j=0;$j<$size;$j++){$score+=abs($values[$i]-$values[$j]);}
      $score=$score/$total;
      $line.="\t".sprintf("%.2f",$score);
    }
    print "$line\n";
  }
}
