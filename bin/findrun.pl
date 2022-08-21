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
my $program_version="2022/08/15";
############################## OPTIONS ##############################
use vars qw($opt_b $opt_d $opt_f $opt_h $opt_m $opt_o $opt_p $opt_r $opt_s);
getopts('bdfhm:o:p:r:s:');
############################## HELP ##############################
sub help{
  print STDERR "\n";
  print STDERR "Command: $program_name [option] TABLE CASE\n";
  print STDERR "Arguments:\n";
  print STDERR " TABLE  Table from vcftable.pl\n";
  print STDERR "  CASE  Case file/directory\n";
  print STDERR "Options:\n";
  print STDERR "    -b  Output in BED format\n";
  print STDERR "    -d  Exclude insertion/deletion\n";
  print STDERR "    -f  Full option\n";
  print STDERR "    -m  hom/het stretch (default=hom)\n";
  print STDERR "    -o  Output directory (default='out')\n";
  print STDERR "    -p  Pickup number (default>=1)\n";
  print STDERR "    -r  Region size (default=1,000,000bp)\n";
  print STDERR "    -s  Skip count (default=0)\n";
  print STDERR "\n";
  print STDERR "Author: akira.hasegawa\@riken.jp\n";
  print STDERR "Update: $program_version\n";
  print STDERR "\n";
  print STDERR "Example: $program_name -b -d -f -m het -o output -p 2 -r 2000000 -s 2 fulltable.txt case\n";
  print STDERR "  - Output in BED format\n";
  print STDERR "  - Exclude indel from calculation\n";
  print STDERR "  - Output in full format (additional length, hit and detail columns)\n";
  print STDERR "  - het stretch mode\n";
  print STDERR "  - Write files to output directory\n";
  print STDERR "  - Needs at least two matches in between\n";
  print STDERR "  - Region length have to be more than 2000000 bp\n";
  print STDERR "  - Accept region upto two mismatch\n";
  print STDERR "  - fulltable.txt from vcftable.pl output\n";
  print STDERR "  - VCF files under case directory\n";
}
############################## MAIN ##############################
if($ARGV[0]eq"sortsubs"){sortSubs();exit();}
elsif(defined($opt_h)||scalar(@ARGV)<2){help();exit();}
my $tableFile=$ARGV[0];
my @caseNames=split(/,/,$ARGV[1]);
my $noindel=$opt_d;
my $stretchMode=defined($opt_m)?lc($opt_m):"hom";
my $outDir=defined($opt_o)?$opt_o:'out';
my $pickupNumber=defined($opt_p)?$opt_p:1;
my $regionSize=defined($opt_r)?$opt_r:1000000;
my $skipCount=defined($opt_s)?$opt_s:0;
my $bedformat=$opt_b;
my $fullMode=$opt_f;
if(scalar(@caseNames)==0){
  print STDERR "\n";
  print STDERR "ERROR  No VCF files were found.\n";
  print STDERR "ERROR  Make sure you specify correct CASE file(s)/directory in the command line.\n";
  help();
  exit(1);
}
if(!-e $tableFile){
  print STDERR "\n";
  print STDERR "ERROR  Table file '$tableFile' doesn't exist.\n";
  print STDERR "ERROR  Make sure you specify correct TABLE in the command line.\n";
  help();
  exit(1);
}
mkdir($outDir);
my ($handler,$labels)=openTable($tableFile);
my $writers=matchIndex(\@caseNames,$labels,$outDir,$stretchMode,$noindel,$pickupNumber,$regionSize,$skipCount,$fullMode,$bedformat);
my $regionCount=0;
while(!eof($handler->[0])){
  checkRegion(nextTable($handler,$noindel),$writers,$stretchMode,$pickupNumber,$regionSize,$skipCount,$fullMode,$bedformat);
  $regionCount++;
  if($regionCount%10000==0){print STDERR "$regionCount...\n";}
}
print STDERR "Total number of regions: $regionCount\n";
foreach my $writer(@{$writers}){close($writer->[2]);}
print STDERR "Completed\n";
############################## absolutePath ##############################
sub absolutePath{
	my $path=shift();
	my $directory=dirname($path);
	my $filename=basename($path);
	return Cwd::abs_path($directory)."/$filename";
}
############################## checkRegion ##############################
sub checkRegion{
  my $positions=shift();
  my $writers=shift();
  my $stretchMode=shift();
  my $pickupNumber=shift();
  my $regionSize=shift();
  my $skipCount=shift();
  my $fullMode=shift();
  my $bedFormat=shift();
  my $startFlag=($stretchMode eq "hom")?2:1;
  my $endFlag=($stretchMode eq "hom")?1:2;
  my $size=scalar(@{$positions});
  foreach my $handler(@{$writers}){
    my $label=$handler->[0];
    my $index=$handler->[1];
    my $writer=$handler->[2];
    my $lastChromosome;
    my $lastEnd;
    for(my $i=0;$i<$size;$i++){
      my $chromosome=$positions->[$i]->[0];
      my $start=$positions->[$i]->[1];
      my $flag=$positions->[$i]->[$index];
      if(($flag&$startFlag)>0){#start search
        my $string=$flag;
        my $end=$start;
        my $match=0;
        my $mismatch=0;
        for(my $j=$i+1;$j<$size;$j++){
          my $length=$end-$start+1;
          my $flag2=$positions->[$j]->[$index];
          if(($flag2&$startFlag)>0){#continue search
            $end=$positions->[$j]->[1];
            $match++;
          }elsif(($flag2&$endFlag)>0){#end search
            if($mismatch<$skipCount){#search more
              if($flag2>0){$string.=" $flag2";}
              $mismatch++;
              next;
            }
            if($length<$regionSize){last;}
            if($match<=$pickupNumber){last;}
            if($chromosome eq $lastChromosome && $end == $lastEnd){last;}
            if($bedFormat){
              my $s=$start-1;
              my $name="$chromosome:$s-$end";
              if($match>1000){$match=1000;}
              print $writer "$chromosome\t$s\t$end\t$name\t$match\n";
            }else{
              print $writer "$chromosome\t$start\t$end";
              $match--;
              if(defined($fullMode)){print $writer "\t$length\t$match\t$mismatch\t$string";}
              print $writer "\n";
            }
            $lastEnd=$end;
            $lastChromosome=$chromosome;
            last;
          }
          if($flag2>0){$string.=" $flag2";}
        }
      }
    }
  }
}
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
############################## listFiles ##############################
sub listFiles{
	my @input_directories=@_;
	my $file_suffix=shift(@input_directories);
	my @input_files=();
	foreach my $input_directory (@input_directories){
		$input_directory=absolutePath($input_directory);
		if(-f $input_directory){push(@input_files,$input_directory);next;}
		elsif(-l $input_directory){push(@input_files,$input_directory);next;}
		opendir(DIR,$input_directory);
		foreach my $file(readdir(DIR)){
			if($file eq "."){next;}
			if($file eq "..") {next;}
			if($file eq ""){next;}
			$file="$input_directory/$file";
			if(-d $file){next;}
			elsif($file!~/$file_suffix$/){next;}
			push(@input_files,$file);
		}
		closedir(DIR);
	}
	return sort(@input_files);
}
############################## getBasename ##############################
sub getBasename{
  my $file=shift();
  my $basename=basename($file);
  if($basename=~/^(.+)\.g\.vcf$/i){$basename=$1}
  elsif($basename=~/^(.+)\.vcf$/i){$basename=$1}
  elsif($basename=~/^(.+)\.bcf$/i){$basename=$1}
  elsif($basename=~/^(.+)\.avinput$/i){$basename=$1}
  return $basename;
}
############################## matchIndex ##############################
sub matchIndex{
  my $files=shift();
  my $labels=shift();
  my $outDir=shift();
  my $stretchMode=shift();
  my $noindel=shift();
  my $pickupNumber=shift();
  my $regionSize=shift();
  my $skipCount=shift();
  my $fullMode=shift();
  my $bedformat=shift();
  my $names={};
  foreach my $file(@{$files}){
    my $basename=basename($file);
    if($basename=~/^(.+)\.g\.vcf$/i){$basename=$1}
    elsif($basename=~/^(.+)\.vcf$/i){$basename=$1}
    elsif($basename=~/^(.+)\.bcf$/i){$basename=$1}
    elsif($basename=~/^(.+)\.avinput$/i){$basename=$1}
    $names->{$basename}=1;
  }
  my @array=();
  for(my $i=0;$i<scalar(@{$labels});$i++){
    my $label=$labels->[$i];
    if(!exists($names->{$label})){next;}
    my $filename="$label.${stretchMode}_region${regionSize}_pick${pickupNumber}_skip${skipCount}";
    if(defined($noindel)){$filename.="_noindel";}
    if(defined($fullMode)){$filename.="_full";}
    if(defined($bedformat)){$filename.=".bed";}
    else{$filename.=".txt";}
    my $writer=IO::File->new(">$outDir/$filename");
    if(defined($bedformat)){}
    elsif(defined($fullMode)){print $writer "#Chr\tStart\tEnd\tLength\tHit\tUnhit\tGenotype\n";}
    else{print $writer "#Chr\tStart\tEnd\n";}
    push(@array,[$label,$i+2,$writer]);
  }
  return \@array;
}
############################## nextTable ##############################
sub nextTable{
  my $handler=shift();
  my $reader=$handler->[0];
  my $next=$handler->[1];
  if(!defined($next)){return ();}
  my @data=@{$next};
  my $chromosome=shift(@data);
  my $position=shift(@data);
  my @positions=([$chromosome,$position,@data]);
  my $endReached=0;
  my $chr;
  my $pos;
  while(<$reader>){
    chomp();
    ($chr,$pos,@data)=split(/\t/);
    if($chr ne $chromosome){$handler->[1]=[$chr,$pos,@data];last;}
    if(defined($noindel)){#take care of indel
      foreach my $d(@data){if(($d&4)>0||($d&8)>0){$d=0;}}
    }
    push(@positions,[$chr,$pos,@data]);
  }
  if(eof($reader)){close($reader);$handler->[1]=undef;}
  return \@positions;
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
############################## openTable ##############################
sub openTable{
  my $file=shift();
  my $reader;
  if($file=~/\.g(ip)?z$/){$reader=IO::File->new("gzip -cd $file|");}
  elsif($file=~/\.b(ip)?z2$/){$reader=IO::File->new("bzip2 -cd $file|");}
  else{$reader=IO::File->new($file);}
  my $line=<$reader>;
  chomp($line);
  my @names=split(/\t/,$line);
  shift(@names);shift(@names);
  $line=<$reader>;
  chomp($line);
  my @tokens=split(/\t/,$line);
  return ([$reader,\@tokens],\@names);
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
