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
my $program_version="2022/08/27";
############################## OPTIONS ##############################
use vars qw($opt_d $opt_h $opt_m $opt_o $opt_r $opt_t);
getopts('dhm:o:r:t:');
############################## HELP ##############################
sub help{
  print STDERR "\n";
  print STDERR "Command: $program_name [option] TABLE CASE GENOME\n";
  print STDERR "Arguments:\n";
  print STDERR " TABLE  Table from vcftable.pl\n";
  print STDERR "  CASE  Case file/directory\n";
  print STDERR "GENOME  Genome lengths\n";
  print STDERR "Options:\n";
  print STDERR "    -d  Exclude insertion/deletion\n";
  print STDERR "    -m  hom/het stretch (default=hom)\n";
  print STDERR "    -o  outdir (default=STDOUT)\n";
  print STDERR "    -r  Region size (default=10,000bp)\n";
  print STDERR "    -t  Pick top X (default=10)\n";
  print STDERR "\n";
  print STDERR "Author: akira.hasegawa\@riken.jp\n";
  print STDERR "Update: $program_version\n";
  print STDERR "\n";
}
############################## MAIN ##############################
if($ARGV[0]eq"sortsubs"){sortSubs();exit();}
elsif(defined($opt_h)||scalar(@ARGV)<2){help();exit();}
my $tableFile=shift(@ARGV);
my $caseFile=shift(@ARGV);
my $genomeFile=shift(@ARGV);
my $basename=getBasename($caseFile);
my $noIndel=$opt_d;
my $stretchMode=defined($opt_m)?lc($opt_m):"hom";
my $regionSize=defined($opt_r)?$opt_r:10000;
my $topX=defined($opt_t)?$opt_t:10;
my $bedfile=createBed($tableFile,$basename,$genomeFile,$noIndel,$stretchMode,$regionSize);
my $groups=groupByDepth($bedfile);
my $writer=fileWriter($opt_o,$basename,$noIndel,$stretchMode,$regionSize,$topX);
pickTopX($writer,$topX,$groups);
############################## absolutePath ##############################
sub absolutePath{
	my $path=shift();
	my $directory=dirname($path);
	my $filename=basename($path);
	return Cwd::abs_path($directory)."/$filename";
}
############################## createBed ##############################
sub createBed{
  my $tableFile=shift();
  my $basename=shift();
  my $genomeFile=shift();
  my $noIndel=shift();
  my $stretchMode=shift();
  my $regionSize=shift();
  my $reader=openFile($tableFile);
  my $index=getColumnIndex($reader,$basename);
  my $mode=0;
  my $half=int($regionSize/2);
  if($stretchMode eq"hom"){$mode=2;}
  elsif($stretchMode eq"het"){$mode=1;}
  my ($fh,$tmpfile)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".bed");
  while(<$reader>){
    chomp;s/\r//g;
    my @tokens=split(/\t/);
    my $chr=$tokens[0];
    my $pos=$tokens[1];
    my $flag=$tokens[$index];
    if(($flag&$mode)==0){next;}
    if($noIndel&&(($flag&4)>0||($flag&8)>0)){next;}
    my $start=$pos-$half;
    my $end=$pos+$half;
    print $fh "$chr\t$start\t$end\t$chr:$pos\t$flag\t.\n";
  }
  close($reader);
  close($fh);
  my ($fh2,$tmpfile2)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".bed");
  close($fh2);
  system("sort -k1,1 -k2,2n -k3,3n $tmpfile > $tmpfile2");
  #https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
  my ($fh3,$tmpfile3)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".bed");
  close($fh3);
  system("bedtools genomecov -i $tmpfile2 -g $genomeFile -bga > $tmpfile3");
  unlink($tmpfile2);
  return $tmpfile3;
}
############################## fileWriter ##############################
sub fileWriter{
  my $outdir=shift();
  my $basename=shift();
  my $noIndel=shift();
  my $stretchMode=shift();
  my $regionSize=shift();
  my $topX=shift();
  if(!defined($outdir)){return IO::File->new(">&STDOUT");}
  mkdir($outdir);
  my $filename="$outdir/$basename.${stretchMode}_region${regionSize}_top${topX}";
  if(defined($noIndel)){$filename.="_noindel";}
  $filename.=".txt";
  return IO::File->new(">$filename");
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
############################## getColumnIndex ##############################
sub getColumnIndex{
  my $reader=shift();
  my $basename=shift();
  my $line=<$reader>;
  chomp($line);
  $line=~s/\r//g;
  my @tokens=split(/\t/,$line);
  my $index=-1;
  for(my $i=0;$i<scalar(@tokens);$i++){
    if($tokens[$i]eq$basename){$index=$i;}
  }
  return $index;
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
############################## groupByDepth ##############################
sub groupByDepth{
  my $file=shift();
  my $reader=openFile($file);
  my @depths=();
  my $total=0;
  my $currentChr;
  my $currentStart;
  my $currentEnd;
  my $currentDepth;
  my @array=();
  while(<$reader>){
    chomp;s/\r//g;
    my ($chr,$start,$end,$depth)=split(/\t/);
    if($depth==0){
      if(defined($currentChr)&&$currentDepth>1){push(@array,[$currentChr,$currentStart,$currentEnd,$currentDepth]);}
      $currentChr=undef;
      $currentStart=undef;
      $currentEnd=undef;
      $currentDepth=undef;
      next;
    }
    if(!defined($currentChr)){
      $currentChr=$chr;
      $currentStart=$start;
      $currentEnd=$end;
      $currentDepth=$depth;
      next;
    }
    $currentEnd=$end;
    if($currentDepth<$depth){$currentDepth=$depth;}
  }
  if(defined($currentChr)&&$currentDepth>1){push(@array,[$currentChr,$currentStart,$currentEnd,$currentDepth]);}
  close($reader);
  return \@array;
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
############################## pickTopX ##############################
sub pickTopX{
  my $writer=shift();
  my $topX=shift();
  my $groups=shift();
  my $outdir=shift();
  my $hash={};
  for(my $i=0;$i<scalar(@{$groups});$i++){
    my ($chr,$start,$end,$depth)=@{$groups->[$i]};
    push(@{$hash->{$depth}},$i);
  }
  my ($fh,$tmpfile)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".bed");
  my @keys=sort{$b<=>$a}keys(%{$hash});
  my $total=0;
  foreach my $key(@keys){
    my @indeces=@{$hash->{$key}};
    foreach my $index(@indeces){
      my $group=$groups->[$index];
      print $fh join("\t",@{$group})."\n";
    }
    $total+=scalar(@indeces);
    if($total>=$topX){last;}
  }
  close($fh);
  my ($fh2,$tmpfile2)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".bed");
  close($fh2);
  system("sort -k1,1 -k2,2n -k3,3n $tmpfile > $tmpfile2");
  my $reader=openFile($tmpfile2);
  print $writer "#Chr\tStart\tEnd\tDepth\n";
  while(<$reader>){chomp;s/\r//g;print $writer "$_\n";}
  close($reader);
  close($writer);
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
