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
use vars qw($opt_d $opt_e $opt_f $opt_F $opt_g $opt_h $opt_o $opt_r $opt_t);
getopts('d:ef:F:ho:r:t:');
############################## HELP ##############################
sub help{
  print STDERR "\n";
  print STDERR "Command: $program_name [option] TABLE GENOME CASE\n";
  print STDERR "Arguments:\n";
  print STDERR " TABLE  Table from vcftable.pl\n";
  print STDERR "GENOME  Genome lengths\n";
  print STDERR "  CASE  Case file/directory\n";
  print STDERR "Options:\n";
  print STDERR "    -c  Use (c)ytoband file for removing telomere candidates\n";
  print STDERR "    -d  threshold (d)epth (default=10)\n";
  print STDERR "    -e  (e)xclude insertion/deletion\n";
  print STDERR "    -f  filter hom/het/N (default=hom)\n";
  print STDERR "    -F  unfilter hom/het/N (default=none)\n";
  print STDERR "    -g  genomecov mode (default='bga')\n";
  print STDERR "    -o  (o)utdir (default=STDOUT)\n";
  print STDERR "    -r  (r)egion size (default=10,000bp)\n";
  print STDERR "    -t  pick (t)op X (default=10)\n";
  print STDERR "\n";
  print STDERR "Flag:\n";
  print STDERR "     0  wild\n";
  print STDERR "     1  (het)erozygous\n";
  print STDERR "     2  (hom)ozygous\n";
  print STDERR "     4  deletion\n";
  print STDERR "     8  insertion\n";
  print STDERR "    16  multi allelic (column5.count(',')>0)\n";
  print STDERR "    32  low quality (column6<QV50)\n";
  print STDERR "    64  segmental deletion candidate (ALT/REF<1/threshold=0.5)\n";
  print STDERR "   128  segmental duplication candidate (ALT/REF>=threshold=1.5)\n";
  print STDERR "   256  low coverage candidate (DEPTH/AVG<=1/threshold=0.5)\n";
  print STDERR "   512  high coverage candidate (DEPTH/AVG>=threshold=2.0)\n";
  print STDERR "\n";
  print STDERR "NOTE: hom - look for homozygous stretch\n";
  print STDERR "      het - look for heterozygous stretch\n";
		print STDERR "      N   - look for flag with integer N (look for Flag for more detail)\n";
		print STDERR "            Flag can be summed, for example, hom(2)+insertion(8)=10\n";
  print STDERR "      https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html\n";
  print STDERR "      https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes\n";
  print STDERR "\n";
  print STDERR "Author: akira.hasegawa\@riken.jp\n";
  print STDERR "Update: $program_version\n";
  print STDERR "\n";
}
############################## MAIN ##############################
if($ARGV[0]eq"sortsubs"){sortSubs();exit();}
elsif(defined($opt_h)||scalar(@ARGV)<2){help();exit();}
my @caseFiles=@ARGV;
my $tableFile=shift(@caseFiles);
my $genomeFile=shift(@caseFiles);
if(scalar(@caseFiles)==0){@caseFiles=getSampleNames($tableFile);}
my $noIndel=$opt_e;
my $thresholdDepth=$opt_d;
my $filterMode=modeParser($opt_f,"hom");
my $unfilterMode=modeParser($opt_F);
my $regionSize=defined($opt_r)?$opt_r:10000;
my $genomecovMode=defined($opt_g)?$opt_g:"bga";
my $topX=defined($opt_t)?$opt_t:10;
my $outdir=(scalar(@caseFiles)>1)?"out":$opt_o;
foreach my $caseFile(@caseFiles){
	my $basename=getBasename($caseFile);
	my $bedfile=createBed($tableFile,$basename,$genomeFile,$noIndel,$filterMode,$unfilterMode,$regionSize,$genomecovMode);
	my $groups=groupByDepth($bedfile);
	my $writer=fileWriter($outdir,$basename,$noIndel,$filterMode,$unfilterMode,$regionSize,$topX,$thresholdDepth);
	if(defined($thresholdDepth)){pickByDepth($writer,$thresholdDepth,$groups);}
	elsif(defined($topX)){pickTopX($writer,$topX,$groups);}
}
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
  my $filterMode=shift();
  my $unfilterMode=shift();
  my $regionSize=shift();
  my $genomecovMode=shift();
  my $index=getColumnIndex($tableFile,$basename);
  my $reader=openFile($tableFile);
		my $lengths=readGenomeFile($genomeFile);
  my $half=int($regionSize/2);
  my ($fh,$tmpfile)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".bed");
		my $count=0;
  while(<$reader>){
    chomp;s/\r//g;
    my @tokens=split(/\t/);
    my $chr=$tokens[0];
    my $pos=$tokens[1];
    my $ref=$tokens[2];
    my $alt=$tokens[3];
    my $flag=$tokens[$index];
				if(defined($filterMode)&&($flag&$filterMode)==0){next;}
				if(defined($unfilterMode)&&($flag&$unfilterMode)>0){next;}
				if($noIndel&&(($flag&4)>0||($flag&8)>0)){next;}
				my $start=$pos-$half;
				my $end=$pos+$half;
				if($start<0){$start=0;}
				if($end<0){$end=0;}
				if(!exists($lengths->{$chr})){
					if(exists($lengths->{"chr$chr"})){$chr="chr$chr";}
					else{print STDERR "Unknown chromosome: $chr\n";next;}
				}
				if($start>$lengths->{$chr}){$start=$lengths->{$chr};}
				if($end>$lengths->{$chr}){$end=$lengths->{$chr};}
				print $fh "$chr\t$start\t$end\t$chr:$pos\t$flag\t.\n";
				$count++;
  }
  close($reader);
  close($fh);
		if($count==0){
			unlink($tmpfile);
			print STDERR "No stretch of regions found found\n";
			exit(1);
		}
  my ($fh2,$tmpfile2)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".bed");
  close($fh2);
  system("sort -k1,1 -k2,2n -k3,3n $tmpfile > $tmpfile2");
  #https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
  my ($fh3,$tmpfile3)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".bed");
  close($fh3);
  system("bedtools genomecov -i $tmpfile2 -g $genomeFile -$genomecovMode > $tmpfile3");
  unlink($tmpfile2);
  return $tmpfile3;
}
############################## fileWriter ##############################
sub fileWriter{
  my $outdir=shift();
  my $basename=shift();
  my $noIndel=shift();
  my $filterMode=shift();
  my $unfilterMode=shift();
  my $regionSize=shift();
  my $topX=shift();
  my $thresholdDepth=shift();
  if(!defined($outdir)){return IO::File->new(">&STDOUT");}
  mkdir($outdir);
  my $filename="$outdir/$basename.";
		if(defined($filterMode)){$filename.="filter".$filterMode;}
		if(defined($unfilterMode)){$filename.="_unfilter".$unfilterMode;}
		if(defined($regionSize)){$filename.="_region".$regionSize;}
  if(defined($thresholdDepth)){$filename.="_depth".$thresholdDepth;}
  elsif(defined($topX)){$filename.="_top${topX}";}
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
  my $tableFile=shift();
  my $basename=shift();
		my $reader=openFile($tableFile);
  my $line=<$reader>;
		while($line=~/^##/){$line=<$reader>;}
		if($line!~/^#/){return 2;}
  chomp($line);
  $line=~s/\r//g;
		close($reader);
  my @tokens=split(/\t/,$line);
  my $index=-1;
  for(my $i=4;$i<scalar(@tokens);$i++){
    if($tokens[$i]eq$basename){$index=$i;}
  }
		if($index<0){
			print STDERR	"Column not found: $basename\n";
			exit(1);
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
############################## getSampleNames ##############################
sub getSampleNames{
	my $tableFile=shift();
	my $reader=openFile($tableFile);
	my $line=<$reader>;
	while($line=~/^##/){$line=<$reader>;}
	if($line!~/^#/){return basename($tableFile);}
	close($reader);
	chomp($line);
	$line=~s/\r//g;
	my @tokens=split(/\t/,$line);
	my @samples=();
	for(my $i=2;$i<scalar(@tokens);$i++){push(@samples,$tokens[$i]);}
	return @samples;
}
############################## groupByDepth ##############################
sub groupByDepth{
  my $file=shift();
  my $reader=openFile($file);
  my @depths=();
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
############################## modeParser ##############################
sub modeParser{
	my $mode=shift();
	my $default=shift();
	if(!defined($mode)){
		if(defined($default)){$mode=$default;}
		else{return;}
	}
	if($mode eq"hom"){$mode=2;}
	elsif($mode eq"het"){$mode=1;}
	elsif($mode eq"dup"){$mode=512;}
	elsif($mode=~/^(\d+)$/){$mode=$1+0;}#convert string to integer
	elsif($mode=~/n/i){print STDERR	"Please use integer\n";}
	else{print STDERR "Unknown stretch mode: $mode\n";exit(1);}
	return $mode
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
############################## pickByDepth ##############################
sub pickByDepth{
  my $writer=shift();
  my $thresholdDepth=shift();
  my $groups=shift();
  my $hash={};
  for(my $i=0;$i<scalar(@{$groups});$i++){
    my ($chr,$start,$end,$depth)=@{$groups->[$i]};
    push(@{$hash->{$depth}},$i);
  }
  my ($fh,$tmpfile)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".bed");
  my @keys=sort{$b<=>$a}keys(%{$hash});
  foreach my $key(@keys){
    my @indeces=@{$hash->{$key}};
    foreach my $index(@indeces){
      my $group=$groups->[$index];
      my $depth=$group->[3];
      if($depth<$thresholdDepth){next;}
      print $fh join("\t",@{$group})."\n";
    }
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
############################## pickTopX ##############################
sub pickTopX{
  my $writer=shift();
  my $topX=shift();
  my $groups=shift();
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
  system("sort -k4,4nr -k1,1 -k2,2n -k3,3n $tmpfile > $tmpfile2");
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
############################## readGenomeFile ##############################
sub readGenomeFile{
	my $genomeFile=shift();
	my $reader=openFile($genomeFile);
	my $hash={};
	while(<$reader>){
		chomp;s/\r//g;
		my ($chr,$length)=split(/\t/);
		$hash->{$chr}=$length;
	}
	close($reader);
	return $hash;
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
