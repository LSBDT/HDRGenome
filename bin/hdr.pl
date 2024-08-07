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
use vars qw($opt_d $opt_e $opt_h $opt_i $opt_l $opt_o $opt_s $opt_t);
getopts('de:hi:lo:s:t:');
############################## HELP ##############################
sub help{
  print STDERR "Command: $program_name [OPTIONS] TABLE CASE POS [POS2 ..]\n";
  print STDERR "Arguments:\n";
  print STDERR " TABLE  Table file from vcftable.pl\n";
  print STDERR "  CASE  Case file/directory\n";
  print STDERR "   POS  Position/region file\n";
  print STDERR "Options:\n";
  print STDERR "    -d  excluded/skip indel hom/het\n";
  print STDERR "    -e  end distance (default=50)\n";
  print STDERR "    -i  interval (default=10)\n";
  print STDERR "    -l  exclude low quality hom/het (default='include')\n";
  print STDERR "    -o  outdir (default=\"out\")\n";
  print STDERR "    -s  start distance (default=10)\n";
  print STDERR "    -t  Target mode (default='AR')\n";
  print STDERR "Mode:\n";
  print STDERR "    AR  hom vs non-hom (1/1 vs 1/0,0/0) both HOM=0\n";
  print STDERR "    DD  het vs non-het (1/0 vs 1/1,0/0) both HET=0\n";
  print STDERR "    AD  all vs all     (1/1 vs 1/0,0/0) exact match=0\n";
  print STDERR "Note:\n";
  print STDERR "    If you are using BCF files, please install bcftools\n";
  print STDERR "    http://samtools.github.io/bcftools/bcftools.html\n";
  print STDERR "\n";
  print STDERR "    Position format is 'CHROM POS' (tab delim)\n";
  print STDERR "    Region format is 'CHROM START END' (tab delim)\n";
  print STDERR "    Chromosome can be either number or 'chr'+number.\n";
  print STDERR "    For example: '1' or 'chr1'\n";
  print STDERR "    Multiple case names can be specified by ','\n";
  print STDERR "    Control names are retrieved from vcf table\n";
  print STDERR "\n";
  print STDERR "Author: akira.hasegawa\@riken.jp\n";
  print STDERR "Update: $program_version\n";
  exit(1);
}
############################## MAIN ##############################
if($ARGV[0]eq"sortsubs"){sortSubs();exit();}
elsif(defined($opt_h)||scalar(@ARGV)<3){help();exit();}
my $tableFile=$ARGV[0];
my $caseNames=[split(/,/,$ARGV[1])];
my @posFiles=listFiles("\\.te?xt\$",$ARGV[2]);
my $ctrlNames=retrieveCtrlNamesFromTable($tableFile);
$ctrlNames=removeCaseFromCtrl($caseNames,$ctrlNames);
print STDERR "==================== setting ====================\n";
if(scalar(@{$caseNames})==0){
  print STDERR "ERROR  Number of case is 0.";
  print STDERR "ERROR  Make sure you specify correct case names in the command line.\n";
  exit(1);
}else{
  print STDERR "Cases: ".scalar(@{$caseNames})."\n";
  foreach my $case(@{$caseNames}){print STDERR "  $case\n";}
}
if(scalar(@{$ctrlNames})==0){
  print STDERR "ERROR  Number of control is 0.";
  print STDERR "ERROR  Make sure you specify correct control names in the VCF table.\n";
  exit(1);
}else{
  print STDERR "Controls: ".scalar(@{$ctrlNames})."\n";
  foreach my $ctrl(@{$ctrlNames}){print STDERR "  $ctrl\n";}
}
if(scalar(@posFiles)==0){
  print STDERR "ERROR  Number of position file is 0.";
  print STDERR "ERROR  Make sure you specify correct position files in the command line.\n";
  exit(1);
}else{
  print STDERR "Position files: ".scalar(@posFiles)."\n";
  foreach my $posFile(@posFiles){print STDERR "  $posFile\n";}
}
my $outdir=(defined($opt_o))?$opt_o:"out";
my $noindel=$opt_d;
my $nolowqc=$opt_l;
my $targetMode=(defined($opt_t))?uc($opt_t):"AR";
my $windowInterval=(defined($opt_i))?$opt_i:10;
my $windowStart=(defined($opt_s))?$opt_s:10;
my $windowEnd=(defined($opt_e))?$opt_e:50;
mkdir($outdir);
print STDERR "Output directory: $outdir\n";
if(defined($noindel)){print STDERR "Include indel:false\n";}else{print STDERR "Include indel:true\n";}
if(defined($nolowqc)){print STDERR "Include low quality: false\n";}else{print STDERR "Include low quality: true\n";}
print STDERR "Target mode: $targetMode\n";
print STDERR "Window interval: $windowInterval\n";
print STDERR "Window start: $windowStart\n";
print STDERR "Window end: $windowEnd\n";
my $matchNames=matchNames($caseNames,$ctrlNames);
for(my $i=0;$i<scalar(@{$caseNames});$i++){
  my $caseName=$caseNames->[$i];
  my $matchName=$matchNames->[$i];
  mkdir("$outdir/$caseName");
  foreach my $posFile(@posFiles){
    my ($label,$outFile)=createLabel($caseName,$ctrlNames,$posFile,$windowStart,$windowEnd,$windowInterval);
    if(defined($noindel)){$outFile.="_noindel";}
    if(defined($nolowqc)){$outFile.="_nolowqc";}
    $outFile="$outdir/$caseName/$outFile.txt";
    calculateHDR($matchName,$targetMode,$noindel,$nolowqc,$tableFile,$posFile,$label,$outFile,$windowStart,$windowEnd,$windowInterval);
  }
}
############################## absolutePath ##############################
sub absolutePath{
	my $path=shift();
	my $directory=dirname($path);
	my $filename=basename($path);
	return Cwd::abs_path($directory)."/$filename";
}
############################## calculateHDR ##############################
sub calculateHDR{
  my $matchNames=shift();
  my $targetMode=shift();
  my $noindel=shift();
  my $nolowqc=shift();
  my $tableFile=shift();
  my $posFile=shift();
  my $header1=shift();
  my $outputFile=shift();
  my $windowStart=shift();
  my $windowEnd=shift();
  my $windowInterval=shift();
	print STDERR "==================== $outputFile ====================\n";
  my ($fh,$tmpfile)=tempfile(DIR=>$outdir,TEMPLATE=>'hdrXXXXXX',SUFFIX=>'.txt');
  print STDERR "Calculating position: $tmpfile\n";
  my $handler=openTable($tableFile,$matchNames);
  my ($reader,$type)=openPosition($posFile);
  my $basename=basename($tableFile);
  my $matchCount=scalar(@{$matchNames});
  if($type eq "position"){#AR/AD
    my @positions=();
    my $posCount=0;
    while(!eof($reader)){
      my ($chr,$pos)=nextPosition($reader);
      my $size=int(($windowEnd-$windowStart)/$windowInterval)+1;
      my @windows=();
      my $window;
      for(my $i=0;$i<$size;$i++){
        $window=1000*int($windowStart/2+$windowInterval/2*$i);
        push(@windows,[$pos-$window,$pos+$window]);
      }
      @positions=nextTable($handler,$chr,$pos-$window,$pos+$window,@positions);
      my @hits=();
      for(my $i=0;$i<$matchCount;$i++){
        my ($label,$index1,$index2)=@{$matchNames->[$i]};
        foreach my $position(@positions){
          my ($chrom,$posit,@data)=@{$position};
          my $flag1=$data[$index1];
          my $flag2=$data[$index2];
          if($flag1==0&&$flag2==0){next;}
          if($noindel){#in case of noindel option
            if(($flag1&4)>0||($flag2&4)>0){next;}
            if(($flag1&8)>0||($flag2&8)>0){next;}
          }
          if($nolowqc){#don't include $nolowqc
            if(($flag1&32)>0||($flag1&32)>0){next;}
          }
          push(@{$hits[$i]},[$chrom,$posit,$flag1,$flag2]);
        }
      }
      my $results=[];
      for(my $i=0;$i<scalar(@windows);$i++){
        $results->[$i]=[];
        my $window=$windows[$i];
        for(my $j=0;$j<$matchCount;$j++){
          $results->[$i]->[$j]->[0]=0;
          $results->[$i]->[$j]->[1]=0;
          foreach my $hit(@{$hits[$j]}){
            my ($chrom,$posit,$flag1,$flag2)=@{$hit};
            my $f1=($flag1&3);
            my $f2=($flag2&3);
            if($window->[0]<=$posit&&$posit<=$window->[1]){
              if($targetMode eq "AR"){
                if($f1==2||$f2==2){
                  $results->[$i]->[$j]->[0]++;
                  if($f1!=$f2){$results->[$i]->[$j]->[1]++;}
                }
              }elsif($targetMode eq "DD"){
                if($f1==1||$f2==1){
                  $results->[$i]->[$j]->[0]++;
                  if($f1!=$f2){$results->[$i]->[$j]->[1]++;}
                }
              }elsif($targetMode eq "AD"){
                if($f1>0||$f2>0){
                  $results->[$i]->[$j]->[0]++;
                  if($f1!=$f2){$results->[$i]->[$j]->[1]++;}
                }
              }
            }
          }
        }
      }
      for(my $i=0;$i<scalar(@windows);$i++){
        my $region=$windowInterval*($i+1);
        my $lineRatio="";
        my $lineTotal="";
        for(my $j=0;$j<$matchCount;$j++){
          my $total=$results->[$i]->[$j]->[0];
          my $count=$results->[$i]->[$j]->[1];
          my $ratio=sprintf("%.3f",($total>0)?($count/$total+0.00000000000001):0);
          $lineRatio.="\t$ratio";
          $lineTotal.="\t$total";
        }
        print $fh "$chr\t$pos\t$region\t1$lineRatio\n";
        print $fh "$chr\t$pos\t$region\t2$lineTotal\n";
      }
      $posCount++;
      if($posCount%10000==0){print STDERR "$posCount...\n";}
    }
    close($fh);
    print STDERR "Total position: $posCount\n";
    my ($fh2,$tmpfile2)=tempfile(DIR=>$outdir,TEMPLATE=>'sortXXXXXX',SUFFIX=>'.txt');
    close($fh2);
    print STDERR "Sorting file: $tmpfile2\n";
    system("sort -k4,4n -k3,3n -k1,1n -k2,2 $tmpfile>$tmpfile2");
    unlink($tmpfile);
    my ($fh3,$tmpfile3)=tempfile(DIR=>$outdir,TEMPLATE=>'finalXXXXXX',SUFFIX=>'.txt');
    print STDERR "Summarizing results: $tmpfile3\n";
    my $header2="#Chr\tPosition\tRegionSize(kb)\tHDR:1/2";
    foreach my $matchName(@{$matchNames}){$header2.="\t".$matchName->[0];}
    print $fh3 "$header1\n$header2\n";
    open(IN,$tmpfile2);
    while(<IN>){print $fh3 "$_";}
    close(IN);
    close($fh3);
    unlink($tmpfile2);
    system("mv $tmpfile3 $outputFile");
    print STDERR "Completed: $outputFile\n";
  }else{#DD
    my @positions=();
    my $posCount=0;
    while(!eof($reader)){
      my ($chr,$min,$max)=nextPosition($reader);
      @positions=nextTable($handler,$chr,$min,$max,@positions);
      my @hits=();
      for(my $i=0;$i<$matchCount;$i++){
        my ($label,$index1,$index2)=@{$matchNames->[$i]};
        foreach my $position(@positions){
          my ($chrom,$posit,@data)=@{$position};
          my $flag1=$data[$index1];
          my $flag2=$data[$index2];
          if($flag1==0&&$flag2==0){next;}
          if($noindel){#in case of noindel option
            if(($flag1&4)>0||($flag2&4)>0){next;}
            if(($flag1&8)>0||($flag2&8)>0){next;}
          }
          if($nolowqc){#don't include $nolowqc
            if(($flag1&32)>0||($flag1&32)>0){next;}
          }
          push(@{$hits[$i]},[$chrom,$posit,$flag1,$flag2]);
        }
      }
      my $results=[];
      for(my $i=0;$i<$matchCount;$i++){
        $results->[$i]->[0]=0;
        $results->[$i]->[1]=0;
        foreach my $hit(@{$hits[$i]}){
          my ($chrom,$posit,$flag1,$flag2)=@{$hit};
          my $f1=($flag1&3);
          my $f2=($flag2&3);
          if($targetMode eq "AR"){
            if($f1==2||$f2==2){
              $results->[$i]->[0]++;
              if($f1!=$f2){$results->[$i]->[1]++;}
            }
          }elsif($targetMode eq "DD"){
            if($f1==1||$f2==1){
              $results->[$i]->[0]++;
              if($f1!=$f2){$results->[$i]->[1]++;}
            }
          }elsif($targetMode eq "AD"){
            if($f1>0||$f2>0){
              $results->[$i]->[0]++;
              if($f1!=$f2){$results->[$i]->[1]++;}
            }
          }
        }
      }
      my $lineRatio="";
      my $lineTotal="";
      for(my $i=0;$i<$matchCount;$i++){
        my $total=$results->[$i]->[0];
        my $count=$results->[$i]->[1];
        my $ratio=sprintf("%.3f",($total>0)?($count/$total+0.00000000000001):0);
        $lineRatio.="\t$ratio";
        $lineTotal.="\t$total";
      }
      print $fh "$chr\t$min\t$max\t1$lineRatio\n";
      print $fh "$chr\t$min\t$max\t2$lineTotal\n";
      $posCount++;
      if($posCount%10000==0){print STDERR "$posCount...\n";}
    }
    close($fh);
    print STDERR "Total position: $posCount\n";
    my ($fh2,$tmpfile2)=tempfile(DIR=>$outdir,TEMPLATE=>'sortXXXXXX',SUFFIX=>'.txt');
    close($fh2);
    print STDERR "Sorting file: $tmpfile2\n";
    system("sort -k4,4n -k1,1n -k2,2 -k3,3n $tmpfile>$tmpfile2");
    unlink($tmpfile);
    my ($fh3,$tmpfile3)=tempfile(DIR=>$outdir,TEMPLATE=>'finalXXXXXX',SUFFIX=>'.txt');
    print STDERR "Summarizing results: $tmpfile3\n";
    my $header2="#Chr\tStart\tEnd\tHDR:1/2";
    foreach my $matchName(@{$matchNames}){$header2.="\t".$matchName->[0];}
    print $fh3 "$header1\n$header2\n";
    open(IN,$tmpfile2);
    while(<IN>){print $fh3 "$_";}
    close(IN);
    close($fh3);
    unlink($tmpfile2);
    system("mv $tmpfile3 $outputFile");
    print STDERR "Completed: $outputFile\n";
  }
}
############################## createLabel ##############################
sub createLabel{
  my $caseNames=shift();
  my $ctrlNames=shift();
  my $posFile=shift();
  my $windowStart=shift();
  my $windowEnd=shift();
  my $windowInterval=shift();
  my $caseCount=scalar(@{$caseNames});
  my $ctrlCount=scalar(@{$ctrlNames});
  open(IN,$posFile);
  my $label=<IN>;
  my $count=0;
  while(<IN>){$count++;}
  close(IN);
  chomp($label);
  $label=~s/\r//g;
  my $basename=basename($posFile,".txt");
  my @tokens=split(/\t/,$label);
  my $size;
  my $outFile;
  if(scalar(@tokens)>2){#chr,start,end
    $size=1;
    $outFile="${basename}_${targetMode}";
  }else{#chr,pos
    $size=int(($windowEnd-$windowStart)/$windowInterval)+1;
    $outFile="${basename}_${targetMode}_${windowStart}_${windowEnd}_${windowInterval}";
  }
  return ("$count $ctrlCount $size # 1:chr1",$outFile);
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
############################## matchNames ##############################
sub matchNames{
  my $caseNames=shift();
  my $ctrlNames=shift();
  my $caseSize=scalar(@{$caseNames});
  my $ctrlSize=scalar(@{$ctrlNames});
  my @matchNames=();
  for(my $i=0;$i<$caseSize;$i++){
    my @array=();
    my $caseName=$caseNames->[$i];
    my $caseLabel=($caseSize==1)?"Case":"Case".($i+1);
    for(my $j=0;$j<$ctrlSize;$j++){
      my $ctrlName=$ctrlNames->[$j];
      my $ctrlLabel=($ctrlSize==1)?"Ctrl":"Ctrl".($j+1);
      push(@array,["$caseLabel-$ctrlLabel",$caseName,$ctrlName]);
    }
    for(my $j=0;$j<$ctrlSize;$j++){
      my $ctrlName1=$ctrlNames->[$j];
      my $ctrlLabel1=($ctrlSize==1)?"Ctrl":"Ctrl".($j+1);
      for(my $k=$j+1;$k<$ctrlSize;$k++){
        my $ctrlName2=$ctrlNames->[$k];
        my $ctrlLabel2=($ctrlSize==1)?"Ctrl":"Ctrl".($k+1);
        push(@array,["$ctrlLabel1-$ctrlLabel2",$ctrlName1,$ctrlName2]);
      }
    }
    push(@matchNames,\@array);
  }
  return \@matchNames;
}
############################## nextPosition ##############################
sub nextPosition{
  my $reader=shift();
  my $line=<$reader>;chomp($line);$line=~s/\r//g;
  my @tokens=split(/\t/,$line);
  if($tokens[0]!~/^chr/){$tokens[0]="chr".$tokens[0];}
  elsif($tokens[0]=~/^CHR(.+)$/i){$tokens[0]="chr$1";}
  if($tokens[0] eq "chrMT"){$tokens[0]="chrM";}
  return @tokens;
}
############################## nextTable ##############################
sub nextTable{
  my @positions=@_;
  my $handler=shift(@positions);
  my $chromosome=shift(@positions);
  my $start=shift(@positions);
  my $end=shift(@positions);
  my ($reader,$next)=@{$handler};
  # Remove data out of range
  my $index=0;
  foreach my $position(@positions){
    my ($chr,$pos,@data)=@{$position};
    if($chr cmp $chromosome){$index++;}
    elsif($pos<$start){$index++;}
    else{last;}
  }
  if($index>0){@positions=splice(@positions,$index);}
  if(!defined($next)){return @positions;}
  my $endReached=0;
  my ($chr,$pos,$ref,$alt,@data)=@{$next};
  while(($chr cmp $chromosome)<0||($pos<$start)){
    if(eof($reader)){$endReached=1;last;}
    my $line=<$reader>;chomp($line);$line=~s/\r//g;
    ($chr,$pos,$ref,$alt,@data)=split(/\t/,$line);
  }
  while(($chr eq $chromosome)&&($pos<=$end)){
			push(@positions,[$chr,$pos,@data]);
			if(eof($reader)){$endReached=1;last;}
			my $line=<$reader>;chomp($line);$line=~s/\r//g;
			($chr,$pos,$ref,$alt,@data)=split(/\t/,$line);
  }
  if($endReached){$handler->[1]=undef;}
  else{$handler->[1]=[$chr,$pos,$ref,$alt,@data];}
  return @positions;
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
############################## openPosition ##############################
sub openPosition{
  my $file=shift();
  my $reader;
  if($file=~/\.g(ip)?z$/){$reader=IO::File->new("gzip -cd $file|sort|");}
  elsif($file=~/\.b(ip)?z2$/){$reader=IO::File->new("bzip2 -cd $file|sort|");}
  else{$reader=IO::File->new("sort $file|");}
  my $line=<$reader>;
  chomp($line);
  $line=~s/\r//g;
  my @tokens=split(/\t/,$line);
  my $type=(scalar(@tokens)>2)?"region":"position";
  return ($reader,$type);
}
############################## openTable ##############################
sub openTable{
  my $file=shift();
  my $matchNames=shift();
  my $reader;
  if($file=~/\.g(ip)?z$/){$reader=IO::File->new("gzip -cd $file|");}
  elsif($file=~/\.b(ip)?z2$/){$reader=IO::File->new("bzip2 -cd $file|");}
  else{$reader=IO::File->new($file);}
  my $line=<$reader>;
  chomp($line);
  $line=~s/\r//g;
  my @tokens=split(/\t/,$line);
  my $hash={};
  for(my $i=2;$i<scalar(@tokens);$i++){$hash->{$tokens[$i]}=$i-2;}
  foreach my $name(@{$matchNames}){
    $name->[1]=$hash->{$name->[1]};
    $name->[2]=$hash->{$name->[2]};
  }
  $line=<$reader>;
  chomp($line);
  $line=~s/\r//g;
  my @tokens=split(/\t/,$line);
  return [$reader,\@tokens];
}
############################## removeCaseFromCtrl ##############################
sub removeCaseFromCtrl{
  my $cases=shift();
  my $ctrls=shift();
  my $hash={};
  foreach my $case(@{$cases}){$hash->{$case}=1;}
  my @array=();
  foreach my $ctrl(@{$ctrls}){if(!exists($hash->{$ctrl})){push(@array,$ctrl);}}
  return \@array;
}
############################## retrieveCtrlNamesFromTable ##############################
sub retrieveCtrlNamesFromTable{
  my $tableFile=shift();
  my $reader=openFile($tableFile);
  my $line=<$reader>;
  close($reader);
  chomp($line);
  $line=~s/\r//g;
  my @names=split(/\t/,$line);
  shift(@names);#chr
		shift(@names);#pos
  shift(@names);#ref
  shift(@names);#alt
  return \@names;
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
