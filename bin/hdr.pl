#!/usr/bin/perl
use strict 'vars';
use Cwd;
use File::Basename;
use File::Temp qw/tempfile/;
use Getopt::Std;
use IO::File;
use vars qw($opt_d $opt_e $opt_h $opt_i $opt_l $opt_o $opt_s $opt_t);
getopts('de:hi:lo:s:t:');
my ($prgname,$prgdir,$prgsuffix)=fileparse($0);
if(defined($opt_h)||scalar(@ARGV)<4){
  print STDERR "Command: $prgname [OPTIONS] TABLE CASE CTRL POS\n";
  print STDERR "Arguments:\n";
  print STDERR " TABLE  Table file from vcftable.pl\n";
  print STDERR "  CASE  Case file/directory\n";
  print STDERR "  CTRL  Control file/directory\n";
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
  print STDERR "Author: akira.hasegawa\@riken.jp\n";
  print STDERR "Update: 2020/06/19\n";
  exit(1);
}
my $tableFile=$ARGV[0];
my @caseFiles=listFiles("\\.([vb]cf|avinput)\$",$ARGV[1]);
my @ctrlFiles=listFiles("\\.([vb]cf|avinput)\$",$ARGV[2]);
@ctrlFiles=removeCaseFromCtrl(\@caseFiles,\@ctrlFiles);
my @posFiles=listFiles("\\.te?xt\$",$ARGV[3]);
my $outdir=(defined($opt_o))?$opt_o:"out";
my $noindel=$opt_d;
my $nolowqc=$opt_l;
my $targetMode=(defined($opt_t))?lc($opt_t):"ar";
my $windowInterval=(defined($opt_i))?$opt_i:10;
my $windowStart=(defined($opt_s))?$opt_s:10;
my $windowEnd=(defined($opt_e))?$opt_e:50;
mkdir($outdir);
my ($matchNames,$caseNames,$ctrlNames)=matchFiles(\@caseFiles,\@ctrlFiles);
for(my $i=0;$i<scalar(@{$caseNames});$i++){
  my $caseName=$caseNames->[$i];
  my $matchName=$matchNames->[$i];
  mkdir("$outdir/$caseName");
  foreach my $posFile(@posFiles){
    my $basename=basename($posFile,".txt");
    my $label=createLabel($caseName,$ctrlNames,$posFile,$windowStart,$windowEnd,$windowInterval);
    my $outFile;
    if($targetMode ne "dd"){$outFile="${basename}_${targetMode}_${windowStart}_${windowEnd}_${windowInterval}";}
    else{$outFile="${basename}_${targetMode}";}
    if(defined($noindel)){$outFile.="_noindel";}
    if(defined($nolowqc)){$outFile.="_nolowqc";}
    $outFile.=".txt";
    $outFile="$outdir/$caseName/$outFile";
    calculateHDR($matchName,$targetMode,$noindel,$nolowqc,$tableFile,$posFile,$label,$outFile,$windowStart,$windowEnd,$windowInterval);
  }
}
############################## removeCaseFromCtrl ##############################
sub removeCaseFromCtrl{
  my $caseFiles=shift();
  my $ctrlFiles=shift();
  my $hash={};
  foreach my $case(@{$caseFiles}){$hash->{$case}=1;}
  my @array=();
  foreach my $ctrl(@{$ctrlFiles}){if(!exists($hash->{$ctrl})){push(@array,$ctrl);}}
  return @array;
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
  my ($fh,$tmpfile)=tempfile(DIR=>$outdir,TEMPLATE=>'hdrXXXXXX',SUFFIX=>'.txt');
  my $handler=openTable($tableFile,$matchNames);
  my ($reader,$type)=openPosition($posFile);
  my $basename=basename($tableFile);
  my $matchCount=scalar(@{$matchNames});
  if($type eq "position"){
    my @positions=();
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
              if($targetMode eq "ar"){
                if($f1==2||$f2==2){
                  $results->[$i]->[$j]->[0]++;
                  if($f1!=$f2){$results->[$i]->[$j]->[1]++;}
                }
              }elsif($targetMode eq "dd"){
                if($f1==1||$f2==1){
                  $results->[$i]->[$j]->[0]++;
                  if($f1!=$f2){$results->[$i]->[$j]->[1]++;}
                }
              }elsif($targetMode eq "ad"){
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
        print $fh "1\t$chr\t$region\t$pos$lineRatio\n";
        print $fh "2\t$chr\t$region\t$pos$lineTotal\n";
      }
    }
    close($fh);
    my ($fh2,$tmpfile2)=tempfile(DIR=>$outdir,TEMPLATE=>'sortXXXXXX',SUFFIX=>'.txt');
    close($fh2);
    system("sort -k1,1n -k2,2 -k3,3n -k4,4n $tmpfile>$tmpfile2");
    unlink($tmpfile);
    my ($fh3,$tmpfile3)=tempfile(DIR=>$outdir,TEMPLATE=>'finalXXXXXX',SUFFIX=>'.txt');
    my $header2=($targetMode eq "dd")?"#Chr\tStart\tEnd\tHDR:1/2":"#Chr\tPosition\tRegionSize(kb)\tHDR:1/2";
    foreach my $matchName(@{$matchNames}){$header2.="\t".$matchName->[0];}
    print $fh3 "$header1\n$header2\n";
    open(IN,$tmpfile2);
    while(<IN>){
      my ($label,$chr,$region,$pos,@data)=split(/\t/);
      print $fh3 "$chr\t$pos\t$region\t$label\t".join("\t",@data);
    }
    close(IN);
    close($fh3);
    unlink($tmpfile2);
    system("mv $tmpfile3 $outputFile");
  }else{
    my @positions=();
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
          if($targetMode eq "ar"){
            if($f1==2||$f2==2){
              $results->[$i]->[0]++;
              if($f1!=$f2){$results->[$i]->[1]++;}
            }
          }elsif($targetMode eq "dd"){
            if($f1==1||$f2==1){
              $results->[$i]->[0]++;
              if($f1!=$f2){$results->[$i]->[1]++;}
            }
          }elsif($targetMode eq "ad"){
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
      print $fh "1\t$chr\t$min\t$max$lineRatio\n";
      print $fh "2\t$chr\t$min\t$max$lineTotal\n";
    }
    close($fh);
    my ($fh2,$tmpfile2)=tempfile(DIR=>$outdir,TEMPLATE=>'sortXXXXXX',SUFFIX=>'.txt');
    close($fh2);
    system("sort -k1,1n -k2,2 -k3,3n $tmpfile>$tmpfile2");
    unlink($tmpfile);
    my ($fh3,$tmpfile3)=tempfile(DIR=>$outdir,TEMPLATE=>'finalXXXXXX',SUFFIX=>'.txt');
    my $header2=($targetMode eq "dd")?"#Chr\tStart\tEnd\tHDR:1/2":"#Chr\tPosition\tRegionSize(kb)\tHDR:1/2";
    foreach my $matchName(@{$matchNames}){$header2.="\t".$matchName->[0];}
    print $fh3 "$header1\n$header2\n";
    open(IN,$tmpfile2);
    while(<IN>){
      my ($label,$chr,$min,$max,@data)=split(/\t/);
      print $fh3 "$chr\t$min\t$max\t$label\t".join("\t",@data);
    }
    close(IN);
    close($fh3);
    unlink($tmpfile2);
    system("mv $tmpfile3 $outputFile");
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
  my @tokens=split(/\t/,$label);
  my $size=(scalar(@tokens)>2)?1:int(($windowEnd-$windowStart)/$windowInterval)+1;
  return "$count $ctrlCount $size # 1:chr1";
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
############################## matchFiles ##############################
sub matchFiles{
  my $caseFiles=shift();
  my $ctrlFiles=shift();
  my @caseNames=();
  foreach my $caseFile(@{$caseFiles}){
    my $basename=basename($caseFile);
    if($basename=~/^(.+)\.g\.vcf$/i){$basename=$1}
    elsif($basename=~/^(.+)\.vcf$/i){$basename=$1}
    elsif($basename=~/^(.+)\.bcf$/i){$basename=$1}
    elsif($basename=~/^(.+)\.avinput$/i){$basename=$1}
    push(@caseNames,$basename);
  }
  my @ctrlNames=();
  foreach my $ctrlFile(@{$ctrlFiles}){
    my $basename=basename($ctrlFile);
    if($basename=~/^(.+)\.g\.vcf$/i){$basename=$1}
    elsif($basename=~/^(.+)\.vcf$/i){$basename=$1}
    elsif($basename=~/^(.+)\.bcf$/i){$basename=$1}
    elsif($basename=~/^(.+)\.avinput$/i){$basename=$1}
    push(@ctrlNames,$basename);
  }
  my @matchNames=();
  for(my $i=0;$i<scalar(@caseNames);$i++){
    my @array=();
    my $caseName=$caseNames[$i];
    my $caseLabel=(scalar(@caseNames)==1)?"Case":"Case".($i+1);
    for(my $j=0;$j<scalar(@ctrlNames);$j++){
      my $ctrlName=$ctrlNames[$j];
      my $ctrlLabel=(scalar(@ctrlNames)==1)?"Ctrl":"Ctrl".($j+1);
      push(@array,["$caseLabel-$ctrlLabel",$caseName,$ctrlName]);
    }
    for(my $j=0;$j<scalar(@ctrlNames);$j++){
      my $ctrlName1=$ctrlNames[$j];
      my $ctrlLabel1=(scalar(@ctrlNames)==1)?"Ctrl":"Ctrl".($j+1);
      for(my $k=$j+1;$k<scalar(@ctrlNames);$k++){
        my $ctrlName2=$ctrlNames[$k];
        my $ctrlLabel2=(scalar(@ctrlNames)==1)?"Ctrl":"Ctrl".($k+1);
        push(@array,["$ctrlLabel1-$ctrlLabel2",$ctrlName1,$ctrlName2]);
      }
    }
    push(@matchNames,\@array);
  }
  return (\@matchNames,\@caseNames,\@ctrlNames);
}
############################## nextPosition ##############################
sub nextPosition{
  my $reader=shift();
  my $line=<$reader>;chomp($line);
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
  my ($chr,$pos,@data)=@{$next};
  while(($chr cmp $chromosome)<0||($pos<$start)){
    if(eof($reader)){$endReached=1;last;}
    my $line=<$reader>;chomp($line);
    ($chr,$pos,@data)=split(/\t/,$line);
  }
  while(($chr eq $chromosome)&&($pos<=$end)){
    push(@positions,[$chr,$pos,@data]);
    if(eof($reader)){$endReached=1;last;}
    my $line=<$reader>;chomp($line);
    ($chr,$pos,@data)=split(/\t/,$line);
  }
  if($endReached){$handler->[1]=undef;}
  else{$handler->[1]=[$chr,$pos,@data];}
  return @positions;
}
############################## openPosition ##############################
sub openPosition{
  my $file=shift();
  my $reader;
  if($file=~/\.g(ip)?z$/){$reader=IO::File->new("gzip -cd $file|");}
  elsif($file=~/\.b(ip)?z2$/){$reader=IO::File->new("bzip2 -cd $file|");}
  else{$reader=IO::File->new($file);}
  my $line=<$reader>;
  chomp($line);
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
  my @tokens=split(/\t/,$line);
  my $hash={};
  for(my $i=2;$i<scalar(@tokens);$i++){$hash->{$tokens[$i]}=$i-2;}
  foreach my $name(@{$matchNames}){
    $name->[1]=$hash->{$name->[1]};
    $name->[2]=$hash->{$name->[2]};
  }
  $line=<$reader>;
  chomp($line);
  my @tokens=split(/\t/,$line);
  return [$reader,\@tokens];
}
