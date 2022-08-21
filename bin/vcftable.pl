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
use vars qw($opt_h $opt_o $opt_q $opt_t);
getopts('ho:qt:');
############################## HELP ##############################
sub help{
  print "\n";
  print "Command: $program_name [option] VCF [VCF2 ..]\n";
  print "Arguments:\n";
  print "   VCF  variant call format (VCF or GVCF or avinput) files or under directories\n";
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
  print "Note:\n";
  print "    If you are using BCF files, please install bcftools\n";
  print "    http://samtools.github.io/bcftools/bcftools.html\n";
  print "\n";
  print "File Suffix:\n";
  print "      .vcf  VCF file\n";
  print "    .g.vcf  GVCF file\n";
  print "  .avinput  AVINPUT file\n";
  print "\n";
  print "Author: akira.hasegawa\@riken.jp\n";
  print "Update: $program_version\n";
  print "\n";
  print "Note:\n";
  print "  If only one VCF file is specified, vcftable.pl assumes VCF file has multiple entries and parse\n";
  print "\n";
}
############################## MAIN ##############################
if($ARGV[0]eq"sortsubs"){sortSubs();exit();}
elsif($ARGV[0]eq"test"){test();exit();}
elsif(defined($opt_h)||scalar(@ARGV)<1){help();exit(1);}
my $threshold=defined($opt_t)?$opt_t:50;
my $tableFile=defined($opt_o)?$opt_o:'-';
my $tmpDir=(-e "/tmp")?"/tmp":"tmp";
mkdir($tmpDir);
my @vcfFiles=listFiles("\\.([vb]cf|avinput)\$",@ARGV);
if(scalar(@vcfFiles)>1){
  my ($splitFiles,$chromosomes,$basenames)=splitInputs($tmpDir,$threshold,@vcfFiles);
  createTable($tableFile,$splitFiles,$chromosomes,$basenames);
  foreach my $chr(keys(%{$splitFiles})){foreach my $splitFile(values(%{$splitFiles->{$chr}})){unlink($splitFile);}}
}elsif(scalar(@vcfFiles)==1){
  parseMultiVCF($tableFile,$vcfFiles[0]);
}
if($tmpDir ne "/tmp"){rmdir($tmpDir);}
############################## absolutePath ##############################
sub absolutePath{
	my $path=shift();
	my $directory=dirname($path);
	my $filename=basename($path);
	return Cwd::abs_path($directory)."/$filename";
}
############################## createFile ##############################
sub createFile{
	my @lines=@_;
	my $path=shift(@lines);
	mkdirs(dirname($path));
	open(OUT,">$path");
	foreach my $line(@lines){print OUT "$line\n";}
	close(OUT);
}
############################## createTable ##############################
sub createTable{
  my $tableFile=shift();
  my $splitFiles=shift();
  my $chromosomes=shift();
  my $basenames=shift();
  my $writer=IO::File->new(">$tableFile");
  my @labels=();
  print $writer "#chromosome\tposition\t".join("\t",@{$basenames})."\n";
  foreach my $chromosome(@{$chromosomes}){
    my @readers=();
    foreach my $basename(@{$basenames}){push(@readers,openTable($splitFiles->{$chromosome}->{$basename},$threshold));}
    my $count=0;
    while(nextTable($writer,@readers)){$count++;}
    if(!defined($opt_q)){print STDERR "#$chromosome count:\t$count\n";}
  }
  close($writer);
}
############################## equals ##############################
sub equals{
	my $obj1=shift();
	my $obj2=shift();
	my $ref1=ref($obj1);
	my $ref2=ref($obj2);
	if($ref1 ne $ref2){return;}
	if($ref1 eq "ARRAY"){
		my $len1=scalar(@{$obj1});
		my $len2=scalar(@{$obj2});
		if($len1!=$len2){return;}
		for(my $i=0;$i<$len1;$i++){if(!equals($obj1->[$i],$obj2->[$i])){return;}}
		return 1;
	}elsif($ref1 eq "HASH"){
		my @keys1=keys(%{$obj1});
		my @keys2=keys(%{$obj2});
		my $len1=scalar(@keys1);
		my $len2=scalar(@keys2);
		if($len1!=$len2){return;}
		foreach my $key(@keys1){
			if(!exists($obj2->{$key})){return;}
			my $val1=$obj1->{$key};
			my $val2=$obj2->{$key};
			if(!equals($val1,$val2)){return;}
		}
		return 1;
	}
	if($obj1 eq $obj2){return 1;}
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
############################## mkdirs ##############################
sub mkdirs{
	my @directories=@_;
	foreach my $directory(@directories){
		if(-d $directory){next;}
		my @tokens=split(/[\/\\]/,$directory);
		if(($tokens[0] eq "")&&(scalar(@tokens)>1)){
			shift(@tokens);
			my $token=shift(@tokens);
			unshift(@tokens,"/$token");
		}
		my $string="";
		foreach my $token(@tokens){
			$string.=(($string eq "")?"":"/").$token;
			if(-d $string){next;}
			if(!mkdir($string)){return 0;}
		}
	}
	return 1;
}
############################## nextLine ##############################
sub nextLine{
  my $reader=shift();
  my $in=$reader->[0];
  my $prev=$reader->[1];
  if(eof($in)){
    close($in);
    $reader->[0]=undef;
    $reader->[1]=[];
  }else{
    my $line=<$in>;
    chomp($line);
    my @tokens=split(/\t/,$line);
    $reader->[1]=\@tokens;
  }
}
############################## nextTable ##############################
sub nextTable{
  my @readers=@_;
  my $writer=shift(@readers);
  my $currentChr;
  my $currentPos;
  foreach my $reader(@readers){
    if(!defined($reader->[0])){next;}
    my ($chr,$pos,$na,$ref,$alt,$ins,$del)=@{$reader->[1]};
    if(!defined($currentChr)||($chr cmp $currentChr)<0){$currentChr=$chr;$currentPos=$pos;}
    elsif(($chr eq $currentChr)&&($pos<$currentPos)){$currentPos=$pos;}
  }
  if(!defined($currentChr)){return 0;}
  my $line="$currentChr\t$currentPos";
  foreach my $reader(@readers){
    if(!defined($reader->[0])){$line.="\t0";next;}
    my ($chr,$pos,$flag)=@{$reader->[1]};
    if($chr ne $currentChr||$pos != $currentPos){$line.="\t0";next;}
    $line.="\t$flag";
    nextLine($reader);
  }
  print $writer "$line\n";
  return 1;
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
  my $reader=IO::File->new($file);
  my $line=<$reader>;
  chomp($line);
  my @tokens=split(/\t/,$line);
  return [$reader,\@tokens];
}
############################## parseAVINPUT ##############################
sub parseAVINPUT{
  my $basename=shift();
  my $threshold=shift();
  my $reader=shift();
  my $writers=shift();
  my $splitFiles=shift();
  my $previous_chr;
  while(<$reader>){
     if(/^#/){next;}#skip comment
    chomp();#remove \n
    my ($chr,$start,$end,$ref,$alt,$type,$score,$qual,$chr2,$start2,$id,$ref2,$alt2,$score2,$qual2,$info,$format,$na)=split(/\t/);
    if($chr!~/^chr/){$chr="chr$chr";}#add chr
    elsif($chr=~/^CHR(.+)$/i){$chr="chr$1";}#lowercase
    if($chr eq "chrMT"){$chr="chrM";}#handle mitochondria
    my ($a,$b)=split(/[:\|\/]/,$na);#example: 0/1,1/1,1/2
    if($a eq "."){$a=0;}#./2 => 0/2
    if($b eq "."){$b=0;}#1/. => 1/0, just in case
    my $flag=0;
    if($a==0&&$b==0){next;}#wild
    elsif($a==$b&&$a>0){$flag=2;}#homo
    else{$flag=1;}#het
    if($alt eq "-"){$alt=".";}
    elsif($alt eq "0"){$alt="";for(my $i=$start;$i<=$end;$i++){$alt.=".";}}
    if($ref eq "-"){$ref=".";}
    elsif($ref eq "0"){$ref="";for(my $i=$start;$i<=$end;$i++){$ref.=".";}}
    my $multiCount=$alt=~tr/,//;#A,B,C=>3, but haven't seen this happen
    if($multiCount<$b){$multiCount=$b;}#1/2=>2, 3/3=>3, 3/4=>4
    if($multiCount>1){$flag+=16;}#multi allelic
    my $c1=length($alt);# 'TTT'=>3
    my $c2=length($ref);# 'AA'=>2 '.'=>1
    if($c1<$c2){$flag+=4;}#deletion
    elsif($c1>$c2){$flag+=8;}#insertion
    elsif($alt eq "." && $ref ne "."){$flag+=4;}
    elsif($ref eq "." && $alt ne "."){$flag+=8;}
    if($qual eq "PASS"||$qual2 eq "PASS"){}
    elsif($qual=~/^\d+$/){if($qual<$threshold){$flag+=32;}}
    if(!defined($previous_chr)||$previous_chr ne $chr){#new chromosome
      if(exists($writers->{$previous_chr})){close($writers->{$previous_chr});}#close
      if(exists($splitFiles->{$chr}->{$basename})){#append to already existing file
        my $tmpfile=$splitFiles->{$chr}->{$basename};
        $writers->{$chr}=IO::File->new(">>$tmpfile");
      }else{#crete new temporary writer
        my ($fh,$tmpfile)=tempfile(DIR=>$tmpDir,TEMPLATE=>"$basename.$chr.XXXXXX",SUFFIX=>".txt");
        $writers->{$chr}=$fh;
        if(!exists($splitFiles->{$chr})){$splitFiles->{$chr}={};}
        $splitFiles->{$chr}->{$basename}=$tmpfile;
      }
      $previous_chr=$chr;#update chromosome
    }
    my $writer=$writers->{$chr};#current writer
    print $writer "$chr\t$start\t$flag\n";#write
  }
  if(exists($writers->{$previous_chr})){close($writers->{$previous_chr});}#close last writer
}
############################## parseGVCF ##############################
sub parseGVCF{
  my $basename=shift();
  my $threshold=shift();
  my $reader=shift();
  my $writers=shift();
  my $splitFiles=shift();
  my $previous_chr;
  while(<$reader>){
    chomp;
    if(/^#/){}
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Pt2402_BL1811.Bt07
    my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$na)=split(/\t/);
    if($chr!~/^chr/){$chr="chr$chr";}
    elsif($chr=~/^CHR(.+)$/i){$chr="chr$1";}
    if($chr eq "chrMT"){$chr="chrM";}
    my ($a,$b)=split(/[:\|\/]/,$na);
    if($a eq "."){$a=0;}
    if($b eq "."){$b=0;}
    my $flag=0;
    if($a==0&&$b==0){next;}#wild
    elsif($a==$b&&$a>0){$flag=2;}#homo
    else{$flag=1;}#het
    $alt=~s/\,?\<NON_REF\>//g;
    my $multiCount=$alt=~tr/,//;#A,B,C=>3, but haven't seen this happen
    if($multiCount<$b){$multiCount=$b;}#1/2=>2, 3/3=>3, 3/4=>4
    if($multiCount>1){$flag+=16;}#multi allelic
    my $c1=length($alt);
    my $c2=length($ref);# 'AA'=>2 '.'=>1
    if($c1<$c2){$flag+=4;}
    elsif($c1>$c2){$flag+=8;}
    elsif($alt eq "." && $ref ne "."){$flag+=4;}
    elsif($ref eq "." && $alt ne "."){$flag+=8;}
    if($qual eq "PASS"){}
    elsif($qual=~/^\d+$/){if($qual<$threshold){$flag+=32;}}
    if(!defined($previous_chr)||$previous_chr ne $chr){
      if(exists($writers->{$previous_chr})){close($writers->{$previous_chr});}
      if(exists($splitFiles->{$chr}->{$basename})){
        my $tmpfile=$splitFiles->{$chr}->{$basename};
        $writers->{$chr}=IO::File->new(">>$tmpfile");
      }else{
        my ($fh,$tmpfile)=tempfile(DIR=>$tmpDir,TEMPLATE=>"$basename.$chr.XXXXXX",SUFFIX=>".txt");
        $writers->{$chr}=$fh;
        if(!exists($splitFiles->{$chr})){$splitFiles->{$chr}={};}
        $splitFiles->{$chr}->{$basename}=$tmpfile;
      }
    }
    $previous_chr=$chr;
    my $writer=$writers->{$chr};
    print $writer "$chr\t$pos\t$flag\n";
  }
  if(exists($writers->{$previous_chr})){close($writers->{$previous_chr});}
}
############################## parseMultiVCF ##############################
sub parseMultiVCF{
  my $tableFile=shift();
  my $file=shift();
  my $reader=openFile($file);
  my @samples=();
  while(<$reader>){
    if(/^##/){next;}
    if(/^#/){
      chomp;
      my @tokens=split(/\t/);
      @samples=splice(@tokens,9);
      last;
    }
  }
  my $writer=IO::File->new(">$tableFile");
  print $writer "#chromosome\tposition\t".join("\t",@samples)."\n";
  while(<$reader>){
    chomp;
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA18939
    my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@nas)=split(/\t/);
    if($chr!~/^chr/){$chr="chr$chr";}
    elsif($chr=~/^CHR(.+)$/i){$chr="chr$1";}
    if($chr eq "chrMT"){$chr="chrM";}
    my @flags=();
    my $printFlag=0;
    foreach my $na(@nas){
      my ($a,$b)=split(/[:\|\/]/,$na);
      if($a eq "."){$a=0;}
      if($b eq "."){$b=0;}
      my $flag=0;
      if($a==0&&$b==0){}#wild
      elsif($a==$b&&$a>0){$flag=2;}#homo
      else{$flag=1;}#het
      my $multiCount=$alt=~tr/,//;#A,B,C=>3, but haven't seen this happen
      if($multiCount<$b){$multiCount=$b;}#1/2=>2, 3/3=>3, 3/4=>4
      if($multiCount>1){$flag+=16;}#multi allelic
      my $c1=length($alt);
      my $c2=length($ref);# 'AA'=>2 '.'=>1
      if($c1<$c2){$flag+=4;}
      elsif($c1>$c2){$flag+=8;}
      elsif($alt eq "." && $ref ne "."){$flag+=4;}
      elsif($ref eq "." && $alt ne "."){$flag+=8;}
      if($qual eq "PASS"){}
      elsif($qual=~/^\d+$/){if($qual<$threshold){$flag+=32;}}
      push(@flags,$flag);
      if($flag>0){$printFlag=1;}
    }
    if($printFlag){print $writer "$chr\t$pos\t".join("\t",@flags)."\n";}
  }
  close($reader);
  close($writer);
}
############################## parseVCF ##############################
sub parseVCF{
  my $basename=shift();
  my $threshold=shift();
  my $reader=shift();
  my $writers=shift();
  my $splitFiles=shift();
  my $previous_chr;
   while(<$reader>){
    if(/^#/){next;}
    chomp;
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA18939
    my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$na)=split(/\t/);
    if($chr!~/^chr/){$chr="chr$chr";}
    elsif($chr=~/^CHR(.+)$/i){$chr="chr$1";}
    if($chr eq "chrMT"){$chr="chrM";}
    my ($a,$b)=split(/[:\|\/]/,$na);
    if($a eq "."){$a=0;}
    if($b eq "."){$b=0;}
    my $flag=0;
    if($a==0&&$b==0){next;}#wild
    elsif($a==$b&&$a>0){$flag=2;}#homo
    else{$flag=1;}#het
    my $multiCount=$alt=~tr/,//;#A,B,C=>3, but haven't seen this happen
    if($multiCount<$b){$multiCount=$b;}#1/2=>2, 3/3=>3, 3/4=>4
    if($multiCount>1){$flag+=16;}#multi allelic
    my $c1=length($alt);
    my $c2=length($ref);# 'AA'=>2 '.'=>1
    if($c1<$c2){$flag+=4;}
    elsif($c1>$c2){$flag+=8;}
    elsif($alt eq "." && $ref ne "."){$flag+=4;}
    elsif($ref eq "." && $alt ne "."){$flag+=8;}
    if($qual eq "PASS"){}
    elsif($qual=~/^\d+$/){if($qual<$threshold){$flag+=32;}}
    if(!defined($previous_chr)||$previous_chr ne $chr){
      if(exists($writers->{$previous_chr})){close($writers->{$previous_chr});}
      if(exists($splitFiles->{$chr}->{$basename})){
        my $tmpfile=$splitFiles->{$chr}->{$basename};
        $writers->{$chr}=IO::File->new(">>$tmpfile");
      }else{
        my ($fh,$tmpfile)=tempfile(DIR=>$tmpDir,TEMPLATE=>"$basename.$chr.XXXXXX",SUFFIX=>".txt");
        $writers->{$chr}=$fh;
        if(!exists($splitFiles->{$chr})){$splitFiles->{$chr}={};}
        $splitFiles->{$chr}->{$basename}=$tmpfile;
      }
      $previous_chr=$chr;
    }
    my $writer=$writers->{$chr};
    print $writer "$chr\t$pos\t$flag\n";
  }
  if(exists($writers->{$previous_chr})){close($writers->{$previous_chr});}
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
############################## readText ##############################
sub readText{
	my $file=shift();
	my $text="";
	open(IN,$file);
	while(<IN>){s/\r//g;$text.=$_;}
	close(IN);
	return $text;
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
############################## splitInputs ##############################
sub splitInputs{
  my @inFiles=@_;
  my $tmpDir=shift(@inFiles);
  my $threshold=shift(@inFiles);
  my @basenames=();
  my $splitFiles={};
  my $chrHash={};
  if(!defined($opt_q)){foreach my $inFile(@inFiles){print STDERR "#Input file: $inFile\n";}}
  foreach my $inFile(@inFiles){
    my $reader;
    my $basename;
    if($inFile=~/^(.+)\.g(ip)?z$/i){$basename=$1;$reader=IO::File->new("gzip -cd $inFile|");}
    elsif($inFile=~/^(.+)\.b(ip)?z2$/i){$basename=$1;$reader=IO::File->new("bzip2 -cd $inFile|");}
    elsif($inFile=~/^(.+)\.bcf$/i){$basename=$1;$reader=IO::File->new("bcftools view $inFile|");}
    else{$basename=$inFile;$reader=IO::File->new($inFile);}
    my $type;
    if($basename=~/^(.+)\.g\.vcf$/i){$type="gvcf";$basename=$1;}
    elsif($basename=~/^(.+)\.vcf$/i){$type="vcf";$basename=$1;}
    elsif($basename=~/^(.+)\.bcf$/i){$type="vcf";$basename=$1;}
    elsif($basename=~/^(.+)\.avinput$/i){$type="avinput";$basename=$1;}
    $basename=basename($basename);
    push(@basenames,basename($basename));
    my $writers={};
    if($type eq "vcf"){parseVCF($basename,$threshold,$reader,$writers,$splitFiles);}
    elsif($type eq "gvcf"){parseGVCF($basename,$threshold,$reader,$writers,$splitFiles);}
    elsif($type eq "avinput"){parseAVINPUT($basename,$threshold,$reader,$writers,$splitFiles);}
    close($reader);
    foreach my $chr(keys(%{$writers})){
      $chrHash->{$chr}=1;
      my $tmpfile=$splitFiles->{$chr}->{$basename};
      my ($fh,$tmpfile2)=tempfile(DIR=>$tmpDir,TEMPLATE=>"$basename.$chr.sort.XXXXXX",SUFFIX=>".txt");
      close($fh);
      if(!defined($opt_q)){print STDERR "#Sorting file: $tmpfile2\n";}
      system("sort -k2,2n $tmpfile>$tmpfile2");
      unlink($tmpfile);
      $splitFiles->{$chr}->{$basename}=$tmpfile2;
    }
  }
  my @chromosomes=sort{$a cmp $b}keys(%{$chrHash});
  return ($splitFiles,\@chromosomes,\@basenames);
}
############################## splitVCFLine ##############################
sub splitVCFLine{
  my $line=shift();
  my $threshold=shift();
  my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$na)=split(/\t/,$line);
  if($chr!~/^chr/){$chr="chr$chr";}
  elsif($chr=~/^CHR(.+)$/i){$chr="chr$1";}
  if($chr eq "chrMT"){$chr="chrM";}
  my ($a,$b)=split(/[:\|\/]/,$na);
  if($a eq "."){$a=0;}
  if($b eq "."){$b=0;}
  my $flag=0;
  if($a==$b&&$a>0){$flag=2;}#homo
  elsif($a==0&&$b==0){$flag=0;}#wild
  else{$flag=1;}#het
  my $count=$alt=~tr/,//;#A,B,C A,B=>3
  if($count>0){$flag+=16;}#multi allelic
  my $c1=length($alt);
  my $c2=length($ref);# 'AA' => 2 '.'=>1
  if($c1<$c2){$flag+=4;}
  elsif($c1>$c2){$flag+=8;}
  elsif($alt eq "." && $ref ne "."){$flag+=4;}
  elsif($ref eq "." && $alt ne "."){$flag+=8;}
  if($qual eq "PASS"){}
  elsif($qual=~/^\d+$/){if($qual<$threshold){$flag+=32;}}
  return ($chr,$pos,$flag);
}
############################## test ##############################
sub test{
  mkdirs("test");
  createFile(
    "test/input1.avinput",
    "chr1	10440	10440	C	-	hom	24.08	3	chr1	10439	rs112766696	AC	A	24.08	PASS	AC=2;AF=0.077;AN=26;BaseQRankSum=-4.310e-01;ClippingRankSum=0.00;DB;DP=72;ExcessHet=0.1703;FS=0.000;InbreedingCoeff=0.0754;MLEAC=2;MLEAF=0.077;MQ=15.37;MQRankSum=-9.670e-01;QD=8.03;ReadPosRankSum=0.967;SOR=0.223	GT:AD:DP:FT:GQ:PL	1/1:1,2:3:LowDP:4:64,4,0",
    "chr1	15903	15903	-	C	hom	4945.94	3	chr1	15903	rs557514207	G	GC	4945.94	PASS	AC=40;AF=0.870;AN=46;BaseQRankSum=-7.920e-01;ClippingRankSum=0.00;DB;DP=259;ExcessHet=0.0087;FS=24.539;InbreedingCoeff=0.6275;MLEAC=42;MLEAF=0.913;MQ=7.50;MQRankSum=0.00;QD=34.24;ReadPosRankSum=0.674;SOR=3.585	GT:AD:DP:FT:GQ:PL	1/1:0,3:3:LowDP:9:107,9,0",
    "chr1	20317	20317	A	-	het	3958.34	38	chr1	20316	.	GA	G	3958.34	PASS	AC=16;AF=0.348;AN=46;BaseQRankSum=-1.090e-01;ClippingRankSum=0.00;DP=984;ExcessHet=20.9140;FS=0.000;InbreedingCoeff=-0.5274;MLEAC=16;MLEAF=0.348;MQ=7.48;MQRankSum=-9.310e-01;QD=5.54;ReadPosRankSum=-4.720e-01;SOR=0.703	GT:AD:DP:GQ:PL	0/1:23,15:38:99:320,0,587",
    "chr1	81447	81447	A	G	het	492.85	41	chr1	81447	.	A	G	492.85	PASS	AC=1;AF=0.022;AN=46;BaseQRankSum=1.42;ClippingRankSum=0.00;DP=810;ExcessHet=3.0103;FS=0.000;InbreedingCoeff=-0.0223;MLEAC=1;MLEAF=0.022;MQ=44.01;MQRankSum=-1.032e+00;QD=12.02;ReadPosRankSum=0.450;SOR=0.616	GT:AD:DP:FT:GQ:PL	0/1:21,20:41:PASS:99:530,0,571",
    "chr1	82136	82137	AA	-	het	13691.92	22	chr1	82133	.	CAAAA	CAA	13691.92	PASS	AC=11,27,6;AF=0.239,0.587,0.130;AN=46;BaseQRankSum=1.70;ClippingRankSum=0.00;DP=635;ExcessHet=3.1079;FS=1.124;InbreedingCoeff=-0.0403;MLEAC=11,25,6;MLEAF=0.239,0.543,0.130;MQ=10.37;MQRankSum=0.260;QD=26.69;ReadPosRankSum=0.805;SOR=0.878	GT:AD:DP:FT:GQ:PL	2/3:0,0,10,12:22:PASS:99:766,694,676,325,322,261,268,283,0,255",
    "chr1	95011	95011	-	G	het	1320.14	40	chr1	95011	.	T	TG	1320.14	PASS	AC=7,2;AF=0.152,0.043;AN=46;BaseQRankSum=1.10;ClippingRankSum=0.00;DP=777;ExcessHet=7.2151;FS=1.793;InbreedingCoeff=-0.2425;MLEAC=6,2;MLEAF=0.130,0.043;MQ=19.55;MQRankSum=-4.795e+00;QD=4.26;ReadPosRankSum=0.999;SOR=0.866	GT:AD:DP:GQ:PGT:PID:PL	0/1:34,6,0:40:99:0|1:94986_C_T:150,0,2303,252,2321,2573",
    "chr1	104160	104160	-	ACAC	het	11577.69	19	chr1	104160	rs372078516	A	AACAC	11577.69	PASS	AC=26,16;AF=0.565,0.348;AN=46;BaseQRankSum=-1.490e-01;ClippingRankSum=0.00;DB;DP=584;ExcessHet=3.1125;FS=50.214;InbreedingCoeff=0.1776;MLEAC=24,13;MLEAF=0.522,0.283;MQ=7.79;MQRankSum=-1.549e+00;QD=29.00;ReadPosRankSum=0.942;SOR=5.670	GT:AD:DP:FT:GQ:PL	1/2:3,8,7:19:PASS:99:581,166,257,225,0,358",
    "chr1	790136	790136	A	G	het	5687.30	10	chr1	790136	rs6696240	A	G	5687.30	PASS	AC=5,7,2,5;AF=0.109,0.152,0.043,0.109;AN=46;BaseQRankSum=0.269;ClippingRankSum=0.00;DB;DP=590;ExcessHet=51.2979;FS=7.965;InbreedingCoeff=-0.7242;MLEAC=5,7,2,5;MLEAF=0.109,0.152,0.043,0.109;MQ=13.39;MQRankSum=0.00;QD=20.10;ReadPosRankSum=0.693;SOR=0.799	GT:AD:DP:FT:GQ:PGT:PID:PL	0/2:5,0,5,0,0:10:PASS:99:.:.:307,262,605,0,327,310,262,605,327,605,262,605,327,605,605",
    "chr1	791101	791101	T	G	hom	22194.47	23	chr1	791101	rs3131980	T	G	22194.47	PASS	AC=41,5;AF=0.891,0.109;AN=46;BaseQRankSum=1.02;ClippingRankSum=0.00;DB;DP=679;ExcessHet=3.0103;FS=0.000;InbreedingCoeff=0.0000;MLEAC=41,5;MLEAF=0.891,0.109;MQ=12.26;MQRankSum=0.00;QD=34.90;ReadPosRankSum=2.02;SOR=0.526	GT:AD:DP:GQ:PGT:PID:PL	1/1:0,23,0:23:69:.:.:741,69,0,741,69,741",
    "chr1	822428	822498	CCTGGCCAGCAGATCCACCCTGTCTATACTACCTGCCTGGCCAGCAGATCCACCCTGTCTATACTACCTGA	0	hom	18262.34	25	chr1	822428	.	CCTGGCCAGCAGATCCACCCTGTCTATACTACCTGCCTGGCCAGCAGATCCACCCTGTCTATACTACCTGA	0	18262.34	PASS	AC=1,37;AF=0.022,0.804;AN=46;BaseQRankSum=0.00;ClippingRankSum=0.00;DP=838;ExcessHet=6.1884;FS=2.276;InbreedingCoeff=-0.2105;MLEAC=1,37;MLEAF=0.022,0.804;MQ=8.79;MQRankSum=0.00;QD=29.74;ReadPosRankSum=-7.330e-01;SOR=0.870	GT:AD:DP:GQ:PL	2/2:0,0,25:25:78:1022,1022,1022,78,78,0",
    "chr1	840411	840411	A	-	hom	14651.44	22	chr1	840409	rs755461528	TAA	TA	14651.44	PASS	AC=9,8,29;AF=0.196,0.174,0.630;AN=46;BaseQRankSum=-4.920e-01;ClippingRankSum=0.00;DB;DP=775;ExcessHet=3.0103;FS=0.000;InbreedingCoeff=0.0007;MLEAC=8,7,29;MLEAF=0.174,0.152,0.630;MQ=11.74;MQRankSum=-3.780e-01;QD=21.33;ReadPosRankSum=0.406;SOR=0.722	GT:AD:DP:GQ:PL	3/3:0,0,0,22:22:65:463,463,463,463,463,463,65,65,65,0",
    "chr1	984611	984611	-	T	hom	19540.64	12	chr1	984611	.	CTTAT	CTTTAT	19540.64	PASS	AC=3,32,1,2,3,2;AF=0.065,0.696,0.022,0.043,0.065,0.043;AN=46;BaseQRankSum=1.37;ClippingRankSum=0.00;DP=526;ExcessHet=3.3099;FS=1.943;InbreedingCoeff=-0.0688;MLEAC=3,32,1,2,3,2;MLEAF=0.065,0.696,0.022,0.043,0.065,0.043;MQ=12.55;MQRankSum=0.00;QD=26.25;ReadPosRankSum=0.199;SOR=0.647	GT:AD:DP:FT:GQ:PGT:PID:PL	2/2:0,0,12,0,0,0,0:12:PASS:39:1|1:984611_C_CT:584,585,585,39,39,0,585,585,39,585,585,585,39,585,585,585,585,39,585,585,585,585,585,39,585,585,585,585",

    "chr1	1010440	1010440	C	-	hom	24.08	3	chr1	1010439	rs112766696	AC	A	24.08	10	AC=2;AF=0.077;AN=26;BaseQRankSum=-4.310e-01;ClippingRankSum=0.00;DB;DP=72;ExcessHet=0.1703;FS=0.000;InbreedingCoeff=0.0754;MLEAC=2;MLEAF=0.077;MQ=15.37;MQRankSum=-9.670e-01;QD=8.03;ReadPosRankSum=0.967;SOR=0.223	GT:AD:DP:FT:GQ:PL	1/1:1,2:3:LowDP:4:64,4,0",
    "chr1	1015903	1015903	-	C	hom	4945.94	3	chr1	1015903	rs557514207	G	GC	4945.94	11	AC=40;AF=0.870;AN=46;BaseQRankSum=-7.920e-01;ClippingRankSum=0.00;DB;DP=259;ExcessHet=0.0087;FS=24.539;InbreedingCoeff=0.6275;MLEAC=42;MLEAF=0.913;MQ=7.50;MQRankSum=0.00;QD=34.24;ReadPosRankSum=0.674;SOR=3.585	GT:AD:DP:FT:GQ:PL	1/1:0,3:3:LowDP:9:107,9,0",
    "chr1	1020317	1020317	A	-	het	3958.34	38	chr1	1020316	.	GA	G	3958.34	12	AC=16;AF=0.348;AN=46;BaseQRankSum=-1.090e-01;ClippingRankSum=0.00;DP=984;ExcessHet=20.9140;FS=0.000;InbreedingCoeff=-0.5274;MLEAC=16;MLEAF=0.348;MQ=7.48;MQRankSum=-9.310e-01;QD=5.54;ReadPosRankSum=-4.720e-01;SOR=0.703	GT:AD:DP:GQ:PL	0/1:23,15:38:99:320,0,587",
    "chr1	1081447	1081447	A	G	het	492.85	41	chr1	1081447	.	A	G	492.85	13	AC=1;AF=0.022;AN=46;BaseQRankSum=1.42;ClippingRankSum=0.00;DP=810;ExcessHet=3.0103;FS=0.000;InbreedingCoeff=-0.0223;MLEAC=1;MLEAF=0.022;MQ=44.01;MQRankSum=-1.032e+00;QD=12.02;ReadPosRankSum=0.450;SOR=0.616	GT:AD:DP:FT:GQ:PL	0/1:21,20:41:PASS:99:530,0,571",
    "chr1	1082136	1082137	AA	-	het	13691.92	22	chr1	1082133	.	CAAAA	CAA	13691.92	14	AC=11,27,6;AF=0.239,0.587,0.130;AN=46;BaseQRankSum=1.70;ClippingRankSum=0.00;DP=635;ExcessHet=3.1079;FS=1.124;InbreedingCoeff=-0.0403;MLEAC=11,25,6;MLEAF=0.239,0.543,0.130;MQ=10.37;MQRankSum=0.260;QD=26.69;ReadPosRankSum=0.805;SOR=0.878	GT:AD:DP:FT:GQ:PL	2/3:0,0,10,12:22:PASS:99:766,694,676,325,322,261,268,283,0,255",
    "chr1	1095011	1095011	-	G	het	1320.14	40	chr1	1095011	.	T	TG	1320.14	15	AC=7,2;AF=0.152,0.043;AN=46;BaseQRankSum=1.10;ClippingRankSum=0.00;DP=777;ExcessHet=7.2151;FS=1.793;InbreedingCoeff=-0.2425;MLEAC=6,2;MLEAF=0.130,0.043;MQ=19.55;MQRankSum=-4.795e+00;QD=4.26;ReadPosRankSum=0.999;SOR=0.866	GT:AD:DP:GQ:PGT:PID:PL	0/1:34,6,0:40:99:0|1:94986_C_T:150,0,2303,252,2321,2573",
    "chr1	10104160	10104160	-	ACAC	het	11577.69	19	chr1	10104160	rs372078516	A	AACAC	11577.69	16	AC=26,16;AF=0.565,0.348;AN=46;BaseQRankSum=-1.490e-01;ClippingRankSum=0.00;DB;DP=584;ExcessHet=3.1125;FS=50.214;InbreedingCoeff=0.1776;MLEAC=24,13;MLEAF=0.522,0.283;MQ=7.79;MQRankSum=-1.549e+00;QD=29.00;ReadPosRankSum=0.942;SOR=5.670	GT:AD:DP:FT:GQ:PL	1/2:3,8,7:19:PASS:99:581,166,257,225,0,358",
    "chr1	10790136	10790136	A	G	het	5687.30	10	chr1	10790136	rs6696240	A	G	5687.30	17	AC=5,7,2,5;AF=0.109,0.152,0.043,0.109;AN=46;BaseQRankSum=0.269;ClippingRankSum=0.00;DB;DP=590;ExcessHet=51.2979;FS=7.965;InbreedingCoeff=-0.7242;MLEAC=5,7,2,5;MLEAF=0.109,0.152,0.043,0.109;MQ=13.39;MQRankSum=0.00;QD=20.10;ReadPosRankSum=0.693;SOR=0.799	GT:AD:DP:FT:GQ:PGT:PID:PL	0/2:5,0,5,0,0:10:PASS:99:.:.:307,262,605,0,327,310,262,605,327,605,262,605,327,605,605",
    "chr1	10791101	10791101	T	G	hom	22194.47	23	chr1	10791101	rs3131980	T	G	22194.47	18	AC=41,5;AF=0.891,0.109;AN=46;BaseQRankSum=1.02;ClippingRankSum=0.00;DB;DP=679;ExcessHet=3.0103;FS=0.000;InbreedingCoeff=0.0000;MLEAC=41,5;MLEAF=0.891,0.109;MQ=12.26;MQRankSum=0.00;QD=34.90;ReadPosRankSum=2.02;SOR=0.526	GT:AD:DP:GQ:PGT:PID:PL	1/1:0,23,0:23:69:.:.:741,69,0,741,69,741",
    "chr1	10822428	10822498	CCTGGCCAGCAGATCCACCCTGTCTATACTACCTGCCTGGCCAGCAGATCCACCCTGTCTATACTACCTGA	0	hom	18262.34	25	chr1	10822428	.	CCTGGCCAGCAGATCCACCCTGTCTATACTACCTGCCTGGCCAGCAGATCCACCCTGTCTATACTACCTGA	0	18262.34	19	AC=1,37;AF=0.022,0.804;AN=46;BaseQRankSum=0.00;ClippingRankSum=0.00;DP=838;ExcessHet=6.1884;FS=2.276;InbreedingCoeff=-0.2105;MLEAC=1,37;MLEAF=0.022,0.804;MQ=8.79;MQRankSum=0.00;QD=29.74;ReadPosRankSum=-7.330e-01;SOR=0.870	GT:AD:DP:GQ:PL	2/2:0,0,25:25:78:1022,1022,1022,78,78,0",
    "chr1	10840411	10840411	A	-	hom	14651.44	22	chr1	10840409	rs755461528	TAA	TA	14651.44	20	AC=9,8,29;AF=0.196,0.174,0.630;AN=46;BaseQRankSum=-4.920e-01;ClippingRankSum=0.00;DB;DP=775;ExcessHet=3.0103;FS=0.000;InbreedingCoeff=0.0007;MLEAC=8,7,29;MLEAF=0.174,0.152,0.630;MQ=11.74;MQRankSum=-3.780e-01;QD=21.33;ReadPosRankSum=0.406;SOR=0.722	GT:AD:DP:GQ:PL	3/3:0,0,0,22:22:65:463,463,463,463,463,463,65,65,65,0",
    "chr1	10984611	10984611	-	T	hom	19540.64	12	chr1	10984611	.	CTTAT	CTTTAT	19540.64	21	AC=3,32,1,2,3,2;AF=0.065,0.696,0.022,0.043,0.065,0.043;AN=46;BaseQRankSum=1.37;ClippingRankSum=0.00;DP=526;ExcessHet=3.3099;FS=1.943;InbreedingCoeff=-0.0688;MLEAC=3,32,1,2,3,2;MLEAF=0.065,0.696,0.022,0.043,0.065,0.043;MQ=12.55;MQRankSum=0.00;QD=26.25;ReadPosRankSum=0.199;SOR=0.647	GT:AD:DP:FT:GQ:PGT:PID:PL	2/2:0,0,12,0,0,0,0:12:PASS:39:1|1:984611_C_CT:584,585,585,39,39,0,585,585,39,585,585,585,39,585,585,585,585,39,585,585,585,585,585,39,585,585,585,585",
  );
  testCommand(
    "perl $program_directory/vcftable.pl -q test/input1.avinput",
    "#chromosome	position	input1",
    "chr1	10440	6",
    "chr1	15903	10",  
    "chr1	20317	5",
    "chr1	81447	1",
    "chr1	82136	21",
    "chr1	95011	9",
    "chr1	104160	25",
    "chr1	790136	17",
    "chr1	791101	2",
    "chr1	822428	18",
    "chr1	840411	22",
    "chr1	984611	26",
    "chr1	1010440	38",
    "chr1	1015903	42",
    "chr1	1020317	37",
    "chr1	1081447	33",
    "chr1	1082136	53",
    "chr1	1095011	41",
    "chr1	10104160	57",
    "chr1	10790136	49",
    "chr1	10791101	34",
    "chr1	10822428	50",
    "chr1	10840411	54",
    "chr1	10984611	58",
  );
  unlink("test/input1.avinput");

  createFile(
    "test/input1.vcf",
    "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA18939",
    "1	887560	var_452	A	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	0/1:3,9:8.37:-43.31,-8.99,-17.36",
    "1	948921	var_803	T	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=134.03	GT:COV:GT_CONF:GL	1/1:0,13:9.01:-64.60,-16.09,-7.08",
    "1	10887560	var_452	AT	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	0/1:3,9:8.37:-43.31,-8.99,-17.36",
    "1	10948921	var_803	T	CA	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=134.03	GT:COV:GT_CONF:GL	1/1:0,13:9.01:-64.60,-16.09,-7.08",
    "1	20887560	var_452	AT	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	1/3:3,9:8.37:-43.31,-8.99,-17.36",
    "1	20948921	var_803	T	CA	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=134.03	GT:COV:GT_CONF:GL	2/2:0,13:9.01:-64.60,-16.09,-7.08",
  );
  testCommand(
    "perl $program_directory/vcftable.pl -q test/input1.vcf",
    "#chromosome	position	input1",
    "chr1	887560	1",
    "chr1	948921	2",
    "chr1	10887560	5",
    "chr1	10948921	10",
    "chr1	20887560	21",
    "chr1	20948921	26",
  );
  unlink("test/input1.vcf");
  system(
    "perl $program_directory/vcftable.pl -q testdata/input/NA18939_v2.vcf testdata/input/NA18940_v2.vcf testdata/input/NA18941_v2.vcf testdata/input/NA18942_v2.vcf > test/table.txt");
  rmdir("test");
  testCommand("diff test/table.txt testdata/table.txt","");
}
############################## testCommand ##############################
sub testCommand{
	my @values=@_;
	my $command=shift(@values);
	my $value2=join("\n",@values);
	my ($writer,$file)=tempfile();
	close($writer);
	if(system("$command > $file")){
		print STDERR ">$command\n";
		print STDERR "Command failed...\n";
		return 1;
	}
	my $value1=readText($file);
	chomp($value1);
	if($value2 eq""){if($value1 eq""){return 0;}}
	if($value1 eq $value2){return 0;}
	print STDERR ">$command\n";
	print STDERR "$value1\n";
	print STDERR "$value2\n";
}
############################## testSub ##############################
sub testSub{
	my $command=shift();
	my $value2=shift();
	my $value1=eval($command);
	if(equals($value1,$value2)){return 0;}
  print STDERR ">$command\n";
	printTable($value1);
	printTable($value2);
}
