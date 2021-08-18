#!/usr/bin/perl
use strict 'vars';
use Cwd;
use File::Basename;
use File::Temp qw/tempfile/;
use IO::File;
use Getopt::Std;
use vars qw($opt_h $opt_o $opt_t);
getopts('ho:t:');
my ($prgname,$prgdir,$prgsuffix)=fileparse($0);
if(defined($opt_h)||scalar(@ARGV)<1){help();exit(1);}
my $threshold=defined($opt_t)?$opt_t:50;
my $tableFile=defined($opt_o)?$opt_o:'-';
my $tmpDir=(-e "/tmp")?"/tmp":"tmp";
my @vcfFiles=listFiles("\\.([vb]cf|avinput)\$",@ARGV);
mkdir($tmpDir);
my ($splitFiles,$chromosomes,$basenames)=splitInputs($tmpDir,$threshold,@vcfFiles);
createTable($tableFile,$splitFiles,$chromosomes,$basenames);
foreach my $chr(keys(%{$splitFiles})){foreach my $splitFile(values(%{$splitFiles->{$chr}})){unlink($splitFile);}}
rmdir($tmpDir);
############################## help ##############################
sub help{
  print STDERR "\n";
  print STDERR "Command: $prgname [option] VCF [VCF2 ..]\n";
  print STDERR "Arguments:\n";
  print STDERR "   VCF  variant call format (VCF or GVCF) or avinput files or under directories\n";
  print STDERR "Options:\n";
  print STDERR "    -o  Output file (default='STDOUT')\n";
  print STDERR "    -t  Quality value threshold for low quality (default<50)\n";
  print STDERR "Flag:\n";
  print STDERR "     0  wild\n";
  print STDERR "     1  hetero\n";
  print STDERR "     2  homo\n";
  print STDERR "     4  deletion\n";
  print STDERR "     8  insertion\n";
  print STDERR "    16  multi allelic (column5.count(',')>0)\n";
  print STDERR "    32  low quality (column6<QV50)\n";
  print STDERR "Note:\n";
  print STDERR "    If you are using BCF files, please install bcftools\n";
  print STDERR "    http://samtools.github.io/bcftools/bcftools.html\n";
  print STDERR "\n";
  print STDERR "File Suffix:\n";
  print STDERR "      .vcf  VCF file\n";
  print STDERR "    .g.vcf  GVCF file\n";
  print STDERR "  .avinput  AVINPUT file\n";
  print STDERR "\n";
  print STDERR "Author: akira.hasegawa\@riken.jp\n";
  print STDERR "Update: 2020/08/08\n";
  print STDERR "\n";
}
############################## absolutePath ##############################
sub absolutePath{
	my $path=shift();
	my $directory=dirname($path);
	my $filename=basename($path);
	return Cwd::abs_path($directory)."/$filename";
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
    print STDERR "$chromosome\t$count\n";
  }
  close($writer);
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
  my $line;
   while(<$reader>){
     chomp;
     if(/^#/){}
     else{$line=$_;last;}
  }
  my $previous_chr;
  while(defined($line)){
    chomp($line);
    my ($chr,$start,$end,$ref,$alt,$type,$score,$qual,$chr2,$start2,$id,$ref2,$alt2,$score2,$qual2,$info,$format,$na)=split(/\t/,$line);
    if($chr!~/^chr/){$chr="chr$chr";}
    elsif($chr=~/^CHR(.+)$/i){$chr="chr$1";}
    if($chr eq "chrMT"){$chr="chrM";}
    my ($a,$b)=split(/[:\|\/]/,$na);
    if($a eq "."){$a=0;}
    if($b eq "."){$b=0;}
    if($a==0&&$b==0){$line=<$reader>;next;}
    my $flag=0;
    if($a==$b&&$a>0){$flag=2;}#homo
    elsif($a==0&&$b==0){$flag=0;}#wild
    else{$flag=1;}#het
    if($alt eq "-"){$alt=".";}
    elsif($alt eq "0"){$alt="";for(my $i=$start;$i<=$end;$i++){$alt.=".";}}
    if($ref eq "-"){$ref=".";}
    elsif($ref eq "0"){$ref="";for(my $i=$start;$i<=$end;$i++){$ref.=".";}}
    my $count=$alt=~tr/,//;#A,B,C A,B=>3
    if($a>1&&$count<$a){$count=$a;}
    if($b>1&&$count<$b){$count=$b;}
    if($count>0){$flag+=16;}#multi alleric
    else{#others
      my $c1=length($alt);
      my $c2=length($ref);# 'AA'=>2 '.'=>1
      if($c1<$c2){$flag+=4;}
      elsif($c1>$c2){$flag+=8;}
      elsif($alt eq "." && $ref ne "."){$flag+=4;}
      elsif($ref eq "." && $alt ne "."){$flag+=8;}
    }
    if($qual eq "PASS"||$qual2 eq "PASS"){}
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
    print $writer "$chr\t$start\t$flag\n";
    $line=<$reader>;
  }
  if(exists($writers->{$previous_chr})){close($writers->{$previous_chr});}
}
############################## parseGVCF ##############################
sub parseGVCF{
  my $basename=shift();
  my $threshold=shift();
  my $reader=shift();
  my $writers=shift();
  my $splitFiles=shift();
  my $line;
  while(<$reader>){
    chomp;
    if(/^#/){}
    else{$line=$_;last;}
  }
  my $previous_chr;
  while(defined($line)){
    chomp($line);
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Pt2402_BL1811.Bt07
    my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$na)=split(/\t/,$line);
    if($chr!~/^chr/){$chr="chr$chr";}
    elsif($chr=~/^CHR(.+)$/i){$chr="chr$1";}
    if($chr eq "chrMT"){$chr="chrM";}
    my ($a,$b)=split(/[:\|\/]/,$na);
    if($a eq "."){$a=0;}
    if($b eq "."){$b=0;}
    if($a==0&&$b==0){$line=<$reader>;next;}
    my $flag=0;
    if($a==$b&&$a>0){$flag=2;}#homo
    elsif($a==0&&$b==0){$flag=0;}#wild
    else{$flag=1;}#het
    $alt=~s/\,?\<NON_REF\>//g;
    my $count=$alt=~tr/,//;#A,B,C A,B=>3
    if($count>0){$flag+=16;}#multi alleric
    else{#others
      my $c1=length($alt);
      my $c2=length($ref);# 'AA'=>2 '.'=>1
      if($c1<$c2){$flag+=4;}
      elsif($c1>$c2){$flag+=8;}
      elsif($alt eq "." && $ref ne "."){$flag+=4;}
      elsif($ref eq "." && $alt ne "."){$flag+=8;}
    }
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
    $line=<$reader>;
  }
  if(exists($writers->{$previous_chr})){close($writers->{$previous_chr});}
}
############################## parseVCF ##############################
sub parseVCF{
  my $basename=shift();
  my $threshold=shift();
  my $reader=shift();
  my $writers=shift();
  my $splitFiles=shift();
  my $line;
   while(<$reader>){
     chomp;
     if(/^#/){}
     else{$line=$_;last;}
  }
  my $previous_chr;
  while(defined($line)){
    chomp($line);
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA18939
    my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$na)=split(/\t/,$line);
    if($chr!~/^chr/){$chr="chr$chr";}
    elsif($chr=~/^CHR(.+)$/i){$chr="chr$1";}
    if($chr eq "chrMT"){$chr="chrM";}
    my ($a,$b)=split(/[:\|\/]/,$na);
    if($a eq "."){$a=0;}
    if($b eq "."){$b=0;}
    if($a==0&&$b==0){$line=<$reader>;next;}
    my $flag=0;
    if($a==$b&&$a>0){$flag=2;}#homo
    elsif($a==0&&$b==0){$flag=0;}#wild
    else{$flag=1;}#het
    my $count=$alt=~tr/,//;#A,B,C A,B=>3
    if($count>0){$flag+=16;}#multi alleric
    else{#others
      my $c1=length($alt);
      my $c2=length($ref);# 'AA'=>2 '.'=>1
      if($c1<$c2){$flag+=4;}
      elsif($c1>$c2){$flag+=8;}
      elsif($alt eq "." && $ref ne "."){$flag+=4;}
      elsif($ref eq "." && $alt ne "."){$flag+=8;}
    }
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
    $line=<$reader>;
  }
  if(exists($writers->{$previous_chr})){close($writers->{$previous_chr});}
}
############################## splitInputs ##############################
sub splitInputs{
  my @inFiles=@_;
  my $tmpDir=shift(@inFiles);
  my $threshold=shift(@inFiles);
  my @basenames=();
  my $splitFiles={};
  my $chrHash={};
  foreach my $inFile(@inFiles){
    print STDERR "Input file: $inFile\n";
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
      print STDERR "Sorting file: $tmpfile2\n";
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
  if($count>0){$flag+=16;}#multi alleric
  else{#others
    my $c1=length($alt);
    my $c2=length($ref);# 'AA' => 2 '.'=>1
    if($c1<$c2){$flag+=4;}
    elsif($c1>$c2){$flag+=8;}
    elsif($alt eq "." && $ref ne "."){$flag+=4;}
    elsif($ref eq "." && $alt ne "."){$flag+=8;}
  }
  if($qual eq "PASS"){}
  elsif($qual=~/^\d+$/){if($qual<$threshold){$flag+=32;}}
  return ($chr,$pos,$flag);
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