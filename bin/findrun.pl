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
use vars qw($opt_b $opt_d $opt_f $opt_h $opt_m $opt_o $opt_p $opt_r $opt_R $opt_s);
getopts('bdfhm:o:p:r:R:s:');
############################## HELP ##############################
sub help{
	print STDERR "\n";
	print STDERR "Command: $program_name [option] TABLE CASE\n";
	print STDERR "Arguments:\n";
	print STDERR " TABLE  Table from vcftable.pl\n";
	print STDERR "  CASE  Case file/directory\n";
	print STDERR "Options:\n";
	print STDERR "    -b  Output in (b)ed format\n";
	print STDERR "    -c  cytoband file for removing telomere\n";
	print STDERR "    -d  Exclude insertion/(d)eletion\n";
	print STDERR "    -f  (f)ull option\n";
	print STDERR "    -h  Display this (h)elp\n";
	print STDERR "    -m  hom/het/dup (m)ode stretch (default=hom)\n";
	print STDERR "    -o  (o)utput directory (default='out')\n";
	print STDERR "    -p  (p)ickup number (default>=1)\n";
	print STDERR "    -r  Minimum (r)egion size (default=1,000,000bp)\n";
	print STDERR "    -R  Maximum (R)egion size\n";
	print STDERR "    -s  (s)kip count (default=0)\n";
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
	print STDERR "\n";
}
############################## MAIN ##############################
if($ARGV[0]eq"sortsubs"){sortSubs();exit();}
elsif($ARGV[0]eq"test"){test();exit();}
elsif(defined($opt_h)||scalar(@ARGV)<2){help();exit();}
my @caseNames=@ARGV;
my $filepath=shift(@caseNames);
my $noindel=$opt_d;
my $stretchMode=defined($opt_m)?lc($opt_m):"hom";
my $outDir=defined($opt_o)?$opt_o:'out';
my $pickupNumber=defined($opt_p)?$opt_p:1;
my $minRegionSize=defined($opt_r)?$opt_r:1000000;
my $maxRegionSize=defined($opt_R)?$opt_R:undef;
my $skipCount=defined($opt_s)?$opt_s:0;
my $bedformat=$opt_b;
my $fullMode=$opt_f;
my $outdir=(defined($opt_o))?$opt_o:"out";
if(scalar(@caseNames)==0){
	print STDERR "\n";
	print STDERR "ERROR  No VCF files were found.\n";
	print STDERR "ERROR  Make sure you specify correct CASE file(s)/directory in the command line.\n";
	help();
	exit(1);
}
if(!-e $filepath){
	print STDERR "\n";
	print STDERR "ERROR  Table file '$filepath' doesn't exist.\n";
	print STDERR "ERROR  Make sure you specify correct TABLE in the command line.\n";
	help();
	exit(1);
}
mkdir($outDir);
#go through files
my $regionCount=0;
my @tableFiles=();
if(-d $filepath){@tableFiles=getFiles($filepath);}
else{@tableFiles=($filepath);}
my $writers;
foreach my $tableFile(@tableFiles){
	my ($handler,$labels)=openTable($tableFile);
	if(!defined($writers)){
		$writers=matchIndex(\@caseNames,$labels,$outDir,$stretchMode,$noindel,$pickupNumber,$minRegionSize,$maxRegionSize,$skipCount,$fullMode,$bedformat);
	}
	while(!eof($handler->[0])){
		$regionCount+=checkRegion(nextTable($handler,$writers,$minRegionSize,$noindel),$writers,$stretchMode,$pickupNumber,$minRegionSize,$maxRegionSize,$skipCount,$fullMode,$bedformat);
		if($regionCount%10000==0){print STDERR "$regionCount...\n";}
	}
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
############################## checkRegion ##############################
# centromereは後の工程で省いている。途中のfindrun/HDRでは除いていない。pickupNumberを指定して行っていた
#1,000,000 big deletion　これまで　それぞれのSNPの重要性が下がる
#500-2000 small deletion IN/DEL　これから　一つのSNPが重要になってくる　目視で
#hom =     INS DEL 個別のチェックこれが必要？
#hom =         DEL 個別のチェックこれが必要？
#hom =     INS     個別のチェックこれも必要？
#hom = SNP/INS/DEL noindel=F　現在の設定
#hom = SNP         noindel=T　現在の設定（これを使っていた）
#hom hom hom hom het
#|-------------|
#het het het het hom
#|-------------|
sub checkRegion{
	my $positions=shift();
	my $writers=shift();
	my $stretchMode=shift();
	my $pickupNumber=shift();
	my $minRegionSize=shift();
	my $maxRegionSize=shift();
	my $skipCount=shift();
	my $fullMode=shift();
	my $bedFormat=shift();
	my $startFlag=($stretchMode eq "hom")?2:1;
	my $endFlag=($stretchMode eq "hom")?1:2;
	if($stretchMode eq "dup"){
		$startFlag=512;
		$endFlag=2;
	}
	my $size=scalar(@{$positions});
	my $count=0;
	my $index=0;
	foreach my $handler(@{$writers}){
		my $label=$handler->[0];
		my $writer=$handler->[2];
		my $lastChromosome;
		my $lastEnd;
		for(my $i=0;$i<$size;$i++){
			printTable($positions->[$i]);
			my $chromosome=$positions->[$i]->[0];
			my $start=$positions->[$i]->[1];
			my $end=$positions->[$i]->[2];
			my $flag=$positions->[$i]->[2+$index];
			if(($flag&$startFlag)>0){#start search
				my $string=$flag;
				my $end=$start;
				my $match=0;
				my $mismatch=0;
				for(my $j=$i+1;$j<$size;$j++){
					my $length=$end-$start+1;
					my $flag2=$positions->[$j]->[2+$index];
					if(($flag2&$startFlag)>0){#continue search
						$end=$positions->[$j]->[1];
						$match++;
					}elsif(($flag2&$endFlag)>0){#end search
						if($mismatch<$skipCount){#search more
							if($flag2>0){$string.=" $flag2";}
							$mismatch++;
							next;
						}
						if($length<$minRegionSize){last;}
						if(defined($maxRegionSize)&&$length>$maxRegionSize){last;}
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
						$count++;
						last;
					}
					if($flag2>0){$string.=" $flag2";}
				}
			}
		}
		$index+=1;
	}
	return $count;
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
############################## getFiles ##############################
sub getFiles{
	my $directory=shift();
	my $grep=shift();
	my @files=();
	opendir(DIR,$directory);
	if(ref($grep)eq"ARRAY"){
		foreach my $file(readdir(DIR)){
			if($file=~/^\./){next;}
			if($file eq""){next;}
			my $path="$directory/$file";
			if(-d $path){next;}
			my $hit=0;
			foreach my $g(@{$grep}){if($path=~/$g/){$hit=1;}}
			if($hit){push(@files,$path);}
		}
	}else{
		foreach my $file(readdir(DIR)){
			if($file=~/^\./){next;}
			if($file eq""){next;}
			my $path="$directory/$file";
			if(-d $path){next;}
			if(defined($grep)&&$path!~/$grep/){next;}
			push(@files,$path);
		}
	}
	closedir(DIR);
  return sort{$a cmp $b}@files;
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
############################## matchIndex ##############################
sub matchIndex{
	my $files=shift();
	my $labels=shift();
	my $outDir=shift();
	my $stretchMode=shift();
	my $noindel=shift();
	my $pickupNumber=shift();
	my $minRegionSize=shift();
	my $maxRegionSize=shift();
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
			my $filename="$label.${stretchMode}_pick${pickupNumber}_skip${skipCount}";
			if(defined($maxRegionSize)){$filename.="_min${minRegionSize}_max${maxRegionSize}";}
			else{$filename.="_min${minRegionSize}";}
			if(defined($noindel)){$filename.="_noindel";}
			if(defined($fullMode)){$filename.="_full";}
			if(defined($bedformat)){$filename.=".bed";}
			else{$filename.=".txt";}
			my $writer=IO::File->new(">$outDir/$filename");
			if(defined($bedformat)){}
			elsif(defined($fullMode)){print $writer "#Chr\tStart\tEnd\tLength\tHit\tUnhit\tGenotype\n";}
			else{print $writer "#Chr\tStart\tEnd\n";}
			push(@array,[$label,$i,$writer]);
	}
	return \@array;
}
############################## nextTable ##############################
sub nextTable{
	my $handler=shift();
	my $writers=shift();
	my $minRegionSize=shift();
	my $noindel=shift();
	my @indeces=();
	#These lines are added to reduce memory usage
	foreach my $handler(@{$writers}){push(@indeces,$handler->[1]);}
	my $reader=$handler->[0];
	my $next=$handler->[1];
	if(!defined($next)){return ();}
	my @data=@{$next};
	my $chromosome=shift(@data);
	my $position=shift(@data);
	my @positions=([$chromosome,$position,@data]);
	my $endReached=0;
	my $previousPoint;
	while(<$reader>){
			chomp;s/\r//g;
			my ($chr,$pos,$ref,$alt,@data)=split(/\t/);
			my @temp=();
			foreach my $index(@indeces){push(@temp,$data[$index]);}
			if($chr ne $chromosome){$handler->[1]=[$chr,$pos,@temp];last;}
			if($pos-$position>$minRegionSize){$handler->[1]=[$chr,$pos,@temp];last;}
			if(defined($noindel)){foreach my $d(@data){if(($d&4)>0||($d&8)>0){$d=0;}}}#take care of indel
			push(@positions,[$chr,$pos,@temp]);
			$previousPoint=$pos;
	}
	if(eof($reader)){close($reader);$handler->[1]=undef;}
	return \@positions;
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
	$line=~s/\r//g;
	my @names=split(/\t/,$line);
	shift(@names);#chr
	shift(@names);#pos
	shift(@names);#ref
	shift(@names);#alt
	$line=<$reader>;
	chomp($line);
	$line=~s/\r//g;
	my @tokens=split(/\t/,$line);
	return ([$reader,\@tokens],\@names);
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
############################## test ##############################
sub test{
	createFile("test/test.txt","#chromosome	position	ref	alt	NA18939	NA18940	NA18941	NA18942",
	"chr1	1	A	C	1	1	0	0",
	"chr1	2	T	C	2	0	0	0",
	"chr1	3	G	A	2	0	2	0",
	"chr1	4	G	A	2	1	1	0",
	"chr1	4	G	T	2	0	1	0",
	"chr1	5	A	C	2	0	0	1",
	"chr1	6	G	A	2	0	0	1",
	"chr1	10	G	A	2	1	0	1",
	"chr1	11	G	A	0	1	0	0",
	"chr1	13	G	A	0	1	0	0",
	"chr1	14	G	A	0	1	0	0",
	"chr1	35	G	A	0	1	0	1",
	"chr1	36	G	A	0	1	0	0",
	"chr1	37	G	A	0	1	0	0",
	"chr1	38	G	A	0	1	0	0",
	"chr1	59	G	A	0	1	0	0",
	"chr1	60	G	A	0	1	0	1");
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
  $value1=~s/\r//g;
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
