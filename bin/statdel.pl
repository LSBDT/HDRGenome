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
use vars qw($opt_h $opt_o);
getopts('ho:');
############################## HELP ##############################
sub help{
	print "\n";
	print "############################## HELP ##############################\n";
	print "\n";
	print "Program: statdel and maxstatRS wrapper program.\n";
	print "Version: $program_version\n";
	print "Author: Akira Hasegawa (akira.hasegawa\@riken.jp)\n";
	print "\n";
	print "Usage: perl statdel.pl -o OUTDIR INPUT\n";
	print "\n";
	print "  INPUT  output directory/file from hdr.pl\n";
	print " OUTDIR  Directory to store statdel results\n";
	print "\n";
}
############################## MAIN ##############################
if($ARGV[0]eq"sortsubs"){sortSubs();exit();}
elsif(defined($opt_h)||scalar(@ARGV)<1){help();exit();}
my @inputs=@ARGV;
my @regionSizes=();#for AR/D mode,
my $outdir=(defined($opt_o))?$opt_o:"out";
mkdir($outdir);
@inputs=listFiles(".txt",@inputs);
my $result=`uname`;
my $os=($result=~/Darwin/)?"mac":"linux";
foreach my $input(@inputs){
	my ($tmpInput,$program)=prepareInput($input);
	my $tmpOutput=prepareOutput();
	my $output="$outdir/".basename($input,".txt").".out";
	my $param=($program eq "statdel")?prepareParamStatdel($tmpInput,$tmpOutput):prepareParamMaxstatRS($tmpInput,$tmpOutput);
	print STDERR "Input file: $input\n";
	print STDERR "Temporary input: $tmpInput\n";
	print STDERR "Parameter file: $param\n";
	print STDERR "Temporary output: $tmpOutput\n";
	print STDERR "Output file: $output\n";
	print STDERR "Command: $program_directory/$os/$program\n";
	system("$program_directory/$os/$program $param");
	handleOutput($tmpOutput,$output);
}
############################## absolutePath ##############################
sub absolutePath {
	my $path=shift();
	my $directory=dirname($path);
	my $filename=basename($path);
	my $path=Cwd::abs_path($directory)."/$filename";
	$path=~s/\/\.\//\//g;
	$path=~s/\/\.$//g;
	return $path
}
############################## chomp2 ##############################
sub chomp2{
	my $line=shift();
	chomp($line);
	if($line=~/^\s+(.+)$/){$line=$1;}
	return $line;
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
############################## handleOutput ##############################
sub handleOutput{
	my $tmpfile=shift();
	my $output=shift();
	my ($writer,$tmpfile2)=tempfile();
	open(IN,$tmpfile);
	my $program;
	while(<IN>){
		chomp;
		print $writer "$_\n";
		if(/Mean HDR/){
			my $line=<IN>;
			if($line=~/Var/){$program="maxstatRS";}
			else{$program="statdel"}
			print $writer $line;
			last;
		}
	}
	if($program eq"statdel"){
		while(<IN>){
			if(/^\s+$/){print $writer $_;next;}
			if(/Note:/){print $writer $_;last;}
			my $line=chomp2($_);
			my @token=split(/\s+/,$line);
			if(scalar(@token)!=13){print STDERR "ERROR not 13:".join("\t",@token)."\n";next;}
			if($token[7]==23){$token[7]="chrX";}
			elsif($token[7]==24){$token[7]="chrY";}
			elsif($token[7]==25){$token[7]="chrM";}
			else{$token[7]="chr".$token[7];}
			$token[8]=~s/,//g;
			$token[9]=~s/,//g;
			print $writer join("\t",@token)."\n";
		}
	}elsif($program eq"maxstatRS"){
		while(<IN>){
			if(/^\s+$/){print $writer $_;next;}
			if(/Note:/){print $writer $_;last;}
			my $line=chomp2($_);
			my @token=split(/\s+/,$line);
			if(scalar(@token)!=13){print STDERR "ERROR not 13:".join("\t",@token)."\n";next;}
			if($token[8]==23){$token[8]="chrX";}
			elsif($token[8]==24){$token[8]="chrY";}
			elsif($token[8]==25){$token[8]="chrM";}
			else{$token[8]="chr".$token[8];}
			$token[9]=~s/,//g;
			print $writer join("\t",@token)."\n";
		}
	}
	close(IN);
	close($writer);
	system("mv $tmpfile2 $output");
}
############################## listFiles ##############################
sub listFiles{
	my @input_directories=@_;
	my $file_suffix=shift(@input_directories);
	my @input_files=();
	foreach my $input_directory (@input_directories){
		$input_directory=absolutePath($input_directory);
		if(-f $input_directory){push(@input_files,$input_directory);next;}# It's a file, so process file
		elsif(-l $input_directory){push(@input_files,$input_directory);next;}# It's a file, so process file
		opendir(DIR,$input_directory);
		foreach my $file(readdir(DIR)){# go through input directory
			if($file eq "."){next;}
			if($file eq "..") {next;}
			if($file eq ""){next;}
			$file="$input_directory/$file";
			if(-d $file){next;}# skip directory element
			elsif($file!~/$file_suffix$/){next;}
			push(@input_files,$file);
		}
		closedir(DIR);
	}
	return sort{$a cmp $b}@input_files;
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
############################## prepareInput ##############################
sub prepareInput{
	my $input=shift();
	my ($writer,$tmp)=tempfile(UNLINK=>1);
	open(IN,$input);
	<IN>;
	my $label=<IN>;
	chomp($label);
	my @labels=split(/\t/,$label);
	my $ctrlsize=0;
	foreach my $l(@labels){if($l=~/^Case/){$ctrlsize++;}}
	my $program=($labels[1]eq"Position")?"maxstatRS":"statdel";
	my @lines=();
	my $replaces={};
	while(<IN>){
		chomp;
		my @token=split(/\t/);
		if($token[0]=~/^chr\d+$/){}#chr1 chr2 chr3 ... chr22
		elsif($token[0]=~/^chrX$/){$token[0]="chr23";$replaces->{"chrX"}="chr23";}#chrX
		elsif($token[0]=~/^chrY$/){$token[0]="chr24";$replaces->{"chrY"}="chr24";}#chrY
		elsif($token[0]=~/^chrM$/){$token[0]="chr25";$replaces->{"chrM"}="chr25";}#chrM
		else{next;}#chrX_random_XXXXX
		push(@lines,\@token);
	}
	close(IN);
	if($program eq "statdel"){
		my $totalsize=scalar(@lines)/2;
		my $diffsize=1;
		print $writer "$totalsize $ctrlsize $diffsize";
		print $writer "\n";
	}else{
		my $hash={};
		foreach my $token(@lines){$hash->{$token->[2]}++;}
		@regionSizes=sort {$a<=>$b}keys(%{$hash});
		my $diffsize=scalar(keys(%{$hash}));
		my $totalsize=scalar(@lines)/$diffsize/2;
		my $casesize=1;
		print $writer "$totalsize $casesize $ctrlsize $diffsize";
		print $writer "\n";
	}
	print $writer "$label\n";
	foreach my $token(@lines){print $writer join("\t",@{$token})."\n";}
	close($writer);
	return ($tmp,$program);
}
############################## prepareOutput ##############################
sub prepareOutput{
	my ($writer,$tmp)=tempfile(UNLINK=>1);
	close($writer);
	unlink($tmp);
	return $tmp;
}
############################## prepareParamMaxstatRS ##############################
sub prepareParamMaxstatRS{
	my $input=shift();
	my $output=shift();
	my ($writer,$tmp)=tempfile();#UNLINK=>1
	print $writer "maxstatRS: Auto generated parameter file\n";
	print $writer "-9 0 0 0\n";
	print $writer "-12 1 1\n";
	my $line="";
	my $min=$regionSizes[0];
	my $max=$regionSizes[scalar(@regionSizes)-1];
	print $writer "$min $max 3.0\n";
	print $writer "$input\n";
	print $writer "$output\n";
	foreach my $regionSize(@regionSizes){print $writer "$regionSize\n";}
	close($writer);
	return $tmp;
}
############################## prepareParamStatdel ##############################
sub prepareParamStatdel{
	my $input=shift();
	my $output=shift();
	my ($writer,$tmp)=tempfile();#UNLINK=>1
	print $writer "Statdel: Auto generated parameter file\n";
	print $writer "-9 0 0 0\n";
	print $writer "-12 1 1\n";
	print $writer "1 3 0.5 0.8\n";
	print $writer "$input\n";
	print $writer "$output\n";
	close($writer);
	return $tmp;
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
