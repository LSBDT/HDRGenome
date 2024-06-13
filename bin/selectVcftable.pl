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
use vars qw($opt_c $opt_C $opt_f $opt_g $opt_h $opt_o $opt_r $opt_s $opt_S $opt_u $opt_v);
getopts('cCfgho:rsSuv');
############################## HELP ##############################
sub help{
	print STDERR "\n";
	print STDERR "Command: $program_name [option] FILE BED\n";
	print STDERR "Arguments:\n";
	print STDERR "  FILE  Table from vcftable.pl or position/region file from findrun.pl or genomecov.pl\n";
	print STDERR "   BED  Bed file to include/exclude\n";
	print STDERR "Options:\n";
	print STDERR "	-c Reporting the number of overlapping features\n";
	print STDERR "	-C Reporting the number of overlapping features for each database file\n";
	print STDERR "	-f Requiring a minimal overlap fraction\n";
	print STDERR "	-g Define an alternate chromosome sort order via a genome file.\n";
	print STDERR "	-o Options for 'wa,wb,loj,wo,wao,abam,ubam,bed,split,sorted,header'\n";
	print STDERR "	-r, and -f Requiring reciprocal minimal overlap fraction\n";
	print STDERR "	-s Enforcing same strandedness\n";
	print STDERR "	-S Enforcing opposite \"strandedness\"\n";
	print STDERR "	-u (unique) Reporting the mere presence of any overlapping features\n";
	print STDERR "	-v Reporting the absence of any overlapping features\n";
	print STDERR "\n";
	print STDERR "Author: akira.hasegawa\@riken.jp\n";
	print STDERR "Update: $program_version\n";
	print STDERR "\n";
}
############################## MAIN ##############################
if($ARGV[0]eq"sortsubs"){sortSubs();exit();}
elsif($ARGV[0]eq"test"){test();exit();}
elsif(defined($opt_h)||scalar(@ARGV)<2){help();exit();}
my $vcfFile=shift(@ARGV);
my $bedFileB=shift(@ARGV);
my ($label,$bedFileA)=convertToBed($vcfFile);
my $intersectBed=intersectBed($bedFileA,$bedFileB);
my $intersectFile=convertFromBed($label,$intersectBed);
system("cat $intersectFile");
############################## convertFromBed ##############################
sub convertFromBed{
	my $label=shift();
	my $bedFile=shift();
	my $format=getFormatFromFilename($bedFile);
	my $reader=openFile($bedFile);
	my ($writer,$tmpfile)=tempfile(TEMPLATE=>"${format}_XXXXXX",DIR=>"/tmp",SUFFIX=>".txt",UNLINK=>1);
	print $writer "$label\n";
	if($format eq "position"){
		while(<$reader>){
			chomp;s/\r//g;
			my ($chr,$start,$end)=split(/\t/);
			print $writer "$chr\t$start\n";
		}
	}elsif($format eq "region"){
		while(<$reader>){
			chomp;s/\r//g;
			my ($chr,$start,$end)=split(/\t/);
			print $writer "$chr\t$start\t$end\n";
		}
	}elsif($format eq "vcftable"){
		while(<$reader>){
			chomp;s/\r//g;
			my ($chr,$start,$end,$ref,$alt,@counts)=split(/\t/);
			print $writer "$chr\t$start\t$ref\t$alt\t".join("\t",@counts)."\n";
		}
	}
	close($writer);
	close($reader);
	return $tmpfile;
}
############################## convertToBed ##############################
sub convertToBed{
	my $file=shift();
	my $reader=openFile($file);
	#CHROM	POS (position)
	#Chr	Start	End (region)
	#chromosome	position reference alternative	NA18939_v2 .. (vcftable)
	my $label=<$reader>;
	chomp($label);
    my $format;
    if($label eq "#CHROM\tPOS"){$format="position";}
    elsif($label eq "#Chr\tStart\tEnd"){$format="region";}
    elsif($label=~/^#chromosome\tposition\t/){$format="vcftable";}
    elsif($label!~/^#/){print STDERR "Error: $file is not a valid file\n";exit(1);}
	my @tokens=split(/\t/,$label);
	if(scalar(@tokens)<=1){
		print STDERR "Error: $file is not a valid file\n";
		exit(1);
	}
		if(!defined($format)){
			if(scalar(@tokens)==2){$format="position";}
			elsif(scalar(@tokens)==3){$format="region";}
			elsif(scalar(@tokens>=4)){$format="vcftable";}
		}
	my ($writer,$tmpfile)=tempfile(TEMPLATE=>"${format}_XXXXXX",DIR=>"/tmp",SUFFIX=>".bed",UNLINK=>1);
	if($format eq "position"){
		while(<$reader>){
			chomp;s/\r//g;
			my ($chr,$position)=split(/\t/);
			print $writer "$chr\t$position\t".($position+1)."\n";
		}
	}elsif($format eq "region"){
		while(<$reader>){
			chomp;s/\r//g;
			my ($chr,$start,$end)=split(/\t/);
			print $writer "$chr\t$start\t$end\n";
		}
	}elsif($format eq "vcftable"){
		while(<$reader>){
			chomp;s/\r//g;
			my ($chr,$position,$ref,$alt,@counts)=split(/\t/);
			print $writer "$chr\t$position\t".($position+1)."\t$ref\t$alt\t".join("\t",@counts)."\n";
		}
	}
	close($reader);
	return ($label,$tmpfile);
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
############################## getFormatFromFilename ##############################
sub getFormatFromFilename{
	my $file=shift();
	my $basename=basename($file,".bed");
	my @tokens=split(/_/,$basename);
	my $format=shift(@tokens);
	return $format;
}
############################## intersectBed ##############################
#https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
sub intersectBed{
	my $bedFileA=shift();
	my $bedFileB=shift();
	my $format=getFormatFromFilename($bedFileA);
	my $command="bedtools intersect -wa";
	my $booleanV=0;
	if(defined($opt_u)){$command.=" -u";}
	if(defined($opt_c)){$command.=" -c";}
	if(defined($opt_C)){$command.=" -C";}
	if(defined($opt_v)){$command.=" -v";$booleanV=1;}
	if(defined($opt_f)){$command.=" -f";}
	if(defined($opt_r)){$command.=" -r";}
	if(defined($opt_s)){$command.=" -s";}
	if(defined($opt_S)){$command.=" -S";}
	if(defined($opt_g)){$command.=" -g";}
	if(defined($opt_o)){
		my @options=split(/,/,$opt_o);
		foreach my $option(@options){
			if($option eq "v"){$booleanV=1;}
			$command.=" -$option";
		}
	}
	if($booleanV==0){$command.=" -u";}
	$command.=" -a $bedFileA -b $bedFileB";
	my ($writer,$tmpfile)=tempfile(TEMPLATE=>"${format}_XXXXXX",DIR=>"/tmp",SUFFIX=>".bed",UNLINK=>1);
	close($writer);
	$command.=" > $tmpfile";
	system($command);
	return $tmpfile;
}
############################## mkdirs ##############################
sub mkdirs{
	my @directories=@_;
	foreach my $directory(@directories){
		if($directory=~/^(.+)\@(.+)\:(.+)/){system("ssh $1\@$2 'mkdir -p $3'");}
		else{system("mkdir -p $directory");}
	}
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
		if(ref($_)eq"ARRAY"){
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
############################## test ##############################
sub test{
	#testing basic input and output functionality
	createFile("test/select.bed","chr1\t100\t200");
    #position
	createFile("test/position.txt","#CHROM\tPOS","chr1\t1","chr1\t99","chr1\t100","chr1\t101","chr1\t199","chr1\t200","chr1\t201");
	testCommand("perl $program_directory/selectVcftable.pl test/position.txt test/select.bed","#CHROM\tPOS","chr1\t100","chr1\t101","chr1\t199");
	testCommand("perl $program_directory/selectVcftable.pl -o 'v' test/position.txt test/select.bed","#CHROM\tPOS","chr1\t1","chr1\t99","chr1\t200","chr1\t201");
    #region
	createFile("test/region.txt","#Chr\tStart\tEnd","chr1\t1\t2","chr1\t99\t100","chr1\t100\t101","chr1\t101\t102","chr1\t199\t200","chr1\t200\t201","chr1\t201\t202");
    testCommand("perl $program_directory/selectVcftable.pl test/region.txt test/select.bed","#Chr\tStart\tEnd","chr1\t100\t101","chr1\t101\t102","chr1\t199\t200");
	testCommand("perl $program_directory/selectVcftable.pl -o 'v' test/region.txt test/select.bed","#Chr\tStart\tEnd","chr1\t1\t2","chr1\t99\t100","chr1\t200\t201","chr1\t201\t202");
    #vcftable
	createFile("test/vcftable.txt","#chromosome\tposition\tA\tB\tC","chr1\t1\t0\t0\t1","chr1\t99\t0\t1\t1","chr1\t100\t1\t1\t1","chr1\t101\t1\t1\t2","chr1\t199\t1\t2\t2","chr1\t200\t2\t2\t2","chr1\t201\t2\t2\t3");
    testCommand("perl $program_directory/selectVcftable.pl test/vcftable.txt test/select.bed","#chromosome\tposition\tA\tB\tC","chr1\t100\t1\t1\t1","chr1\t101\t1\t1\t2","chr1\t199\t1\t2\t2");
    testCommand("perl $program_directory/selectVcftable.pl -o 'v' test/vcftable.txt test/select.bed","#chromosome\tposition\tA\tB\tC","chr1\t1\t0\t0\t1","chr1\t99\t0\t1\t1","chr1\t200\t2\t2\t2","chr1\t201\t2\t2\t3");
}
############################## testCommand ##############################
sub testCommand{
	my @values=@_;
	my $command=shift(@values);
	my $value2=join("\n",@values);
	my ($writer,$file)=tempfile(DIR=>"/tmp",SUFFIX=>".txt",UNLINK=>1);
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
	if(compareValues($value1,$value2)==0){return;}
	print STDERR ">$command\n";
	if(ref($value2)eq"ARRAY"||ref($value2)eq"ARRAY"){printTable($value1);printTable($value2);}
	else{print STDERR "'$value1' != '$value2'\n";}
}
