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
use vars qw($opt_h);
getopts('h');
############################## HELP ##############################
sub help{
  print STDERR "\n";
  print STDERR "Command: $program_name [option] STATS BED GENOME\n";
  print STDERR "Arguments:\n";
  print STDERR " STATS  Result from statdeRSl or maxstat\n";
  print STDERR "   BED  BED file with gene locations\n";
  print STDERR "GENOME  GENOME length file with gene locations\n";
  print STDERR "\n";
  print STDERR "Author: akira.hasegawa\@riken.jp\n";
  print STDERR "Update: $program_version\n";
  print STDERR "\n";
}
############################## MAIN ##############################
if($ARGV[0]eq"sortsubs"){sortSubs();exit();}
elsif(defined($opt_h)||scalar(@ARGV)<2){help();exit();}
my $statsFile=shift(@ARGV);
my $telomereFile=shift(@ARGV);
my $annotationFile=shift(@ARGV);
my $genomeFile=shift(@ARGV);
my $jsFile=shift(@ARGV);
my $assembly=basename($genomeFile,".genome");
my ($program,$bedFile,$results)=readStats($statsFile);
intersectTelomere($results,$bedFile,$telomereFile);
findClosest($results,$bedFile,$annotationFile,$genomeFile);
printResult($results,$program,$assembly,$jsFile);
############################## absolutePath ##############################
sub absolutePath{
	my $path=shift();
	my $directory=dirname($path);
	my $filename=basename($path);
	return Cwd::abs_path($directory)."/$filename";
}
############################## findClosest ##############################
sub findClosest{
	my $results=shift();
	my $bedA=shift();
	my $bedB=shift();
	my $genomeFile=shift();
	my ($fh,$tmpfile)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".txt");
	close($fh);
	system("bedtools closest -g $genomeFile -d -a $bedA -b $bedB>$tmpfile");
	#system("cat $tmpfile");
	my $reader=openFile($tmpfile);
	while(<$reader>){
		chomp;
		my @token=split(/\t/);
		my $id=$token[3];
		my $gene=$token[9];
		my $distance=$token[12];
		if($id=~/^id(\d+)$/){
		$id=$1;
		push(@{$results->[$id]->[13]},[$gene,$distance]);
		}
	}
	close($reader);
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
############################## intersectTelomere ##############################
sub intersectTelomere{
	my $results=shift();
	my $bedA=shift();
	my $bedB=shift();
	if(scalar(@{$results})==0){return;}
	my ($fh,$tmpfile)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".txt");
	close($fh);
	system("bedtools intersect -a $bedA -b $bedB > $tmpfile");
	#system("cat $tmpfile");
	my $reader=openFile($tmpfile);
	while(<$reader>){
		chomp;
		my @token=split(/\t/);
		my $id=$token[3];
		if($id=~/^id(\d+)$/){
		$id=$1;
		push(@{$results->[$id]->[13]},["telomere",0]);
		}
	}
	close($reader);
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
############################## printResult ##############################
sub printResult{
	my $results=shift();
	my $program=shift();
	my $assembly=shift();
	my $jsFile=shift();
	#https://neil.fraser.name/software/tablesort/
	if(defined($jsFile)){
		print "<script>\n";
		my $reader=openFile($jsFile);
		while(<$reader>){print;}
		print "</script>\n";
	}
	print "<style>\n";
	print "table,td,th{border: 2px #000000 solid;}\n";
	print "</style>\n";
	print "<table>\n";
	print "<thead>\n<tr><th class=\"num\">Group_1</th><th class=\"num\">Group_2</th><th class=\"num\">n1</th><th class=\"num\">n2</th><th class=\"num\">Teststat</th><th class=\"num\">n</th><th class=\"num\">p</th><th class=\"case\">genome</th><th class=\"case\">mH</th><th class=\"num\">rank</th><th class=\"num\">Line</th><th class=\"case\">annotation</th></tr>\n</thead>\n";
	print "<tbody>\n";
	foreach my $result(@{$results}){
		my $group_1=$result->[0];
		my $group_2=$result->[1];
		my $n1=$result->[2];
		my $n2=$result->[3];
		my $teststat=$result->[4];
		my $n=$result->[5];
		my $p=$result->[6];
		my $chr=($program eq "statdel")?$result->[7]:$result->[8];
		my $start_bp=($program eq "statdel")?$result->[8]:$result->[9]-100;
		my $end_bp=($program eq "statdel")?$result->[9]:$result->[9]+100;
		my $genome=" <a href=\"https://genome.ucsc.edu/cgi-bin/hgTracks?db=$assembly&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=$chr%3A${start_bp}%2D${end_bp}\">$chr:$start_bp-$end_bp</a>";
		my $mH=$result->[10];
		my $rank=$result->[11];
		my $line=$result->[12];
		my $annotation=$result->[13];
		print "<tr><td>$group_1</td><td>$group_2</td><td>$n1</td><td>$n2</td><td>$teststat</td><td>$n</td><td>$p</td><td>$genome</td><td>$mH</td><td>$rank</td><td>$line</td>";
		my $line=shift();
		foreach my $ann(@{$annotation}){
			my ($gene,$dist)=@{$ann};
			if(defined($line)){$line.=".";}
			if($gene ne "telomere"){$gene=" <a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=$gene\">$gene</a>"}
			$line.="$gene($dist)";
		}
		print "<td>$line</td>";
		print "</tr>\n";
	}
	print "<tbody>\n";
	print "<table>\n";
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
############################## readStats ##############################
sub readStats{
	my $file=shift();
	my $reader=openFile($file);
	my @columns;
	my $program;
	while(<$reader>){
		#Var Group_1 Group_2    n1    n2 Teststat      R      p chr   position median HDR   rank
		#Group_1	Group_2	n1	n2	Teststat	n	p	chr	start_bp	end_bp	mH	rank	Line
		if(/Var Group_1/){chomp;@columns=split(/\s+/);$program="maxstatRS";last;}
		if(/Group_1/){chomp;@columns=split(/\t/);$program="statdel";last;}
	}
	my ($fh,$tmpfile)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".bed");
	my $index=0;
	my @results=();
	if($program eq "statdel"){
		while(<$reader>){
			if(/^\s+$/){last;}
			if(/Note:/){last;}
			chomp;
			my @tokens=split(/\t/);
			push(@results,\@tokens);
			my $chr=$tokens[7];
			my $start=$tokens[8];
			my $end=$tokens[9];
			my $name="id$index";
			print $fh "$chr\t$start\t$end\t$name\t0\t.\n";
			$index++;
		}
	}elsif($program eq "maxstatRS"){
		while(<$reader>){
			if(/^\s+$/){last;}
			if(/Note:/){last;}
			chomp;
			my @tokens=split(/\t/);
			push(@results,\@tokens);
			my $chr=$tokens[8];
			my $start=$tokens[9];
			my $name="id$index";
			print $fh "$chr\t$start\t$start\t$name\t0\t.\n";
			$index++;
		}
	}
	close($fh);
	my ($fh2,$tmpfile2)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".bed");
	close($fh2);
	system("sort -k1,1 -k2,2n -k3,3n $tmpfile > $tmpfile2");
	return ($program,$tmpfile2,\@results);
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
