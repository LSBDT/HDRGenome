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
my $program_version="2023/08/30";
############################## OPTIONS ##############################
use vars qw($opt_b $opt_h);
getopts('b:h');
############################## HELP ##############################
sub help{
  print STDERR "\n";
  print STDERR "Command: $program_name [option] STATS TELOM ANNOT BED GENOME\n";
  print STDERR "Arguments:\n";
  print STDERR " STATS  Result from statdeRSl or maxstat\n";
  print STDERR " TELOM  Telomere BED file\n";
  print STDERR " ANNOT  Annotation BED file\n";
  print STDERR "GENOME  GENOME length file\n";
  print STDERR "\n";
  print STDERR "Options:\n";
  print STDERR "     -b  bed file used to limit VCF table regions\n";
  print STDERR "Author: akira.hasegawa\@riken.jp\n";
  print STDERR "Update: $program_version\n";
  print STDERR "\n";
}
############################## MAIN ##############################
if($ARGV[0]eq"sortsubs"){sortSubs();exit();}
elsif(defined($opt_h)||scalar(@ARGV)<4){help();exit();}
my $statsFile=shift(@ARGV);
my $telomereFile=shift(@ARGV);
my $annotationFile=shift(@ARGV);
my $genomeFile=shift(@ARGV);
my $jsFile=shift(@ARGV);
my $filename=basename($statsFile,".out");
my $assembly=basename($genomeFile,".chrom.sizes");
my ($program,$bedFile,$results)=readStats($statsFile,readBedRegions($opt_b));
intersectTelomere($results,$bedFile,$telomereFile);
findClosest($results,$bedFile,$annotationFile,$genomeFile);
printResult($results,$program,$assembly,$jsFile,$filename);
############################## absolutePath ##############################
sub absolutePath{
	my $path=shift();
	my $directory=dirname($path);
	my $filename=basename($path);
	return Cwd::abs_path($directory)."/$filename";
}
############################## checkBedRegions ##############################
sub checkBedRegions{
	my $regions=shift();
	if(!defined($regions)){return 0;}#continue
	my $chr=shift();
	if(!exists($regions->{$chr})){return 1;}#skip
	my $start=shift();
	my $end=shift();
	foreach my $array(@{$regions->{$chr}}){
		my ($s,$e)=@{$array};
		if($s<=$start&&$start<$e){return 0;}#intersect
		if(!defined($end)){next;}
		if($s<=$end&&$end<$e){return 0;}#intersect
	}
	return 1;#skip
}
############################## findClosest ##############################
sub findClosest{
	my $results=shift();
	my $bedA=shift();
	my $bedB=shift();
	my $genomeFile=shift();
	my ($fh,$tmpfile)=tempfile(DIR=>"/tmp",TEMPLATE=>"XXXXXX",SUFFIX=>".txt");
	close($fh);
	system("bedtools closest -g $genomeFile -D ref -a $bedA -b $bedB>$tmpfile");
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
############################## getMetaDataFromFilename ##############################
sub getMetaDataFromFilename{
	my $filename=shift();
	my $assembly=shift();
	my @metadata=();
	my ($name,$other)=split(/\./,$filename);
	my ($target,@metas)=split(/_/,$other);
	push(@metadata,"target:$target");
	foreach my $meta(@metas){
		if($meta=~/pick(\d+)/){push(@metadata,"pick hits >: $1");next;}
		if($meta=~/skip(\d+)/){push(@metadata,"skip/mismatch: $1");next;}
		if($meta=~/min(\d+)/){push(@metadata,"minimum region size: $1");next;}
		if($meta=~/max(\d+)/){push(@metadata,"maximum region size: $1");next;}
		if($meta=~/top(\d+)/){push(@metadata,"top region selected: $1");next;}
		if($meta=~/full/){push(@metadata,"output in full mode");next;}
		if($meta=~/noindel/){push(@metadata,"skipping insertion deletion");next;}
		if($meta=~/dd/){push(@metadata,"HDR mode: DD");next;}
		if($meta=~/ar/){push(@metadata,"HDR mode: AR");next;}
		if($meta=~/ad/){push(@metadata,"HDR mode: AD");next;}
	}
	return ($name,@metadata);
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
	unlink($tmpfile);
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
	my $filename=shift();
	if($assembly=~/^([^\.]+)\./){$assembly=$1;}
	my ($name,@metadata)=getMetaDataFromFilename($filename);
	print "<html>\n";
	print "<head>\n";
	print "<title>$name</title>\n";
	print "</head>\n";
	#https://neil.fraser.name/software/tablesort/
	if(defined($jsFile)){
		print "<script>\n";
		my $reader=openFile($jsFile);
		while(<$reader>){print;}
		print "</script>\n";
	}
	print "<style>\n";
	print "table,td,th{border: 2px #000000 solid;}\n";
	print ".num{text-align:right;}\n";
	print "</style>\n";
	print "<body>\n";
	print "<h1>$name</h1>\n";
	print "<ul>\n";
	print "<li>program: $program</li>\n";
	print "<li>assembly: $assembly</li>\n";
	foreach my $meta(@metadata){print "<li>$meta</li>\n";}
	print "</ul>\n";
	print "<table>\n";
	if($program eq "statdel"){
		print "<thead>\n<tr><th class=\"num\">Group_1</th><th class=\"num\">Group_2</th><th class=\"num\">n1</th><th class=\"num\">n2</th><th class=\"num\">Teststat</th><th class=\"num\">n</th><th class=\"num\">p</th><th class=\"case\">genome</th><th class=\"case\">mH</th><th class=\"num\">rank</th><th class=\"num\">Line</th><th class=\"case\">annotation</th></tr>\n</thead>\n";
		print "<tbody>\n";
		#Group_1	Group_2	n1	n2	Teststat	n	p	chr	start_bp	end_bp	mH	rank	Line
		foreach my $result(@{$results}){
			my $group_1=$result->[0];
			my $group_2=$result->[1];
			my $n1=$result->[2];
			my $n2=$result->[3];
			my $teststat=$result->[4];
			my $n=$result->[5];
			my $p=$result->[6];
			my $chr=$result->[7];
			my $start_bp=$result->[8];
			my $end_bp=$result->[9];
			my $mH=$result->[10];
			my $rank=$result->[11];
			my $line=$result->[12];
			my $annotation=$result->[13];
			my $genome=" <a href=\"https://genome.ucsc.edu/cgi-bin/hgTracks?db=$assembly&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=$chr%3A${start_bp}%2D${end_bp}\">$chr:$start_bp-$end_bp</a>";
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
	}else{
		print "<thead>\n<tr><th class=\"num\">Var</th><th class=\"num\">Group_1</th><th class=\"num\">Group_2</th><th class=\"num\">n1</th><th class=\"num\">n2</th><th class=\"num\">Teststat</th><th class=\"num\">R</th><th class=\"num\">p</th><th class=\"case\">genome</th><th class=\"case\">median</th><th class=\"num\">HDR</th><th class=\"num\">rank</th><th class=\"case\">annotation</th></tr>\n</thead>\n";
		print "<tbody>\n";
		foreach my $result(@{$results}){
			# Var Group_1 Group_2 n1 n2 Teststat R      p chr   position median HDR   rank
			# 130	0.8737	0.1527	3	3	19.9671	30	0.2500	chr1	19499560	ok	0.7500	1
			my $var=$result->[0];
			my $group_1=$result->[1];
			my $group_2=$result->[2];
			my $n1=$result->[3];
			my $n2=$result->[4];
			my $teststat=$result->[5];
			my $r=$result->[6];
			my $p=$result->[7];
			my $chr=$result->[8];
			my $position=$result->[9];
			my $median=$result->[10];
			my $hdr=$result->[11];
			my $rank=$result->[12];
			my $annotation=$result->[13];
			my $start_bp=$position-$r;
			my $end_bp=$position+$r;
			my $genome=" <a href=\"https://genome.ucsc.edu/cgi-bin/hgTracks?db=$assembly&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=$chr%3A${start_bp}%2D${end_bp}\">$chr:$start_bp-$end_bp</a>";
			print "<tr><td>$var</td><td>$group_1</td><td>$group_2</td><td>$n1</td><td>$n2</td><td>$teststat</td><td>$r</td><td>$p</td><td>$genome</td><td>$median</td><td>$hdr</td><td>$rank</td>";
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
	}
	print "</tbody>\n";
	print "</table>\n";
	print "</body>\n";
	print "</html>\n";
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
############################## readBedRegions ##############################
sub readBedRegions{
	my $bedFile=shift();
	if(!defined($bedFile)){return;}
	my $regions={};
	my $reader=openFile($bedFile);
	my $hit=0;
	while(<$reader>){
		chomp;
		my ($chr,$start,$end,@tokens)=split(/\t/);
		if(!exists($regions->{$chr})){$regions->{$chr}=[];}
		push(@{$regions->{$chr}},[$start,$end]);
		$hit++;
	}
	if($hit>0){return $regions;}
}
############################## readStats ##############################
sub readStats{
	my $file=shift();
	my $regions=shift();
	my $reader=openFile($file);
	my @columns;
	my $program;
	while(<$reader>){
		# Var Group_1 Group_2    n1    n2 Teststat      R      p chr   position median HDR   rank
		#Group_1	Group_2	n1	n2	Teststat	n	p	chr	start_bp	end_bp	mH	rank	Line
		if(/^\s+(.+)$/){$_=$1;}
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
			if(checkBedRegions($regions,$chr,$start,$end)){next;}
			my $name="id$index";
			print $fh "$chr\t$start\t$end\t$name\t0\t.\n";
			$index++;
		}
	}elsif($program eq "maxstatRS"){
		# 130	0.8737	0.1527	3	3	19.9671	30	0.2500	chr1	19499560	ok	0.7500	1
		while(<$reader>){
			if(/^\s+$/){last;}
			if(/Note:/){last;}
			chomp;
			my @tokens=split(/\t/);
			push(@results,\@tokens);
			my $margin=$tokens[6];
			if($margin==-9){$margin=0;}
			my $chr=$tokens[8];#chr1
			my $start=$tokens[9]-$margin;#19499560-30
			my $end=$tokens[9]+$margin;#19499560+30
			if(checkBedRegions($regions,$chr,$start,$end)){next;}
			my $name="id$index";
			print $fh "$chr\t$start\t$end\t$name\t0\t.\n";
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
