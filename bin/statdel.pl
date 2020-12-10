#!/usr/bin/perl
use strict 'vars';
use Cwd;
use File::Basename;
use File::Temp;
use FileHandle;
use Getopt::Std;
use File::Temp qw/tempfile tempdir/;
############################## HEADER ##############################
my ($program_name,$prgdir,$program_suffix)=fileparse($0);
$prgdir=Cwd::abs_path($prgdir);
my $program_path="$prgdir/$program_name";
############################## OPTIONS ##############################
use vars qw($opt_h $opt_o);
getopts('ho:');
############################## HELP ##############################
sub help{
	print "\n";
	print "############################## HELP ##############################\n";
	print "\n";
	print "Program: statdel wrapper program.\n";
	print "Version: 2020/12/10\n";
	print "Author: Akira Hasegawa (akira.hasegawa\@riken.jp)\n";
	print "\n";
	print "Usage: perl statdel.pl -o OUTDIR INPUT\n";
	print "\n";
	print "  INPUT  output directory/file from hdr.pl\n";
	print " OUTDIR  Directory to store statdel results\n";
	print "\n";
}
############################## MAIN ##############################
if(defined($opt_h)||scalar(@ARGV)<1){help();exit(1);}
my @inputs=@ARGV;
my $outdir=(defined($opt_o))?$opt_o:"out";
mkdir($outdir);
@inputs=listFiles(".txt",@inputs);
foreach my $input(@inputs){
	my $tmpInput=prepareInput($input);
	my $tmpOutput=prepareOutput();
	my $output="$outdir/".basename($input,".txt").".out";
	my $param=prepareParam($tmpInput,$tmpOutput);
	system("$prgdir/statdel $param");
	handleOutput($tmpOutput,$output);
}
############################## handleOutput ##############################
sub handleOutput{
	my $tmpfile=shift();
	my $output=shift();
	my ($writer,$tmpfile2)=tempfile();
	open(IN,$tmpfile);
	while(<IN>){
		chomp;
		print $writer "$_\n";
		if(/Mean HDR/){
			my $label=<IN>;
			$label=chomp2($label);
			my @token=split(/\s+/,$label);
			print $writer join("\t",@token)."\n";
			last;
		}
	}
	while(<IN>){
		my $line=chomp2($_);
		my @token=split(/\s+/,$line);
		if($token[7]==23){$token[7]="chrX";}
		elsif($token[7]==24){$token[7]="chrY";}
		elsif($token[7]==25){$token[7]="chrM";}
		else{$token[7]="chr".$token[7];}
		$token[8]=~s/,//g;
		$token[9]=~s/,//g;
		print $writer join("\t",@token)."\n";
	}
	close(IN);
	close($writer);
	system("mv $tmpfile2 $output");
}
############################## chomp2 ##############################
sub chomp2{
	my $line=shift();
	chomp($line);
	if($line=~/^\s+(.+)$/){$line=$1;}
	return $line;
}
############################## prepareOutput ##############################
sub prepareOutput{
	my ($writer,$tmp)=tempfile(UNLINK=>1);
	close($writer);
	unlink($tmp);
	return $tmp;
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
	my $type=($labels[1]eq"Position")?"ar":"dd";
	my @lines=();
	my $replaces={};
	while(<IN>){
		chomp;
		my @token=split(/\t/);
		if($token[0]=~/^chr\d+$/){}
		elsif($token[0]=~/^chrX$/){$token[0]="chr23";$replaces->{"chrX"}="chr23";}
		elsif($token[0]=~/^chrY$/){$token[0]="chr24";$replaces->{"chrY"}="chr24";}
		elsif($token[0]=~/^chrM$/){$token[0]="chr25";$replaces->{"chrM"}="chr25";}
		else{next;}
		push(@lines,\@token);
	}
	close(IN);
	if($type ne "dd"){
		my $hash={};
		foreach my $token(@lines){$hash->{$token->[2]}++;}
		my $diffsize=scalar(keys(%{$hash}));
		my $totalsize=scalar(@lines)/$diffsize/2;
		print $writer "$totalsize $ctrlsize $diffsize";
		print $writer "\n";
	}else{
		my $totalsize=scalar(@lines)/2;
		print $writer "$totalsize $ctrlsize 1";
		print $writer "\n";
	}
	print $writer "$label\n";
	foreach my $token(@lines){print $writer join("\t",@{$token})."\n";}
	close($writer);
	return $tmp;
}
############################## prepareParam ##############################
sub prepareParam{
	my $input=shift();
	my $output=shift();
	my ($writer,$tmp)=tempfile(UNLINK=>1);
	print $writer "Statdel: Pt369 with known pathogenic deletions\n";
	print $writer "-9 17 14005439 15217437  17p12  code for missing, chrom, start and end bp position of dis.del, 0 0 0 if unknown\n";
	print $writer "-12 1 1      code for test statistic, exponent, include obs\n";
	print $writer "1 3 0.5 0.8  no. of \"n\" values, min bp HDR, min median(HDR), min. overlap\n";
	print $writer "$input\n";
	print $writer "$output\n";
	close($writer);
	return $tmp;
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