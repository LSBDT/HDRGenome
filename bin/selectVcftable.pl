#!/usr/bin/perl
my $vcfFile=shift(@ARGV);
my $bedFile=shift(@ARGV);
my $regions=readBedFile($bedFile);
selectRegions($vcfFile,$regions);
############################## selectRegions ##############################
sub selectRegions{
    my $vcfFile=shift();
    my $regions=shift();
    my $reader=openFile($vcfFile);
    my $preChr;
    my $region;
    my $label=<$reader>;
    print $label;
    while(<$reader>){
        chomp;
		s/\r//g;
        my ($chr,$position,@counts)=split(/\t/);
        if($preChr ne $chr){$region=exists($regions->{$chr})?$regions->{$chr}:undef;$preChr=$chr;}
        if(!defined($region)){next;}
        my ($start,$end)=@{$region->[0]};
        while($end<$position){
            shift(@{$region});
            if(scalar(@{$region})<1){$region=undef;delete($regions->{$chr});last;}
            ($start,$end)=@{$region->[0]};
        }
        if(!defined($region)){next;}
        if($position<$start){next;}
        print "$chr\t$position\t".join("\t",@counts)."\n";
    }
    close($reader);
}
############################## readBedFile ##############################
sub readBedFile{
    my $file=shift();
    my $reader=openFile($file);
    my $regions={};
    while(<$reader>){
        chomp;
		s/\r//g;
        my ($chr,$start,$end,@tokens)=split(/\t/);
        if(!defined($regions->{$chr})){$regions->{$chr}=[];}
        push(@{$regions->{$chr}},[$start,$end]);
    }
    close($reader);
    return $regions;
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