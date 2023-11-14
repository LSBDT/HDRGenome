my $stats={};
my @labels=("wild","het","hom");
print join("\t",@labels)."\ttotal\n";
while(<STDIN>){
	chomp;s/\r//g;
	my ($chr,$pos,@flags)=split(/\t/);
	my $hash={};
	foreach my $flag(@flags){
		if($flag&1){$flag="het";}
		elsif($flag&2){$flag="hom";}
		else{$flag="wild";}
		$hash->{$flag}++;
	}
	my $string;
	foreach my $key(@labels){
		my $val=exists($hash->{$key})?$hash->{$key}:0;
		if(defined($string)){$string.="\t";}
		$string.="$val";
	}
	$stats->{$string}++;
}
foreach my $key(sort{$stats->{$b}<=>$stats->{$a}}keys(%{$stats})){
	my $val=$stats->{$key};
	print "$key\t$val\n";
}