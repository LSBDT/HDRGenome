#!/bin/sh
if [ $# -ne 1 ]; then
echo "hdr.sh [CONFIG]"
exit
fi

perl bin/dag.pl -d db -q config $1
mkdir -p annotation
mkdir -p download

perl bin/moirai2.pl clear all

perl bin/moirai2.pl -c moirai2/biotools -P "# Creating project directory..." -d db -s 1 -m 1 -i 'root->project->$project' -o 'mkdir->$project/flag->completed' command << 'EOF'
projectdir=$project
mkdir -p $projectdir
mkdir -p $projectdir/log
mkdir -p $projectdir/findrun
mkdir -p $projectdir/genomecov
mkdir -p $projectdir/hdr
mkdir -p $projectdir/stats
mkdir -p $projectdir/html
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Setup cytobands" -d db -s 1 -m 1 -i 'root->genome->$genome' -o '$genome->cytoband->$cytoband' command << 'EOF'
wget -q -P $tmpdir/ http://hgdownload.cse.ucsc.edu/goldenpath/$genome/database/cytoBand.txt.gz
gzip -cd $tmpdir/cytoBand.txt.gz | grep acen | sort -k 1,1 -k 2,2n > $cytoband
rm $tmpdir/cytoBand.txt.gz
mv $cytoband annotation/$genome.cytoband.bed
cytoband=annotation/$genome.cytoband.bed
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Setup gencode" -d db -s 1 -m 1 -i 'root->genome->$genome' -o '$genome->gencode->$gencode' -S bin/selectGenesFromGencode.pl command << 'EOF'
if [ "$genome" = "hg38" ]; then
wget -q -O $tmpdir/annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
elif [ "$genome" = "hg19" ]; then
wget -q -O $tmpdir/annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh37_mapping/gencode.v41lift37.annotation.gtf.gz
fi
gzip -cd $tmpdir/annotation.gtf.gz | selectGenesFromGencode.pl | sort -k 1,1 -k 2,2n > $gencode
rm $tmpdir/annotation.gtf.gz
mv $gencode annotation/$genome.gencode.bed
gencode=annotation/$genome.gencode.bed
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Setup chrominfo" -d db -s 1 -m 1 -i 'root->genome->$genome' -o '$genome->chrominfo->$chrominfo' command 'chrominfo=annotation/$genome.chrom.sizes' << 'EOF'
wget -q -O $chrominfo https://hgdownload.soe.ucsc.edu/goldenPath/$genome/bigZips/$genome.chrom.sizes
cat $chrominfo | grep -v '_' | sort -k 1,1 -k 2,2n> $tmpdir/$genome.sorted
mv $tmpdir/$genome.sorted $chrominfo
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Download from 1000Genomes..." -d db -s 1 -m 1 -i '$project->download#url->$url' -o '$url->download#filepath->$filepath' command << 'EOF'
wget -q -P download $url
filename=`basename $url`
filepath="download/$filename"
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Symbolic link download files" -d db -s 1 -m 1 -i '$project->download#url->$url,$url->download#filepath->$filepath,$project->project#indir->$indir' -o 'root->$project/downloadfile->$downloadfile' command << 'EOF'
filename=`basename $filepath`
abspath=$(cd $(dirname $filepath) && pwd)/$(basename $filepath)
downloadfile="$indir/$filename"
mkdir -p $indir
ln -s $abspath $downloadfile
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# List input directory..." -d db -s 1 -m 1 -i '$project->project#indir->$input' -o '$basename->$project/input->$filepath' ls
perl bin/moirai2.pl -c moirai2/biotools -P "# List position directory..." -d db -s 1 -m 1 -i '$project->project#positiondir->$input' -o '$basename->$project/position->$filepath' ls

perl bin/moirai2.pl -c moirai2/biotools -P "# Creating VCF table..." -d db -s 60 -m 1 -S bin/vcftable.pl -i 'root->project->$project,$project->project#indir->$indir,$project->vcftable#qvThreshold->$qvThreshold' -o 'root->$project/vcftable->$output' command 'output=$project/vcftable.txt' << 'EOF'
output=$tmpdir/vcftable.txt
logFile=$project/log/vcftable.txt
vcftable.pl -t $qvThreshold -o $output $indir/* 2> $logFile
echo "INSERT root->$project/vcftable#log->$logFile"
EOF
perl bin/moirai2.pl -d db -s 60 -m 1 -S bin/vcftable.pl -S bin/getCases.pl -i 'root->$project/vcftable->$vcftable' -o 'root->$project/case->$case' command 'output=$project/vcftable.txt' << 'EOF'
case=(`getCases.pl < $vcftable`)
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Narrow down VCF table regions..." -d db -s 60 -m 1 -S bin/selectVcftable.pl -i 'root->$project/vcftable->$vcftable,$project->vcftable#selectBed->$selectBed' -o 'vcfselect->$project/flag->completed' command << 'EOF'
output=$project/vcftable.select.txt
selectVcftable.pl $vcftable $selectBed > $output
echo "UPDATE root->$project/vcftable->$output"
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Convert VCF table to BED format..." -d db -s 60 -m 1 -S bin/selectVcftable.pl -i 'root->$project/vcftable->$vcftable' -o 'root->$project/vcfbed->$output' command << 'EOF'
output=$project/vcftable.bed
bin/vcftable2bed.pl < $vcftable > $output
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Calculating positions by genomecov..." -d db -s 60 -m 1 -S bin/genomecov.pl -b '$noIndel:-d' -i '$project->hdr#targetMode->DD,root->$project/vcftable->$table,root->$project/case->$case,$project->genomecov#noIndel->$noIndel,$project->genomecov#stretchMode->$stretchMode,$project->genomecov#topX->$topX,$project->genomecov#regionSize->$regionSize,$project->annotation#genome->$genome,$genome->chrominfo->$chrominfo' -o '$case->$project/genomecov->$output' command 'output=$project/genomecov/$filename' << 'EOF'
logFile=$project/log/genomecov.$case.txt
genomecov.pl $noIndel -o $tmpdir -m $stretchMode -t $topX -r $regionSize $table $case $chrominfo > $output 2> $logFile
output=`ls $tmpdir/*.txt`
filename=`basename $output`
echo "INSERT $case->$project/genomecov#log->$logFile"
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Calculating positions by findrun..." -d db -s 60 -m 1 -S bin/findrun.pl -b '$noIndel:-d' -i '$project->hdr#targetMode->DD,root->$project/vcftable->$table,root->$project/case->$case,$project->findrun#noIndel->$noIndel,$project->findrun#stretchMode->$stretchMode,$project->findrun#pickupNumber->$pickupNumber,$project->findrun#regionSize->$regionSize,$project->findrun#skipCount->$skipCount' -o '$case->$project/findrun->$output' command 'output=$project/findrun/$filename' << 'EOF'
logFile=$project/log/findrun.$case.txt
findrun.pl $noIndel -m $stretchMode -p $pickupNumber -r $regionSize -s $skipCount -o $tmpdir $table $case 2> $logFile
output=`ls $tmpdir/*.txt`
filename=`basename $output`
echo "INSERT $case->$project/findrun#log->$logFile"
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Calculating positions by findrun with min+max regions..." -d db -s 60 -m 1 -S bin/findrun.pl -b '$noIndel:-d' -i '$project->hdr#targetMode->DD,root->$project/vcftable->$table,root->$project/case->$case,$project->findrun#noIndel->$noIndel,$project->findrun#stretchMode->$stretchMode,$project->findrun#pickupNumber->$pickupNumber,$project->findrun#regionSizeMin->$regionSizeMin,$project->findrun#regionSizeMax->$regionSizeMax,$project->findrun#skipCount->$skipCount' -o '$case->$project/findrun->$output' command 'output=$project/findrun/$filename' << 'EOF'
logFile=$project/log/findrun.$case.txt
findrun.pl $noIndel -m $stretchMode -p $pickupNumber -r $regionSizeMin -R $regionSizeMax -s $skipCount -o $tmpdir $table $case 2> $logFile
output=`ls $tmpdir/*.txt`
filename=`basename $output`
echo "INSERT $case->$project/findrun#log->$logFile"
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Calculating positions by findrun with min regions..." -d db -s 60 -m 1 -S bin/findrun.pl -b '$noIndel:-d' -i '$project->hdr#targetMode->DD,root->$project/vcftable->$table,root->$project/case->$case,$project->findrun#noIndel->$noIndel,$project->findrun#stretchMode->$stretchMode,$project->findrun#pickupNumber->$pickupNumber,$project->findrun#regionSizeMin->$regionSizeMin,$project->findrun#skipCount->$skipCount' -o '$case->$project/findrun->$output' command 'output=$project/findrun/$filename' << 'EOF'
logFile=$project/log/findrun.$case.txt
findrun.pl $noIndel -m $stretchMode -r $regionSizeMin -R $regionSizeMax -s $skipCount -o $tmpdir $table $case 2> $logFile
output=`ls $tmpdir/*.txt`
filename=`basename $output`
echo "INSERT $case->$project/findrun#log->$logFile"
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Calculating positions by findrun with max regions..." -d db -s 60 -m 1 -S bin/findrun.pl -b '$noIndel:-d' -i '$project->hdr#targetMode->DD,root->$project/vcftable->$table,root->$project/case->$case,$project->findrun#noIndel->$noIndel,$project->findrun#stretchMode->$stretchMode,$project->findrun#pickupNumber->$pickupNumber,$project->findrun#regionSizeMax->$regionSizeMax,$project->findrun#skipCount->$skipCount' -o '$case->$project/findrun->$output' command 'output=$project/findrun/$filename' << 'EOF'
logFile=$project/log/findrun.$case.txt
findrun.pl $noIndel -m $stretchMode -p $pickupNumber -R $regionSizeMax -s $skipCount -o $tmpdir $table $case 2> $logFile
output=`ls $tmpdir/*.txt`
filename=`basename $output`
echo "INSERT $case->$project/findrun#log->$logFile"
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Calculating HDR with DD mode findrun..." -d db -s 60 -m 1 -S bin/hdr.pl -b '$excludeIndel:-d,$excludeLowQuality:-l' -i '$project->hdr#targetMode->DD,root->$project/vcftable->$table,root->$project/case->$case,$case->$project/findrun->$findrun,$project->hdr#excludeIndel->$excludeIndel,$project->hdr#excludeLowQuality->$excludeLowQuality' -o '$findrun->$project/hdr->$output' command 'output=$project/hdr/$filename' << 'EOF'
logFile=$project/log/hdr.$case.txt
hdr.pl $excludeIndel $excludeLowQuality -t DD -o $tmpdir $table $case $findrun 2> $logFile
output=`ls $tmpdir/$case/*.txt`
mv $output $tmpdir/.
rmdir $tmpdir/$case
filename=`basename $output`
output=$tmpdir/$filename
echo "INSERT $findrun->$project/hdr#log->$logFile"
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Calculating HDR with DD mode genomecov..." -d db -s 60 -m 1 -S bin/hdr.pl -b '$excludeIndel:-d,$excludeLowQuality:-l' -i '$project->hdr#targetMode->DD,root->$project/vcftable->$table,root->$project/case->$case,$case->$project/genomecov->$genomecov,$project->hdr#excludeIndel->$excludeIndel,$project->hdr#excludeLowQuality->$excludeLowQuality' -o '$genomecov->$project/hdr->$output' command 'output=$project/hdr/$filename' << 'EOF'
logFile=$project/log/hdr.$case.txt
hdr.pl $excludeIndel $excludeLowQuality -t DD -o $tmpdir $table $case $genomecov 2> $logFile
output=`ls $tmpdir/$case/*.txt`
mv $output $tmpdir/.
rmdir $tmpdir/$case
filename=`basename $output`
output=$tmpdir/$filename
echo "INSERT $genomecov->$project/hdr#log->$logFile"
EOF

perl bin/moirai2.pl -c moirai2/biotools -P "# Calculating HDR with AR/AD mode..." -d db -s 60 -m 1 -S bin/hdr.pl -b '$excludeIndel:-d,$excludeLowQuality:-l' -i '$project->hdr#targetMode->$targetMode,root->$project/vcftable->$table,root->$project/case->$case,$case->$project/position->$position,$project->hdr#interval->$interval,$project->hdr#startDistance->$startDistance,$project->hdr#endDistance->$endDistance,$project->hdr#excludeIndel->$excludeIndel,$project->hdr#excludeLowQuality->$excludeLowQuality' -o '$position->$project/hdr->$output' command 'output=$project/hdr/$filename' << 'EOF'
logFile=$project/log/hdr.$case.txt
hdr.pl $excludeIndel $excludeLowQuality -t $targetMode -s $startDistance -e $endDistance -i $interval -o $tmpdir $table $case $position 2> $logFile
output=`ls $tmpdir/$case/*.txt`
mv $output $tmpdir/.
rmdir $tmpdir/$case
filename=`basename $output`
output=$tmpdir/$filename
echo "INSERT $position->$project/hdr#log->$logFile"
EOF

#AR/AD/DD
perl bin/moirai2.pl -c moirai2/biotools -P "# Calculating statdel/maxStats..." -d db -s 10 -m 1 -i '$input->$project/hdr->$hdr' -o '$hdr->$project/stats->$output' -O "Results in file" command 'output=$project/stats/$filename' << 'EOF'
basename=`basename $hdr .txt`
logFile=$project/log/stats.${basename}.txt
perl bin/statdel.pl -o $tmpdir $hdr 2> $logFile
output=`ls $tmpdir/*.out`
filename=`basename $output`
echo "INSERT $hdr->$project/stats#log->$logFile"
EOF

#AR/AD/DD
perl bin/moirai2.pl -c moirai2/biotools -P "# Parsing statdel/maxStats and creating HTML..." -d db -s 10 -m 1 -S bin/annotation.pl -i '$project->annotation#genome->$genome,$genome->chrominfo->$chrominfo,$genome->cytoband->$cytoband,$genome->gencode->$gencode,$hdr->$project/stats->$stats' -o '$stats->$project/html->$html' command 'html=$project/html/$filename.html' << 'EOF'
annotation.pl $stats $cytoband $gencode $chrominfo js/tablesort-min.js > $html
filename=`basename $stats .out`
EOF
