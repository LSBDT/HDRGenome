#!/bin/sh
if [ $# -ne 1 ]; then
echo "run.sh [CONFIG]"
exit
fi

perl bin/rdf.pl -d db -q config $1
mkdir -p annotation
mkdir -p download

#reset all previous logs
perl bin/moirai2.pl clear all

echo "# Creating project directory..."
perl bin/moirai2.pl -d db -s 1 -m 1 -i 'root->project->$project' -o '$project->project#projectdir->$projectdir' command << 'EOF'
projectdir=$project
mkdir -p $projectdir
mkdir -p $projectdir/log
mkdir -p $projectdir/findrun
mkdir -p $projectdir/hdr
mkdir -p $projectdir/stats
mkdir -p $projectdir/html
EOF

echo "# Getting annotations..."
#cytobnds
perl bin/moirai2.pl -d db -s 1 -m 1 -i 'root->genome->$genome' -o '$genome->cytoband->$cytoband' command << 'EOF'
wget -q -P $tmpdir/ http://hgdownload.cse.ucsc.edu/goldenpath/$genome/database/cytoBand.txt.gz
gzip -cd $tmpdir/cytoBand.txt.gz | grep acen > $cytoband
rm $tmpdir/cytoBand.txt.gz
mv $cytoband annotation/$genome.cytoband.bed
cytoband=annotation/$genome.cytoband.bed
EOF

#gencode
perl bin/moirai2.pl -d db -s 1 -m 1 -i 'root->genome->$genome' -o '$genome->gencode->$gencode' -S bin/selectGenesFromGencode.pl command << 'EOF'
if [ "$genome" = "hg38" ]; then
wget -q -O $tmpdir/annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
elif [ "$genome" = "hg19" ]; then
wget -q -O $tmpdir/annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh37_mapping/gencode.v41lift37.annotation.gtf.gz
fi
gzip -cd $tmpdir/annotation.gtf.gz | selectGenesFromGencode.pl > $gencode
rm $tmpdir/annotation.gtf.gz
mv $gencode annotation/$genome.gencode.bed
gencode=annotation/$genome.gencode.bed
EOF

#chrominfo
perl bin/moirai2.pl -d db -s 1 -m 1 -i 'root->genome->$genome' -o '$genome->chrominfo->$chrominfo' -S bin/selectGenesFromGencode.pl command 'chrominfo=annotation/$genome.genome' << 'EOF'
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from $genome.chromInfo" > $chrominfo
EOF

echo "# Download from 1000Genomes..."
perl bin/moirai2.pl -d db -s 1 -m 1 -i '$project->download#url->$url' -o '$url->download#filepath->$filepath' command << 'EOF'
wget -q -P download $url
filename=`basename $url`
filepath="download/$filename"
EOF

perl bin/moirai2.pl -d db -s 1 -m 1 -i '$project->download#url->$url,$url->download#filepath->$filepath,$project->project#indir->$indir' -o 'root->$project/downloadfile->$downloadfile' command << 'EOF'
filename=`basename $filepath`
abspath=$(cd $(dirname $filepath) && pwd)/$(basename $filepath)
downloadfile="$indir/$filename"
mkdir -p $indir
ln -s $abspath $downloadfile
EOF

echo "# List input and position directory..."
perl bin/moirai2.pl -d db -s 1 -m 1 -i '$project->project#indir->$input' -o '$basename->$project/input->$filepath' ls
perl bin/moirai2.pl -d db -s 1 -m 1 -i '$project->project#positiondir->$input' -o '$basename->$project/position->$filepath' ls

echo "# Creating VCF table..."
perl bin/moirai2.pl -d db -s 60 -m 1 -S bin/vcftable.pl -S bin/getCases.pl -i 'root->project->$project,$project->project#indir->$indir,$project->vcftable#qvThreshold->$qvThreshold' -o 'root->$project/vcftable->$output' command 'output=$project/vcftable.txt' << 'EOF'
output=$tmpdir/vcftable.txt
logFile=$project/log/vcftable.txt
vcftable.pl -t $qvThreshold -o $output $indir/* 2> $logFile
echo "INSERT root->$project/vcftable#log->$logFile"
EOF
perl bin/moirai2.pl -d db -s 60 -m 1 -S bin/vcftable.pl -S bin/getCases.pl -i 'root->$project/vcftable->$vcftable' -o 'root->$project/case->$case' command 'output=$project/vcftable.txt' << 'EOF'
case=(`getCases.pl < $vcftable`)
EOF

echo "# Calculating positions by findrun or genomecov..."
#genomecov.pl
perl bin/moirai2.pl -d db -s 60 -m 1 -S bin/genomecov.pl -b '$noIndel:-d' -i '$project->hdr#targetMode->DD,root->$project/vcftable->$table,root->$project/case->$case,$project->coverage#noIndel->$noIndel,$project->coverage#stretchMode->$stretchMode,$project->coverage#topX->$topX,$project->coverage#regionSize->$regionSize,$project->annotation#genome->$genome,$genome->chrominfo->$chrominfo' -o '$case->$project/position->$output' command 'output=$project/findrun/$filename' 'output=$project/findrun/$filename' << 'EOF'
logFile=$project/log/coverage.$case.txt
genomecov.pl $noIndel -o $tmpdir -m $stretchMode -t $topX -r $regionSize $table $case $chrominfo > $output 2> $logFile
output=`ls $tmpdir/*.txt`
filename=`basename $output`
echo "INSERT $case->$project/position#log->$logFile"
EOF

#findrun.pl
perl bin/moirai2.pl -d db -s 60 -m 1 -S bin/findrun.pl -b '$noIndel:-d' -i '$project->hdr#targetMode->DD,root->$project/vcftable->$table,root->$project/case->$case,$project->findrun#noIndel->$noIndel,$project->findrun#stretchMode->$stretchMode,$project->findrun#pickupNumber->$pickupNumber,$project->findrun#regionSize->$regionSize,$project->findrun#skipCount->$skipCount' -o '$case->$project/position->$output' command 'output=$project/findrun/$filename' << 'EOF'
logFile=$project/log/findrun.$case.txt
findrun.pl $noIndel -m $stretchMode -p $pickupNumber -r $regionSize -s $skipCount -o $tmpdir $table $case 2> $logFile
output=`ls $tmpdir/*.txt`
filename=`basename $output`
echo "INSERT $case->$project/position#log->$logFile"
EOF

echo "# Calculating HDR..."
# DD mode calcualte per case
perl bin/moirai2.pl -d db -s 60 -m 1 -S bin/hdr.pl -b '$excludeIndel:-d,$excludeLowQuality:-l' -i '$project->hdr#targetMode->DD,root->$project/vcftable->$table,root->$project/case->$case,$case->$project/position->$position,$project->hdr#excludeIndel->$excludeIndel,$project->hdr#excludeLowQuality->$excludeLowQuality' -o '$case->$project/hdr->$output' command 'output=$project/hdr/$filename' << 'EOF'
logFile=$project/log/hdr.$case.txt
hdr.pl $excludeIndel $excludeLowQuality -t DD -o $tmpdir $table $case $position 2> $logFile
output=`ls $tmpdir/$case/*.txt`
mv $output $tmpdir/.
rmdir $tmpdir/$case
filename=`basename $output`
output=$tmpdir/$filename
echo "INSERT $case->$project/hdr#log->$logFile"
EOF

#AR/AD calcualte per case
perl bin/moirai2.pl -d db -s 60 -m 1 -S bin/hdr.pl -b '$excludeIndel:-d,$excludeLowQuality:-l' -i '$project->hdr#targetMode->$targetMode,root->$project/vcftable->$table,root->$project/case->$case,$case->$project/position->$position,$project->hdr#interval->$interval,$project->hdr#startDistance->$startDistance,$project->hdr#endDistance->$endDistance,$project->hdr#excludeIndel->$excludeIndel,$project->hdr#excludeLowQuality->$excludeLowQuality' -o '$case->$project/hdr->$output' command 'output=$project/hdr/$filename' << 'EOF'
logFile=$project/log/hdr.$case.txt
hdr.pl $excludeIndel $excludeLowQuality -t $targetMode -s $startDistance -e $endDistance -i $interval -o $tmpdir $table $case $position 2> $logFile
output=`ls $tmpdir/$case/*.txt`
mv $output $tmpdir/.
rmdir $tmpdir/$case
filename=`basename $output`
output=$tmpdir/$filename
echo "INSERT $case->$project/hdr#log->$logFile"
EOF

echo "# Calculating statdel/maxStats..."
#AR/AD/DD
perl bin/moirai2.pl -d db -s 10 -m 1 -i 'root->$project/case->$case,$case->$project/hdr->$hdr' -o '$case->$project/stats->$output' -O "Results in file" command 'output=$project/stats/$filename' << 'EOF'
logFile=$project/log/stats.${case}.txt
perl bin/statdel.pl -o $tmpdir $hdr 2> $logFile
output=`ls $tmpdir/$case*`
filename=`basename $output`
echo "INSERT $case->$project/stats#log->$logFile"
EOF

echo "# Parsing statdel/maxStats..."
#AR/AD/DD
perl bin/moirai2.pl -d db -s 10 -m 1 -S bin/annotation.pl -i '$project->annotation#genome->$genome,$genome->chrominfo->$chrominfo,$genome->cytoband->$cytoband,$genome->gencode->$gencode,root->$project/case->$case,$case->$project/stats->$stats' -o '$case->$project/stats#html->$html' command 'html=$project/html/$filename.html' << 'EOF'
annotation.pl $stats $cytoband $gencode $chrominfo js/tablesort-min.js > $html
filename=`basename $stats .out`
EOF
