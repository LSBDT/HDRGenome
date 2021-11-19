#!/bin/sh
if [ $# -ne 1 ]; then
echo "run.sh [CONFIG]"
exit
fi

bin/rdf.pl -q config $1

echo "# Creating project directory..."
perl bin/moirai2.pl -m 1 -i 'root->project->$project' -o '$project->flag/dircreated->T' command << 'EOF'
mkdir -p $project
mkdir -p $project/log
EOF

echo "# List input directory..."
perl bin/moirai2.pl -x -m 1 -i '$project->indir->$input' -o '$project->case->$basename,$basename->input->$path' ls

echo "# List position directory..."
perl bin/moirai2.pl -x -m 1 -i '$project->positiondir->$input' -o '$basename->$project/position->$path' ls

echo "# Creating VCF table..."
perl bin/moirai2.pl -m 1 -i 'root->project->$project,$project->indir->$directory,$project->qvThreshold->$qvThreshold' -o '$project->vcftable->$output,$project->vcfTableLog->$logFile' command 'output=$project/vcftable.txt' << 'EOF'
output=$tmpdir/vcftable.txt
logFile=$project/log/vcftable.txt
perl bin/vcftable.pl -t $qvThreshold -o $output $directory/* 2> $logFile
EOF

echo "# Calculating Findrun..."
# findrun
perl bin/moirai2.pl -m 1 -b '$noIndel:-d' -i '$project->targetMode->DD,$project->vcftable->$table,$project->case->$case,$case->input->$input,$project->noIndel->$noIndel,$project->stretchMode->$stretchMode,$project->pickupNumber->$pickupNumber,$project->regionSize->$regionSize,$project->skipCount->$skipCount' -o '$case->$project/position->$output,$case->$project/findRunLog->$logFile' command << 'EOF'
outdir=$project/findrun
logFile=$project/log/findrun.txt
perl bin/findrun.pl $noIndel -m $stretchMode -p $pickupNumber -r $regionSize -s $skipCount -o $outdir $table $input 2> $logFile
output=`ls $outdir/$case*`
EOF

echo "# Calculating HDR..."
# DD mode
perl bin/moirai2.pl -m 1 -b '$excludeIndel:-d,$excludeLowQuality:-l' -i '$project->targetMode->DD,$project->vcftable->$table,$project->case->$case,$case->input->$input,$case->$project/position->$position,$project->indir->$control,$project->excludeIndel->$excludeIndel,$project->excludeLowQuality->$excludeLowQuality' -o '$case->$project/hdr->$output,$case->$project/hdrLog->$logFile' command << 'EOF'
outdir=$project/hdr
logFile=$project/log/hdr.txt
perl bin/hdr.pl $excludeIndel $excludeLowQuality -t DD -o $outdir $table $input $control $position 2> $logFile
output=`ls $outdir/$case/*`
EOF

#AR/AD
perl bin/moirai2.pl -m 1 -b '$excludeIndel:-d,$excludeLowQuality:-l' -i '$project->targetMode->$targetMode,$project->vcftable->$table,$project->case->$case,$case->input->$input,$case->$project/position->$position,$project->indir->$control,$project->interval->$interval,$project->startDistance->$startDistance,$project->endDistance->$endDistance,$project->excludeIndel->$excludeIndel,$project->excludeLowQuality->$excludeLowQuality' -o '$case->$project/hdr->$output,$case->$project/statsLog->$logFile' command << 'EOF'
outdir=$project/hdr
logFile=$project/log/hdr.txt
perl bin/hdr.pl $excludeIndel $excludeLowQuality -t $targetMode -s $startDistance -e $endDistance -i $interval -o $outdir $table $input $control $position 2> $logFile
output=`ls $outdir/$case/*`
EOF
echo "# Calculating statdel/maxStats..."
#AR/AD/DD
perl bin/moirai2.pl -m 1 -i '$project->case->$case,$case->$project/hdr->$hdr' -o '$case->$project/stats->$output,$case->$project/statsLog->$logFile' -O "Results in file" command << 'EOF'
outdir=$project/stats
logFile=$project/log/stats.txt
perl bin/statdel.pl -o $outdir $hdr 2> $logFile
output=`ls $outdir/$case*`
EOF
