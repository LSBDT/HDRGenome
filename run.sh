#!/bin/sh
if [ $# -ne 1 ]; then
echo "run.sh [PROJECTNAME]"
exit
fi

# project
projectName=$1
perl bin/moirai2.pl -o "root->project->$projectName" insert

# mode
perl bin/moirai2.pl -o '$project->targetMode->$answer' prompt "\$project=$projectName" '[$project] target mode {AR|AD|DD} [default=AR]?'

# vcftable
perl bin/moirai2.pl -o '$project->indir->$answer' prompt "\$project=$projectName" '[vcftable] Path to input directory [default=testdata/input]?'
perl bin/moirai2.pl -o '$project->qvThreshold->$answer' prompt "\$project=$projectName" '[vcftable] QV threshold for low quality [default=50]?'

# DD findrun
perl bin/moirai2.pl -i '$project->targetMode->DD' -o '$project->noIndel->$answer' prompt "\$project=$projectName" '[findrun] Exclude indel {T|F} [default=F]?'
perl bin/moirai2.pl -i '$project->targetMode->DD' -o '$project->stretchMode->$answer' prompt "\$project=$projectName" '[findrun] stretch mode {hom|het} [default=hom]?'
perl bin/moirai2.pl -i '$project->targetMode->DD' -o '$project->pickupNumber->$answer' prompt "\$project=$projectName" '[findrun] pickup number [default=1]?'
perl bin/moirai2.pl -i '$project->targetMode->DD' -o '$project->regionSize->$answer' prompt "\$project=$projectName" '[findrun] region size [default=1000000]?'
perl bin/moirai2.pl -i '$project->targetMode->DD' -o '$project->skipCount->$answer' prompt "\$project=$projectName" '[findrun] skip count [default=0]?'

# position file
perl bin/moirai2.pl -i '$project->targetMode->AR' -o '$project->positiondir->$answer' prompt "\$project=$projectName" '[position] Path to position files [default=testdata/position/]?'
perl bin/moirai2.pl -i '$project->targetMode->AD' -o '$project->positiondir->$answer' prompt "\$project=$projectName" '[position] Path to position files [default=testdata/position/]?'

# hdr
perl bin/moirai2.pl -o '$project->excludeIndel->$answer' prompt "\$project=$projectName" '[hdr] Exclude indel {T|F} [default=F]?'
perl bin/moirai2.pl -o '$project->excludeLowQuality->$answer' prompt "\$project=$projectName" '[hdr] exclude low quality {T|F} [default=F]?'

# AR hdr
perl bin/moirai2.pl -i '$project->targetMode->AR' -o '$project->interval->$answer' prompt "\$project=$projectName" '[hdr] interval [default=10]?'
perl bin/moirai2.pl -i '$project->targetMode->AR' -o '$project->startDistance->$answer' prompt "\$project=$projectName" '[hdr] start distance [default=10]?'
perl bin/moirai2.pl -i '$project->targetMode->AR' -o '$project->endDistance->$answer' prompt "\$project=$projectName" '[hdr] end distance [default=50]?'

# AD hdr
perl bin/moirai2.pl -i '$project->targetMode->AD' -o '$project->interval->$answer' prompt "\$project=$projectName" '[hdr] interval [default=10]?'
perl bin/moirai2.pl -i '$project->targetMode->AD' -o '$project->startDistance->$answer' prompt "\$project=$projectName" '[hdr] start distance [default=10]?'
perl bin/moirai2.pl -i '$project->targetMode->AD' -o '$project->endDistance->$answer' prompt "\$project=$projectName" '[hdr] end distance [default=50]?'

echo "# Creating project directory..."
perl bin/moirai2.pl -i 'root->project->$project' -o '$project->flag/dircreated->T' command << 'EOF'
mkdir -p $project
EOF

echo "# List input directory..."
perl bin/moirai2.pl -i '$project->indir->$input' -o '$project->case->$basename,$basename->input->$path' ls

echo "# List position directory..."
perl bin/moirai2.pl -i '$project->positiondir->$input' -o '$basename->$project/position->$path' ls

echo "# Creating VCF table..."
perl bin/moirai2.pl -i 'root->project->$project,$project->indir->$directory,$project->qvThreshold->$qvThreshold' -o '$project->vcftable->$output' command 'output=$project/vcftable.txt' << 'EOF'
output=$tmpdir/vcftable.txt
perl bin/vcftable.pl -t $qvThreshold -o $output $directory/*
EOF

echo "# Calculating Findrun..."
# findrun
perl bin/moirai2.pl -b '$noIndel:-d' -i '$project->targetMode->DD,$project->vcftable->$table,$project->case->$case,$case->input->$input,$project->noIndel->$noIndel,$project->stretchMode->$stretchMode,$project->pickupNumber->$pickupNumber,$project->regionSize->$regionSize,$project->skipCount->$skipCount' -o '$case->$project/position->$output' command << 'EOF'
outdir=$project/findrun
perl bin/findrun.pl $noIndel -m $stretchMode -p $pickupNumber -r $regionSize -s $skipCount -o $outdir $table $input
output=`ls $outdir/$case*`
EOF

echo "# Calculating HDR..."
# DD mode
perl bin/moirai2.pl -b '$excludeIndel:-d,$excludeLowQuality:-l' -i '$project->targetMode->DD,$project->vcftable->$table,$project->case->$case,$case->input->$input,$case->$project/position->$position,$project->indir->$control,$project->excludeIndel->$excludeIndel,$project->excludeLowQuality->$excludeLowQuality' -o '$case->$project/hdr->$output' command << 'EOF'
outdir=$project/hdr
perl bin/hdr.pl $excludeIndel $excludeLowQuality -t DD -o $outdir $table $input $control $position
output=`ls $outdir/$case/*`
EOF

#AR/AD
perl bin/moirai2.pl -b '$excludeIndel:-d,$excludeLowQuality:-l' -i '$project->targetMode->$targetMode,$project->vcftable->$table,$project->case->$case,$case->input->$input,$case->$project/position->$position,$project->indir->$control,$project->interval->$interval,$project->startDistance->$startDistance,$project->endDistance->$endDistance,$project->excludeIndel->$excludeIndel,$project->excludeLowQuality->$excludeLowQuality' -o '$case->$project/hdr->$output' command << 'EOF'
outdir=$project/hdr
perl bin/hdr.pl $excludeIndel $excludeLowQuality -t $targetMode -s $startDistance -e $endDistance -i $interval -o $outdir $table $input $control $position
output=`ls $outdir/$case/*`
EOF

echo "# Calculating statdel..."
perl bin/moirai2.pl -i '$case->$project/hdr->$hdr' -o '$case->$project/statdel->$output' -O "Results in file" command << 'EOF'
outdir=$project/statdel
perl bin/statdel.pl -o $outdir $hdr
output=`ls $outdir/$case*`
EOF