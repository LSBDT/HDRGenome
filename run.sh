#!/bin/sh
projectName=$1
if [ "$projectName" == "" ]; then
projectName=hdr
fi
perl bin/moirai2.pl -o "hdr->project->$projectName" assign
perl bin/moirai2.pl -i 'hdr->project->$project' -o '$project->targetMode->[hdr] target mode {AR|AD|DD} [default=AR]? ' prompt
perl bin/moirai2.pl -i 'hdr->project->$project' -o '$project->indir->[vcftable] Path to input directory [default=testdata/input]? ' prompt
perl bin/moirai2.pl -i '$project->indir->$input' -o '$project->input->$path' ls
perl bin/moirai2.pl -i 'hdr->project->$project' -o '$project->qvThreshold->[vcftable] QV threshold for low quality [default=50]? ' prompt

perl bin/moirai2.pl -i '$project->targetMode->DD' -o '$project->noIndel->[findrun] include indel {T|F} [default=F]?' prompt
perl bin/moirai2.pl -i '$project->targetMode->DD' -o '$project->stretchMode->[findrun] stretch mode {hom|het} [default=hom]?' prompt
perl bin/moirai2.pl -i '$project->targetMode->DD' -o '$project->pickupNumber->[findrun] pickup number [default=1]?' prompt
perl bin/moirai2.pl -i '$project->targetMode->DD' -o '$project->regionSize->[findrun] region size [default=1000000]?' prompt
perl bin/moirai2.pl -i '$project->targetMode->DD' -o '$project->skipCount->[findrun] skip count [default=0]?' prompt

echo "# Creating project directory..."
perl bin/moirai2.pl -i 'hdr->project->$project' -o '$project->flag/dircreated->T' command << 'EOF'
mkdir -p $project
EOF

echo "# Creating VCF table..."
perl bin/moirai2.pl -i 'hdr->project->$project,$project->indir->$directory,$project->qvThreshold->$qvThreshold' -o '$project->vcftable->$output' command 'output=$project/vcftable.txt' << 'EOF'
output=$tmpdir/vcftable.txt
perl bin/vcftable.pl -t $qvThreshold -o $output $directory/*
EOF

echo "# Calculating Findrun..."
# findrun
perl bin/moirai2.pl -i '$project->targetMode->DD,$project->vcftable->$table,$project->input->$input,$project->noIndel->$noIndel,$project->stretchMode->$stretchMode,$project->pickupNumber->$pickupNumber,$project->regionSize->$regionSize,$project->skipCount->$skipCount' -o '$input->findrun->$output' command << 'EOF'
outdir=$project/findrun
basename=`basename $input`
basename=${basename%.*}
perl bin/findrun.pl -m $stretchMode -p $pickupNumber -r $regionSize -s $skipCount -o $outdir $table $input
output=`ls $outdir/$basename*`
EOF

echo "# Calculating HDR..."
perl bin/moirai2.pl -i '$project->targetMode->$targetMode,$project->vcftable->$table,$project->input->$case,$project->indir->$control,$case->findrun->$findrun' -o '$case->hdr->$output' command << 'EOF'
outdir=$project/hdr
basename=`basename $case`
basename=${basename%.*}
perl bin/hdr.pl -t $targetMode -o $outdir $table $case $control $findrun
output=`ls $outdir/$basename/*`
EOF

echo "# Calculating statdel..."
perl bin/moirai2.pl -i '$project->input->$input,$input->hdr->$output' -o '$input->statdel->$output' command << 'EOF'
outdir=$project/statdel
basename=`basename $input`
basename=${basename%.*}
perl bin/statdel.pl -o $outdir $input
output=`ls $outdir/$basename*`
EOF


exit

perl bin/moirai2.pl -i 'hdr->mode->AR' -o 'AR->excludeIndel->[HDR] Exclude indel {T|F} [default=F]?' prompt
perl bin/moirai2.pl -i 'hdr->mode->AR' -o 'AR->interval->[HDR] interval [default=10]?' prompt
perl bin/moirai2.pl -i 'hdr->mode->AR' -o 'AR->excludeLowQuality->[HDR] exclude low quality [default=F]?' prompt
perl bin/moirai2.pl -i 'hdr->mode->AR' -o 'AR->startDistance->[HDR] start distance [default=10]?' prompt
perl bin/moirai2.pl -i 'hdr->mode->AR' -o 'AR->endDistance->[HDR] end distance [default=50]?' prompt
perl bin/moirai2.pl -i 'hdr->mode->AR' -o 'AR->positionFile->Path to position file?' prompt

exit
