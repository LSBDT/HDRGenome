# HDRGenome Project
## Abstract
 Hamming Distance Ratio (HDR)

> High-density oligonucleotide arrays have widely been used to detect pathogenic chromosomal deletions. In addition to high-density oligonucleotide arrays, programs using whole-exome sequencing have become available for estimating copy-number variations using depth of coverage. Here, we propose a new statistical method, HDR-del, to prioritize pathogenic chromosomal deletions based on Hamming distance in exome sequencing. In vcf (variant call format) files generated from exome sequencing, hemizygous chromosomal deletion regions lack heterozygous variants and lead to apparent long runs of homozygosity (ROH). In our Hamming distance ratio (HDR)-del approach, we calculate the "difference" in heterozygous status between an affected individual and control individuals using the HDR over all candidate chromosomal deletion regions defined as ROH longer than 1Mbp. Using a suitable test statistic, which is expected to be large for a true pathogenic deletion region, we prioritize candidate chromosomal deletion regions based on this statistic. In our approach, we were able to considerably narrow down true pathogenic chromosomal deletion regions, which were confirmed by high-density oligonucleotide arrays in four mitochondrial disease patients. Our HDR-del approach represents an easy method for detecting chromosomal deletions.

## Member
- Takeya Kasukawa <takeya.kasukawa@riken.jp>
- Akira Hasegawa <akira.hasegawa@riken.jp>	
- Atsuko Imai <atsuko_imai_hp@yahoo.co.jp>
- Atsushi Kondo <atsushi.kondo@riken.jp>
- Akihiro Nakaya <nakaya@edu.k.u-tokyo.ac.jp>
- Yasushi Okazaki <ya-okazaki@juntendo.ac.jp>

## Structure
```
hdrgenome/
├── annotation/ - files for annotation
│   ├── hg19.chromsizes - chrom sizes of hg19
│   ├── hg19.cytoband.bed - telomere information of hg19
│   ├── hg19.gencode.bed - gencode information of hg19
│   ├── hg38.chrom.sizes  - chrom sizes of hg38
│   ├── hg38.cytoband.bed - telomere information of hg38
│   └── hg38.gencode.bed - gencode information of hg38
├── bin/
│   ├── annotation.pl - Annotates statistical HDR results with overlapping telomere and genes.
│   ├── dag.pl -Workflow database script
│   ├── findrun.pl - Uses table from vcftable.pl to compute HDR
│   ├── genomecov.pl - List up candidates with genome coverage
│   ├── getCases.pl - Get patient/controls from vcftable
│   ├── hdr.pl - Uses table from vcftable.pl to compute HDR
│   ├── linux/ - linux binary directory
│   │   ├── maxstatRS - Linux version of stadel (binary)
│   │   └── statdel - Linux version of stadel (binary)
│   ├── mac/ - Mac binary directory
│   │   ├── maxstatRS - MacOS version of stadel (binary)
│   │   └── statdel - MacOS version of stadel (binary)
│   ├── moirai2.pl - workflow script
│   ├── selectGenesFromGencode.pl - Select genes from gencode file
│   ├── selctVcftable.pl - Count number of hits per vcf regions
│   ├── statdel.pl - A wrapper script to run statdel
│   ├── vcftable.pl - Merge multiple VCF files into one table
│   └── vcftable2bed.pl - Convert vcf table to bed format
├── docker.sh - Script to run pipeline with a docker.
├── hdr.sh - Script to run pipeline.
├── LINCENCE - LICENCE of this project
├── README - README of this project
└── testdata - Test data used in test case
    ├── configAD.txt - config file used by a pipeline for test data
    ├── configAD_run_with_vcftable_from_configAD.txt - config file used by a pipeline for test data
    ├── configAR.txt - config file used by a pipeline for test data
    ├── configDD.txt - config file used by a pipeline for test data
    ├── download1000GenomesChromosomeSamples.sh - download sample chromosome VCF files from 1000Genomes
    ├── download1000GenomesExonSamples.sh - download sample exon VCF files from 1000Genomes
    ├── input/ -  sample VCF files of case/control
    │   ├── NA18939_v2.vcf - variant call format of NA18939_v2 patient
    │   ├── NA18940_v2.vcf - variant call format of NA18940_v2
    │   ├── NA18941_v2.vcf - variant call format of NA18941_v2
    │   └── NA18942_v2.vcf - variant call format of NA18942_v2
    └── position/
        ├── NA18939_v2.txt - region to investigate for test case NA18939_v2
        ├── NA18940_v2.txt - region to investigate for test case NA18940_v2
        ├── NA18941_v2.txt - region to investigate for test case NA18941_v2
        └── NA18942_v2.txt - region to investigate for test case NA18942_v2
```
## URL
- Europe PMC: https://europepmc.org/article/med/28722338
- Nature: https://www.ncbi.nlm.nih.gov/pubmed/26143870

## Requirement
- MaxOS or Linux environment
- perl (v5 or more)
- Docker: https://www.docker.com

## Automation Pipeline
- "hdr.sh" is a pipeline which does following processes:
  - vcftable.pl
  - findrun.pl / genomecov.pl [DD Mode]
  - hdr.pl
  - statdel.pl
  - annotation.pl
- Use "docker.sh" instead of "hdr.sh", if you want to run with docker.
- Please pull "moirai2/biotools" from docker hub with
```
docker pull moirai2/biotools
```
- "hdr.sh" is a pipeline which does following processes:
- Place vcf|bcf|avinput under a directory (for example, input/ directory)
```
hdrgenome/
└──input/
    ├── PatientA.avinput
    ├── PatientB.avinput
    └── PatientC.avinput
```
- Under hdrgenome root directory, start pipline with the following command line.
```
bash hdr.sh [CONFIG]
OR
bash docker.sh [CONFIG]
```
- A directory with specified PROJECT_NAME will be created.
- All the results will be stored under the project directory.
- Multiple projects can be created under a work directory.
- To run test samples:
```
bash hdr.sh testdata/configAD.txt
bash hdr.sh testdata/configAR.txt
bash hdr.sh testdata/configDD.txt
```
- Directory will be created with a project name.
- testdataAD for example will create these directories and files:
```
testdataAD/
├──hdr/  Result from HDR computation
├──log/
|   ├── hdr.txt  Log from HDR computation
|   ├── stats.txt  Log from statdel/maxStats computation
|   └── vcftable.txt  Log from VCFTable computation
├──stats/  Result from statdel/maxStats computation
└──vcftable.txt  Result from VCFTable computation
```
- 'configAD_run_with_vcftable_from_configAD.txt' reuses vcftable from configAD.txt by specifying vcftable information in config file.
```
$project->vcftable	testAD/vcftable.txt
```
- If you don't want to recompute vcftable.txt, adding this line in config file will save time.
- You can't RERUN the pipeline with the same project name.
- Moirai directory used for pipeline control and pipeline database.
- When something goes wrong in the pipeline, check HDRGenome/moirai/log/error directory.
- HDRGenome/moirai/db directory is used to store pipeline status in triple database.
```
HDRGenome/
└─-moirai/
    ├──ctrl/  Directory used to control pipeline process.
    ├──db/  Database of pipeline in triple
    └──log/  Result from HDR computation
        ├──check/  Directory to store checking of process.
        ├──error/  Directory to store error processes.
        └──json/  Directory to store commands of pipeline information in json
```
- If you delete moirai/db/PROJECT directory, you can rerun the pipeline again.
- If you just want to run specific step again:
  - Remove moirai/db/PROJECT/vcftable.txt to recompute VCFTable step.
  - Remove moirai/db/PROJECT/hdr.txt to recompute HDR step.
  - Remove moirai/db/PROJECT/stats.txt to recompute statdel/maxStats step.
### Config File
- Line starting with '#' is a comment line.
- Config lines are separated by a tab.
```
#project
$project	test
$project->targetMode	AR
root->project	$project
#VCF table
$project->indir		testdata/input
$project->qvThreshold	50
$project->noIndel	F
#$project->vcftable	testdata/vcftable.txt
#HDR
$project->excludeIndel	F
$project->excludeLowQuality	F
#DD mode
$project->stretchMode	hom
$project->pickupNumber	1
$project->regionSize	1000000
$project->skipCount	0
#AR|AD mode
$project->positiondir	testdata/position
$project->interval	10
$project->startDistance	10
$project->endDistance	50
```
- Example of config files can be found at testdata/
  - testdata/configAD.txt  Example of config for AD mode
  - testdata/configAD_run_with_vcftable_from_configAD.txt  Example of config for AD mode using already calculated VCF table.
  - testdata/configAR.txt  Example of config for AR mode
  - testdata/configDD.txt  Example of config for DD mode
- By specifying path to VCFtable file, VCFtable computation step will be  
```
$project->vcftable	testAD/vcftable.txt
```

## Scripts
- Here are instructions on how to use each scripts.

### vcftable.pl
```
Command: vcftable.pl [option] VCF [VCF2 ..]
Arguments:
   VCF  variant call format directory/files
Options:
  -o  Output file (default='STDOUT')
  -t  Quality value threshold (default='40')
Note:
    If you are using BCF files, please install bcftools
    http://samtools.github.io/bcftools/bcftools.html
Flag:
     0  wild
     1  hetero
     2  homo
     4  deletion
     8  insertion
    16  multi allelic (column5.count(',')>=2)
    32  low quality (column6<QV40)
```

### hdr.pl
```
Command: hdr.pl [OPTIONS] CASE CTRL POS
Arguments:
  CASE  Case file
  CTRL  Control files/directory
   POS  Position/region file
Options:
  -d  Indel included/excluded (default='included')
  -e  end distance (default=50)
  -i  Interval (default=10)
  -o  outdir (default="out")
  -s  start distance (default=10)
  -t  Target mode (default='AR')
Mode:
    AR  hom vs non-hom (1/1 vs 1/0,0/0) both HOM=0
    DD  het vs non-het (1/0 vs 1/1,0/0) both HET=0
    AD  all vs all     (1/1 vs 1/0,0/0) exact match=0
Note:
    If you are using BCF files, please install bcftools
    http://samtools.github.io/bcftools/bcftools.html
```

### statdel.pl
- statdel and maxstatRS needs config files which specify parameters and input/output files.
- statdel.pl is a wrapper script which creates parameters automatically.
  - Removes all non normal chromosomes (random chromosomes are removed)
  - Convert chromosome notation to maxstatRS and statdel programs.
    - chrX => chr23
    - chrY => chr24
    - chrM => chr25
    - These notation changes will be fixed back to normal notation after statdel/maxstatRS steps by the script.
  - statdel.pl checks header of input file and decides which program to use.
    - #Chr, Start, End, HDR:1/2 - DD mode - uses statdel
    - #Chr, Position, RegionSize(kb), HDR:1/2 - AD/AR mode - uses maxstatRS 
- statdel config generated by script looks like this:
```
Statdel: Auto generated parameter file
-9 0 0 0
-12 1 1
1 3 0.5 0.8
$input
$output
```
- maxstatRS config generated by script looks like this:
```
maxstatRS: Auto generated parameter file
-9 0 0 0
-12 1 1
$MIN $MAX 3.0
$input
$output
$REGION1
$REGION2
$REGION3
$REGION4
$REGION5
```

## Example
- To create a table from multiple VCF files.
- You can specify directory with VCF files for argument too.
```
perl bin/vcftable.pl testdata/case testdata/ctrl > testdata/table.txt
```
- output
```
#chromosome	position	NA18939_v2	NA18940_v2	NA18941_v2	NA18942_v2	NA18943_v2	NA18944_v2	NA18945_v2	NA18946_v2	NA18947_v2	NA18948_v2	NA18949_v2	NA18950_v2	NA18951_v2	NA18952_v2	NA18953_v2	NA18954_v2	NA18956_v2	NA18957_v2	NA18959_v2	NA18960_v2	NA18961_v2
chr1	809756	0	0	0	1	0	0	1	0	0	1	0	0	0	0	0	0	0	0	0	2	0
chr1	871215	0	1	2	1	0	0	0	0	0	1	0	0	0	1	0	0	1	0	0	1	0
chr1	871269	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0
chr1	874826	0	0	0	0	0	0	1	0	0	0	0	0	1	0	0	0	0	0	0	0	0
chr1	876499	0	2	2	2	2	2	2	0	2	0	2	0	2	2	2	0	2	0	1	2	0
chr1	877715	0	2	2	2	2	2	2	0	2	2	2	0	2	2	2	0	2	0	1	2	0
chr1	877782	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	2	2	0
chr1	879317	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	2	0	0	0	0	0
chr1	880238	0	2	2	2	2	2	2	2	2	2	2	0	2	2	2	2	2	2	1	2	2
chr1	880716	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
chr1	880906	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
```
- To compute HDR values in different modes:
- AR mode:
```
perl bin/hdr.pl testdata/table.txt testdata/case testdata/ctrl testdata/position/CHR_POS_NA18939_v2.txt
```
- Output:
```
387 3 5 # 1:chr1
#Chr	Position	RegionSize(kb)	HDR:1/2	Case-Ctrl1	Case-Ctrl2	Case-Ctrl3	Ctrl1-Ctrl2	Ctrl1-Ctrl3	Ctrl2-Ctrl3
chr1	948921	10	1	0.000	0.000	0.000	0.000	0.000	0.000
chr1	949654	10	1	0.000	0.000	0.000	0.000	0.000	0.000
chr1	949925	10	1	0.000	0.000	0.000	0.000	0.000	0.000
chr1	981087	10	1	0.667	0.667	0.667	0.000	0.000	0.000
chr1	982941	10	1	0.571	0.571	0.571	0.000	0.000	0.000
chr1	982994	10	1	0.571	0.571	0.571	0.000	0.000	0.000
chr1	1120431	10	1	1.000	1.000	1.000	0.000	0.000	0.000
chr1	1120536	10	1	1.000	1.000	1.000	0.000	0.000	0.000
chr1	1138913	10	1	1.000	1.000	1.000	0.000	0.000	0.000
chr1	1198618	10	1	1.000	0.000	0.000	1.000	1.000	0.000
chr1	1237538	10	1	1.000	0.500	1.000	1.000	0.000	1.000
chr1	1246004	10	1	0.000	0.000	0.000	0.000	0.000	0.000
chr1	1254841	10	1	0.000	0.000	0.000	0.000	0.000	0.000
chr1	1263144	10	1	0.500	0.500	0.500	0.000	0.000	0.000
chr1	1268847	10	1	0.000	0.000	0.000	0.000	0.000	0.000
chr1	1269554	10	1	0.000	0.000	0.000	0.000	0.000	0.000
chr1	1310074	10	1	0.500	0.500	0.500	0.000	0.000	0.000
chr1	1330726	10	1	0.000	0.500	0.000	0.500	0.000	0.500
chr1	1333598	10	1	0.000	0.333	0.333	0.333	0.333	0.667
```
- DD mode:
```
perl bin/hdr.pl -t DD testdata/table.txt testdata/case testdata/ctrl testdata/position/Region_NA18939_v2.txt
```
- Output:
```
6 3 1 # 1:chr1
#Chr	Start	End	HDR:1/2	Case-Ctrl1	Case-Ctrl2	Case-Ctrl3	Ctrl1-Ctrl2	Ctrl1-Ctrl3	Ctrl2-Ctrl3
chr1	1237538	1337334	1	1.000	1.000	1.000	1.000	0.600	1.000
chr1	16342103	17249334	1	1.000	1.000	1.000	0.920	0.708	0.875
chr1	17706297	18023509	1	1.000	1.000	1.000	0.571	0.600	0.783
chr1	22183739	22447746	1	1.000	1.000	1.000	1.000	1.000	0.667
chr1	24413118	24495968	1	1.000	1.000	1.000	1.000	1.000	1.000
chr1	25889539	26378211	1	1.000	1.000	1.000	0.833	0.895	0.500
chr1	1237538	1337334	2	3	2	4	5	5	6
chr1	16342103	17249334	2	20	7	11	25	24	16
chr1	17706297	18023509	2	10	10	18	14	20	23
chr1	22183739	22447746	2	1	1	3	2	4	3
chr1	24413118	24495968	2	1	1	5	2	6	6
chr1	25889539	26378211	2	16	7	7	18	19	8
```
- AD mode:
```
perl bin/hdr.pl -t AD testdata/table.txt testdata/case testdata/ctrl testdata/position/CHR_POS_NA18939_v2.txt
```
- Output:
```
387 3 5 # 1:chr1
#Chr	Position	RegionSize(kb)	HDR:1/2	Case-Ctrl1	Case-Ctrl2	Case-Ctrl3	Ctrl1-Ctrl2	Ctrl1-Ctrl3	Ctrl2-Ctrl3
chr1	948921	10	1	0.000	0.000	0.000	0.000	0.000	0.000
chr1	949654	10	1	0.000	0.000	0.000	0.000	0.000	0.000
chr1	949925	10	1	0.000	0.000	0.000	0.000	0.000	0.000
chr1	981087	10	1	0.667	0.667	0.727	0.000	0.182	0.182
chr1	982941	10	1	0.571	0.571	0.667	0.000	0.222	0.222
chr1	982994	10	1	0.571	0.571	0.667	0.000	0.222	0.222
chr1	1120431	10	1	1.000	1.000	1.000	0.000	0.000	0.000
chr1	1120536	10	1	1.000	1.000	1.000	0.000	0.000	0.000
chr1	1138913	10	1	1.000	1.000	1.000	0.000	0.000	0.000
chr1	1198618	10	1	1.000	0.000	0.000	1.000	1.000	0.000
chr1	1237538	10	1	0.750	0.500	1.000	0.750	0.400	1.000
chr1	1246004	10	1	0.000	0.000	0.500	0.000	0.500	0.500
chr1	1254841	10	1	0.000	0.000	0.000	0.000	0.000	0.000
chr1	1263144	10	1	0.500	0.500	0.500	0.000	0.000	0.000
chr1	1268847	10	1	0.333	0.333	0.000	0.500	0.333	0.333
chr1	1269554	10	1	0.333	0.333	0.000	0.500	0.333	0.333
chr1	1310074	10	1	0.500	0.500	0.500	0.000	0.000	0.000
chr1	1330726	10	1	0.000	0.500	0.000	0.500	0.000	0.500
chr1	1333598	10	1	0.000	0.333	0.333	0.333	0.333	0.667
```

## flag
- Output from vcftable.pl
| flag | description |
----|----
| 0 | wild |
| 1 | hetero |
| 2 | homo |
| 4 | deletion |
| 8 | insertion |
| 16 | multi allelic (column5.count(',')>=2) |
| 32 | low quality (column6<QV40) |

| flag | description |
----|----
| 0 | wild |
| 1 | hetero |
| 2 | homo |
| 3 | |
| 4 | deletion, wild |
| 5 | deletion, hetero |
| 6 | deletion, homo |
| 7 | |
| 8 | insertion, wild |
| 9 | insertion, hetero |
| 10 | insertion, homo |
| 11 | |
| 12 | |
| 13 | |
| 14 | |
| 15 | |
| 16 | multi alleric, wild |
| 17 | multi alleric, hetero |
| 18 | multi alleric, homo |
| 19 | |
| 20 | multi alleric, deletion, wild |
| 21 | multi alleric, deletion, hetero |
| 22 | multi alleric, deletion, homo |
| 23 | |
| 24 | multi alleric, insertion, wild |
| 25 | multi alleric, insertion, hetero |
| 26 | multi alleric, insertion, homo |
| 27 | |
| 28 | |
| 29 | |
| 30 | |
| 31 | |
| 32 | low quality, wild |
| 33 | low quality, hetero |
| 34 | low quality, homo |
| 35 | |
| 36 | low quality, deletion, wild |
| 37 | low quality, deletion, hetero |
| 38 | low quality, deletion, homo |
| 39 | |
| 40 | low quality, insertion, wild |
| 41 | low quality, insertion, hetero |
| 42 | low quality, insertion, homo |
| 43 | |
| 44 | |
| 45 | |
| 46 | |
| 47 | |
| 48 | low quality, multi alleric, wild |
| 49 | low quality, multi alleric, hetero |
| 50 | low quality, multi alleric, homo |
| 51 | |
| 52 | low quality, multi alleric, deletion, wild |
| 53 | low quality, multi alleric, deletion, hetero |
| 54 | low quality, multi alleric, deletion, homo |
| 55 | |
| 56 | low quality, multi alleric, insertion, wild |
| 57 | low quality, multi alleric, insertion, hetero |
| 58 | low quality, multi alleric, insertion, homo |
| 59 | |
| 60 | |
| 61 | |
| 62 | |
| 63 | |

### database
- Database are stored in triple (subject->predicate->object)
```
hdrgenome/
└── moirai/
    ├── ctrl/ - Used by moirai2 to control command processes
    ├── db/ - database of the pipeline in (triple format)
    └── log/ - log of command processes are are kept here
        ├── YYYYMMDD/ - processing/completed Logs will be stored here
        ├── error/ - error logs will be stored here
        └── json/ - json files describing commands will be stored here
```
- You can edit the file under db directly through text editor.
- Database are in triple format (subject->predicate->object).
- Predicate is a filename of a text.
- Subject and object are defined in the text and separted by a tab.
```
$ cat moirai/db/stretchMode.txt
hdr hom
```
- This is a triple of "hdr->stretchMode->hom" ('hdr' is a root of a triple tree).
- By changing hom to het, the tiple becomes "hdr->stretchMode->het".
- If vcftable is already created, you can skip the vcftable.pl process by placing a "vcftable.txt" table under moirai/db/ directory.
```
$ cat moirai/db/stretchMode.txt
hdr	hdr/vcftable.txt
```

## Compiling Pascal Scripts
- statdel and maxstatRS are written in Pascal and there is a need to compile them.
- We prepared compiled statdel and maxstatRS and stored them at bin/max and bin/linux directories.
- If you want to compile by yourself, here are the steps:

### mac
- To compile with Mac, use free pascal compiler
- Mac:https://sourceforge.net/projects/freepascal/files/Mac%20OS%20X/3.2.0/
```
cd maxstatRS/
fpc maxstatRS.pas
cd statdel/
fpc statdel.pas
```

#### compile error
- This is a note on what I did.
- When I was compiling with fpc, error occured:
```
Free Pascal Compiler version 3.2.0 [2020/05/31] for x86_64
Copyright (c) 1993-2020 by Florian Klaempfl and others
Target OS: Darwin for x86_64
Compiling maxstatRS.pas
errtrap.p(3,13) Error: Duplicate identifier "exitsave"
maxstatRS.pas(1245) Fatal: There were 1 errors compiling module, stopping
Fatal: Compilation aborted
Error: /usr/local/bin/ppcx64 returned an error exitcode
```
- In statdel.pas, this line is not declared, so I removed.
```
var exitsave:pointer;
```

### linux
- Compiling with free pascal docker image.
- https://github.com/Docker-Hub-frolvlad/docker-alpine-fpc
```
cd maxstatRS/
docker run --rm -v `pwd`:/tmp frolvlad/alpine-fpc fpc /tmp/maxstatRS.pas
cd statdel/
docker run --rm -v `pwd`:/tmp frolvlad/alpine-fpc fpc /tmp/statdel.pas
```