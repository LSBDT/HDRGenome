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
├── bin/
│   ├── findrun.pl - Uses table from vcftable.pl to compute HDR
│   ├── hdr.pl - Uses table from vcftable.pl to compute HDR
│   ├── moirai2.pl - workflow script
│   ├── rdf.pl - pipeline script
│   ├── statdel - Mac version of stadel (binary)
│   ├── statdel.pl - A wrapper script to run statdel
│   ├── vcftable.pl - Merge multiple VCF files into one table
├── README - this README website
├── run.sh - Pipeline
└── testdata - Test data used in test case
    ├── case/ -  sample VCF files of case
    ├── ctrl/ -  sample VCF files of control
    ├── hdr/ -  cxample of hdr output
    ├── position/
    │   ├── CHR_POS_NA18939_v2.txt - position TSV used by test case
    │   └── Region_NA18939_v2.txt - region TSV used by test case
    └── statdel/ - example of statdel output
```
## URL
  - Own Cloud: https://genomec.gsc.riken.jp/gerg/owncloud/index.php/s/oanbE95OdimomfW
  - Europe PMC: https://europepmc.org/article/med/28722338
  - Nature: https://www.ncbi.nlm.nih.gov/pubmed/26143870
## Script
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
### example
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

### flag
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

### pipeline
- "run.sh" is a pipeline which does following processes:
  - vcftable.pl
  - findrun.pl
  - hdr.pl
  - statdel.pl
#### command
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
bash run.sh [PROJECT_NAME]
```
- A directory with specified PROJECT_NAME will be created.
- All the results will be stored under the project directory.
- Multiple projects can be created under a work directory.

#### parameters
- At the start of a workflow, parameters will be set through prmpts.
- Parameters needed by the pipeline will be entered through prompt questions.
- If you want to change the parameter, type in and then hit return.
- If it's OK with the default value, just hit return.
```
$ bash run.sh
[hdr] target mode {AR|AD|DD} [default=AR]?
[vcftable] Path to input directory [default=input]?
[vcftable] QV threshold for low quality [default=50]?
[hdr] Exclude indel {T|F} [default=F]?T
[hdr] exclude low quality [default=F]?T
```

- Following question will be asked for the DD mode:
```
[findrun] include indel {T|F} [default=F]?
[findrun] stretch mode {hom|het} [default=hom]?
[findrun] pickup number [default=1]?
[findrun] region size [default=1000000]?
[findrun] skip count [default=0]?
```

- Following question will be asked for the AR/AD mode:
```
[hdr] interval [default=10]?
[hdr] start distance [default=10]?50
[hdr] end distance [default=50]?100
```

#### database
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