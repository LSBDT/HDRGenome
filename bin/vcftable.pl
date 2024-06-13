#!/usr/bin/perl
use strict 'vars';
use Cwd;
use File::Basename;
use File::Temp qw/tempfile tempdir/;
use FileHandle;
use Getopt::Std;
use IO::File;
use Time::localtime;
############################## HEADER ##############################
my ($program_name,$program_directory,$program_suffix)=fileparse($0);
$program_directory=Cwd::abs_path($program_directory);
my $program_path="$program_directory/$program_name";
my $program_version="2024/06/13";
############################## OPTIONS ##############################
use vars qw($opt_a $opt_c $opt_C $opt_h $opt_o $opt_q $opt_s $opt_S $opt_t);
getopts('a:c:C:ho:qs:S:t:');
############################## HELP ##############################
sub help{
	print "\n";
	print "Command: $program_name [option] VCF [VCF2 ..]\n";
	print "Arguments:\n";
	print "   VCF  variant call format (VCF or GVCF or avinput) files or under directories\n";
	print "Options:\n";
	print "    -a  Alpha Missense file to add new 1024 flag(default=none)\n";
	print "    -c  Do coverage depth analysis with low/high thresholds (default=0.5,2.0)\n";
	print "    -C  Write coverage depth to this file (default=none)\n";
	print "    -h  Show help\n";
	print "    -o  Output file (default='STDOUT')\n";
	print "    -q  Quiet\n";
	print "    -s  Do segmental duplication analysis with low/high thresholds (default=0.5,1.5)\n";
	print "    -S  Write ref/alt ratio to this file (default=none)\n";
	print "    -t  Quality value threshold for low quality (default<50)\n";
	print "Flag:\n";
	print "     0  wild\n";
	print "     1  (het)erozygous\n";
	print "     2  (hom)ozygous\n";
	print "     4  deletion\n";
	print "     8  insertion\n";
	print "    16  multi allelic (column5.count(',')>0)\n";
	print "    32  low quality (column6<QV50)\n";
	print "    64  segmental deletion candidate (ALT/REF<1/threshold=0.5)\n";
	print "   128  segmental duplication candidate (ALT/REF>=threshold=1.5)\n";
	print "   256  low coverage candidate (DEPTH/AVG<=1/threshold=0.5)\n";
	print "   512  high coverage candidate (DEPTH/AVG>=threshold=2.0)\n";
	print "  1024  'likely pathogenic' from Deepmind AlphaMissense software\n";
	print "\n";
	print "Combinations:\n";
	print "   1 = het(1)\n";
	print "   2 = hom(2)\n";
	print "   5 = het(1)+del(4)\n";
	print "   6 = hom(2)+del(4)\n";
	print "   9 = het(1)+ins(8)\n";
	print "  10 = hom(2)+ins(8)\n";
	print "  17 = het(1)+multi(16)\n";
	print "  18 = hom(2)+multi(16)\n";
	print "  21 = het(1)+del(4)+multi(16)\n";
	print "  22 = hom(2)+del(4)+multi(16)\n";
	print "  25 = het(1)+ins(8)+multi(16)\n";
	print "  26 = hom(2)+ins(8)+multi(16)\n";
	print "  33 = het(1)+lowQual(32)\n";
	print "  34 = hom(2)+lowQual(32)\n";
	print "  37 = het(1)+del(4)+lowQual(32)\n";
	print "  38 = hom(2)+del(4)+lowQual(32)\n";
	print "  41 = het(1)+ins(8)+lowQual(32)\n";
	print "  42 = hom(2)+ins(8)+lowQual(32)\n";
	print "  49 = het(1)+multi(16)+lowQual(32)\n";
	print "  50 = hom(2)+multi(16)+lowQual(32)\n";
	print "  53 = het(1)+del(4)+multi(16)+lowQual(32)\n";
	print "  54 = hom(2)+del(4)+multi(16)+lowQual(32)\n";
	print "  57 = het(1)+ins(8)+multi(16)+lowQual(32)\n";
	print "  58 = hom(2)+ins(8)+multi(16)+lowQual(32)\n";
	print "  65 = het(1)+segDel(64)\n";
	print "  66 = hom(2)+segDel(64)\n";
	print "  69 = het(1)+del(4)+segDel(64)\n";
	print "  70 = hom(2)+del(4)+segDel(64)\n";
	print "  73 = het(1)+ins(8)+segDel(64)\n";
	print "  74 = hom(2)+ins(8)+segDel(64)\n";
	print "  81 = het(1)+multi(16)+segDel(64)\n";
	print "  82 = hom(2)+multi(16)+segDel(64)\n";
	print "  85 = het(1)+del(4)+multi(16)+segDel(64)\n";
	print "  86 = hom(2)+del(4)+multi(16)+segDel(64)\n";
	print "  89 = het(1)+ins(8)+multi(16)+segDel(64)\n";
	print "  90 = hom(2)+ins(8)+multi(16)+segDel(64)\n";
	print "  97 = het(1)+lowQual(32)+segDel(64)\n";
	print "  98 = hom(2)+lowQual(32)+segDel(64)\n";
	print " 101 = het(1)+del(4)+lowQual(32)+segDel(64)\n";
	print " 102 = hom(2)+del(4)+lowQual(32)+segDel(64)\n";
	print " 105 = het(1)+ins(8)+lowQual(32)+segDel(64)\n";
	print " 106 = hom(2)+ins(8)+lowQual(32)+segDel(64)\n";
	print " 113 = het(1)+multi(16)+lowQual(32)+segDel(64)\n";
	print " 114 = hom(2)+multi(16)+lowQual(32)+segDel(64)\n";
	print " 117 = het(1)+del(4)+multi(16)+lowQual(32)+segDel(64)\n";
	print " 118 = hom(2)+del(4)+multi(16)+lowQual(32)+segDel(64)\n";
	print " 121 = het(1)+ins(8)+multi(16)+lowQual(32)+segDel(64)\n";
	print " 122 = hom(2)+ins(8)+multi(16)+lowQual(32)+segDel(64)\n";
	print " 129 = het(1)+segDup(128)\n";
	print " 130 = hom(2)+segDup(128)\n";
	print " 133 = het(1)+del(4)+segDup(128)\n";
	print " 134 = hom(2)+del(4)+segDup(128)\n";
	print " 137 = het(1)+ins(8)+segDup(128)\n";
	print " 138 = hom(2)+ins(8)+segDup(128)\n";
	print " 145 = het(1)+multi(16)+segDup(128)\n";
	print " 146 = hom(2)+multi(16)+segDup(128)\n";
	print " 149 = het(1)+del(4)+multi(16)+segDup(128)\n";
	print " 150 = hom(2)+del(4)+multi(16)+segDup(128)\n";
	print " 153 = het(1)+ins(8)+multi(16)+segDup(128)\n";
	print " 154 = hom(2)+ins(8)+multi(16)+segDup(128)\n";
	print " 161 = het(1)+lowQual(32)+segDup(128)\n";
	print " 162 = hom(2)+lowQual(32)+segDup(128)\n";
	print " 165 = het(1)+del(4)+lowQual(32)+segDup(128)\n";
	print " 166 = hom(2)+del(4)+lowQual(32)+segDup(128)\n";
	print " 169 = het(1)+ins(8)+lowQual(32)+segDup(128)\n";
	print " 170 = hom(2)+ins(8)+lowQual(32)+segDup(128)\n";
	print " 177 = het(1)+multi(16)+lowQual(32)+segDup(128)\n";
	print " 178 = hom(2)+multi(16)+lowQual(32)+segDup(128)\n";
	print " 181 = het(1)+del(4)+multi(16)+lowQual(32)+segDup(128)\n";
	print " 182 = hom(2)+del(4)+multi(16)+lowQual(32)+segDup(128)\n";
	print " 185 = het(1)+ins(8)+multi(16)+lowQual(32)+segDup(128)\n";
	print " 186 = hom(2)+ins(8)+multi(16)+lowQual(32)+segDup(128)\n";
	print " 257 = het(1)+lowCov(256)\n";
	print " 258 = hom(2)+lowCov(256)\n";
	print " 261 = het(1)+del(4)+lowCov(256)\n";
	print " 262 = hom(2)+del(4)+lowCov(256)\n";
	print " 265 = het(1)+ins(8)+lowCov(256)\n";
	print " 266 = hom(2)+ins(8)+lowCov(256)\n";
	print " 273 = het(1)+multi(16)+lowCov(256)\n";
	print " 274 = hom(2)+multi(16)+lowCov(256)\n";
	print " 277 = het(1)+del(4)+multi(16)+lowCov(256)\n";
	print " 278 = hom(2)+del(4)+multi(16)+lowCov(256)\n";
	print " 281 = het(1)+ins(8)+multi(16)+lowCov(256)\n";
	print " 282 = hom(2)+ins(8)+multi(16)+lowCov(256)\n";
	print " 289 = het(1)+lowQual(32)+lowCov(256)\n";
	print " 290 = hom(2)+lowQual(32)+lowCov(256)\n";
	print " 293 = het(1)+del(4)+lowQual(32)+lowCov(256)\n";
	print " 294 = hom(2)+del(4)+lowQual(32)+lowCov(256)\n";
	print " 297 = het(1)+ins(8)+lowQual(32)+lowCov(256)\n";
	print " 298 = hom(2)+ins(8)+lowQual(32)+lowCov(256)\n";
	print " 305 = het(1)+multi(16)+lowQual(32)+lowCov(256)\n";
	print " 306 = hom(2)+multi(16)+lowQual(32)+lowCov(256)\n";
	print " 309 = het(1)+del(4)+multi(16)+lowQual(32)+lowCov(256)\n";
	print " 310 = hom(2)+del(4)+multi(16)+lowQual(32)+lowCov(256)\n";
	print " 313 = het(1)+ins(8)+multi(16)+lowQual(32)+lowCov(256)\n";
	print " 314 = hom(2)+ins(8)+multi(16)+lowQual(32)+lowCov(256)\n";
	print " 321 = het(1)+segDel(64)+lowCov(256)\n";
	print " 322 = hom(2)+segDel(64)+lowCov(256)\n";
	print " 325 = het(1)+del(4)+segDel(64)+lowCov(256)\n";
	print " 326 = hom(2)+del(4)+segDel(64)+lowCov(256)\n";
	print " 329 = het(1)+ins(8)+segDel(64)+lowCov(256)\n";
	print " 330 = hom(2)+ins(8)+segDel(64)+lowCov(256)\n";
	print " 337 = het(1)+multi(16)+segDel(64)+lowCov(256)\n";
	print " 338 = hom(2)+multi(16)+segDel(64)+lowCov(256)\n";
	print " 341 = het(1)+del(4)+multi(16)+segDel(64)+lowCov(256)\n";
	print " 342 = hom(2)+del(4)+multi(16)+segDel(64)+lowCov(256)\n";
	print " 345 = het(1)+ins(8)+multi(16)+segDel(64)+lowCov(256)\n";
	print " 346 = hom(2)+ins(8)+multi(16)+segDel(64)+lowCov(256)\n";
	print " 353 = het(1)+lowQual(32)+segDel(64)+lowCov(256)\n";
	print " 354 = hom(2)+lowQual(32)+segDel(64)+lowCov(256)\n";
	print " 357 = het(1)+del(4)+lowQual(32)+segDel(64)+lowCov(256)\n";
	print " 358 = hom(2)+del(4)+lowQual(32)+segDel(64)+lowCov(256)\n";
	print " 361 = het(1)+ins(8)+lowQual(32)+segDel(64)+lowCov(256)\n";
	print " 362 = hom(2)+ins(8)+lowQual(32)+segDel(64)+lowCov(256)\n";
	print " 369 = het(1)+multi(16)+lowQual(32)+segDel(64)+lowCov(256)\n";
	print " 370 = hom(2)+multi(16)+lowQual(32)+segDel(64)+lowCov(256)\n";
	print " 373 = het(1)+del(4)+multi(16)+lowQual(32)+segDel(64)+lowCov(256)\n";
	print " 374 = hom(2)+del(4)+multi(16)+lowQual(32)+segDel(64)+lowCov(256)\n";
	print " 377 = het(1)+ins(8)+multi(16)+lowQual(32)+segDel(64)+lowCov(256)\n";
	print " 378 = hom(2)+ins(8)+multi(16)+lowQual(32)+segDel(64)+lowCov(256)\n";
	print " 385 = het(1)+segDup(128)+lowCov(256)\n";
	print " 386 = hom(2)+segDup(128)+lowCov(256)\n";
	print " 389 = het(1)+del(4)+segDup(128)+lowCov(256)\n";
	print " 390 = hom(2)+del(4)+segDup(128)+lowCov(256)\n";
	print " 393 = het(1)+ins(8)+segDup(128)+lowCov(256)\n";
	print " 394 = hom(2)+ins(8)+segDup(128)+lowCov(256)\n";
	print " 401 = het(1)+multi(16)+segDup(128)+lowCov(256)\n";
	print " 402 = hom(2)+multi(16)+segDup(128)+lowCov(256)\n";
	print " 405 = het(1)+del(4)+multi(16)+segDup(128)+lowCov(256)\n";
	print " 406 = hom(2)+del(4)+multi(16)+segDup(128)+lowCov(256)\n";
	print " 409 = het(1)+ins(8)+multi(16)+segDup(128)+lowCov(256)\n";
	print " 410 = hom(2)+ins(8)+multi(16)+segDup(128)+lowCov(256)\n";
	print " 417 = het(1)+lowQual(32)+segDup(128)+lowCov(256)\n";
	print " 418 = hom(2)+lowQual(32)+segDup(128)+lowCov(256)\n";
	print " 421 = het(1)+del(4)+lowQual(32)+segDup(128)+lowCov(256)\n";
	print " 422 = hom(2)+del(4)+lowQual(32)+segDup(128)+lowCov(256)\n";
	print " 425 = het(1)+ins(8)+lowQual(32)+segDup(128)+lowCov(256)\n";
	print " 426 = hom(2)+ins(8)+lowQual(32)+segDup(128)+lowCov(256)\n";
	print " 433 = het(1)+multi(16)+lowQual(32)+segDup(128)+lowCov(256)\n";
	print " 434 = hom(2)+multi(16)+lowQual(32)+segDup(128)+lowCov(256)\n";
	print " 437 = het(1)+del(4)+multi(16)+lowQual(32)+segDup(128)+lowCov(256)\n";
	print " 438 = hom(2)+del(4)+multi(16)+lowQual(32)+segDup(128)+lowCov(256)\n";
	print " 441 = het(1)+ins(8)+multi(16)+lowQual(32)+segDup(128)+lowCov(256)\n";
	print " 442 = hom(2)+ins(8)+multi(16)+lowQual(32)+segDup(128)+lowCov(256)\n";
	print " 513 = het(1)+highCov(512)\n";
	print " 514 = hom(2)+highCov(512)\n";
	print " 517 = het(1)+del(4)+highCov(512)\n";
	print " 518 = hom(2)+del(4)+highCov(512)\n";
	print " 521 = het(1)+ins(8)+highCov(512)\n";
	print " 522 = hom(2)+ins(8)+highCov(512)\n";
	print " 529 = het(1)+multi(16)+highCov(512)\n";
	print " 530 = hom(2)+multi(16)+highCov(512)\n";
	print " 533 = het(1)+del(4)+multi(16)+highCov(512)\n";
	print " 534 = hom(2)+del(4)+multi(16)+highCov(512)\n";
	print " 537 = het(1)+ins(8)+multi(16)+highCov(512)\n";
	print " 538 = hom(2)+ins(8)+multi(16)+highCov(512)\n";
	print " 545 = het(1)+lowQual(32)+highCov(512)\n";
	print " 546 = hom(2)+lowQual(32)+highCov(512)\n";
	print " 549 = het(1)+del(4)+lowQual(32)+highCov(512)\n";
	print " 550 = hom(2)+del(4)+lowQual(32)+highCov(512)\n";
	print " 553 = het(1)+ins(8)+lowQual(32)+highCov(512)\n";
	print " 554 = hom(2)+ins(8)+lowQual(32)+highCov(512)\n";
	print " 561 = het(1)+multi(16)+lowQual(32)+highCov(512)\n";
	print " 562 = hom(2)+multi(16)+lowQual(32)+highCov(512)\n";
	print " 565 = het(1)+del(4)+multi(16)+lowQual(32)+highCov(512)\n";
	print " 566 = hom(2)+del(4)+multi(16)+lowQual(32)+highCov(512)\n";
	print " 569 = het(1)+ins(8)+multi(16)+lowQual(32)+highCov(512)\n";
	print " 570 = hom(2)+ins(8)+multi(16)+lowQual(32)+highCov(512)\n";
	print " 577 = het(1)+segDel(64)+highCov(512)\n";
	print " 578 = hom(2)+segDel(64)+highCov(512)\n";
	print " 581 = het(1)+del(4)+segDel(64)+highCov(512)\n";
	print " 582 = hom(2)+del(4)+segDel(64)+highCov(512)\n";
	print " 585 = het(1)+ins(8)+segDel(64)+highCov(512)\n";
	print " 586 = hom(2)+ins(8)+segDel(64)+highCov(512)\n";
	print " 593 = het(1)+multi(16)+segDel(64)+highCov(512)\n";
	print " 594 = hom(2)+multi(16)+segDel(64)+highCov(512)\n";
	print " 597 = het(1)+del(4)+multi(16)+segDel(64)+highCov(512)\n";
	print " 598 = hom(2)+del(4)+multi(16)+segDel(64)+highCov(512)\n";
	print " 601 = het(1)+ins(8)+multi(16)+segDel(64)+highCov(512)\n";
	print " 602 = hom(2)+ins(8)+multi(16)+segDel(64)+highCov(512)\n";
	print " 609 = het(1)+lowQual(32)+segDel(64)+highCov(512)\n";
	print " 610 = hom(2)+lowQual(32)+segDel(64)+highCov(512)\n";
	print " 613 = het(1)+del(4)+lowQual(32)+segDel(64)+highCov(512)\n";
	print " 614 = hom(2)+del(4)+lowQual(32)+segDel(64)+highCov(512)\n";
	print " 617 = het(1)+ins(8)+lowQual(32)+segDel(64)+highCov(512)\n";
	print " 618 = hom(2)+ins(8)+lowQual(32)+segDel(64)+highCov(512)\n";
	print " 625 = het(1)+multi(16)+lowQual(32)+segDel(64)+highCov(512)\n";
	print " 626 = hom(2)+multi(16)+lowQual(32)+segDel(64)+highCov(512)\n";
	print " 629 = het(1)+del(4)+multi(16)+lowQual(32)+segDel(64)+highCov(512)\n";
	print " 630 = hom(2)+del(4)+multi(16)+lowQual(32)+segDel(64)+highCov(512)\n";
	print " 633 = het(1)+ins(8)+multi(16)+lowQual(32)+segDel(64)+highCov(512)\n";
	print " 634 = hom(2)+ins(8)+multi(16)+lowQual(32)+segDel(64)+highCov(512)\n";
	print " 641 = het(1)+segDup(128)+highCov(512)\n";
	print " 642 = hom(2)+segDup(128)+highCov(512)\n";
	print " 645 = het(1)+del(4)+segDup(128)+highCov(512)\n";
	print " 646 = hom(2)+del(4)+segDup(128)+highCov(512)\n";
	print " 649 = het(1)+ins(8)+segDup(128)+highCov(512)\n";
	print " 650 = hom(2)+ins(8)+segDup(128)+highCov(512)\n";
	print " 657 = het(1)+multi(16)+segDup(128)+highCov(512)\n";
	print " 658 = hom(2)+multi(16)+segDup(128)+highCov(512)\n";
	print " 661 = het(1)+del(4)+multi(16)+segDup(128)+highCov(512)\n";
	print " 662 = hom(2)+del(4)+multi(16)+segDup(128)+highCov(512)\n";
	print " 665 = het(1)+ins(8)+multi(16)+segDup(128)+highCov(512)\n";
	print " 666 = hom(2)+ins(8)+multi(16)+segDup(128)+highCov(512)\n";
	print " 673 = het(1)+lowQual(32)+segDup(128)+highCov(512)\n";
	print " 674 = hom(2)+lowQual(32)+segDup(128)+highCov(512)\n";
	print " 677 = het(1)+del(4)+lowQual(32)+segDup(128)+highCov(512)\n";
	print " 678 = hom(2)+del(4)+lowQual(32)+segDup(128)+highCov(512)\n";
	print " 681 = het(1)+ins(8)+lowQual(32)+segDup(128)+highCov(512)\n";
	print " 682 = hom(2)+ins(8)+lowQual(32)+segDup(128)+highCov(512)\n";
	print " 689 = het(1)+multi(16)+lowQual(32)+segDup(128)+highCov(512)\n";
	print " 690 = hom(2)+multi(16)+lowQual(32)+segDup(128)+highCov(512)\n";
	print " 693 = het(1)+del(4)+multi(16)+lowQual(32)+segDup(128)+highCov(512)\n";
	print " 694 = hom(2)+del(4)+multi(16)+lowQual(32)+segDup(128)+highCov(512)\n";
	print " 697 = het(1)+ins(8)+multi(16)+lowQual(32)+segDup(128)+highCov(512)\n";
	print " 698 = hom(2)+ins(8)+multi(16)+lowQual(32)+segDup(128)+highCov(512)\n";
	print "1025 = het(1)+pathogenic(1024)\n";
	print "1026 = hom(2)+pathogenic(1024)\n";
	print "1029 = het(1)+del(4)+pathogenic(1024)\n";
	print "1030 = hom(2)+del(4)+pathogenic(1024)\n";
	print "1033 = het(1)+ins(8)+pathogenic(1024)\n";
	print "1034 = hom(2)+ins(8)+pathogenic(1024)\n";
	print "1041 = het(1)+multi(16)+pathogenic(1024)\n";
	print "1042 = hom(2)+multi(16)+pathogenic(1024)\n";
	print "1045 = het(1)+del(4)+multi(16)+pathogenic(1024)\n";
	print "1046 = hom(2)+del(4)+multi(16)+pathogenic(1024)\n";
	print "1049 = het(1)+ins(8)+multi(16)+pathogenic(1024)\n";
	print "1050 = hom(2)+ins(8)+multi(16)+pathogenic(1024)\n";
	print "1057 = het(1)+lowQual(32)+pathogenic(1024)\n";
	print "1058 = hom(2)+lowQual(32)+pathogenic(1024)\n";
	print "1061 = het(1)+del(4)+lowQual(32)+pathogenic(1024)\n";
	print "1062 = hom(2)+del(4)+lowQual(32)+pathogenic(1024)\n";
	print "1065 = het(1)+ins(8)+lowQual(32)+pathogenic(1024)\n";
	print "1066 = hom(2)+ins(8)+lowQual(32)+pathogenic(1024)\n";
	print "1073 = het(1)+multi(16)+lowQual(32)+pathogenic(1024)\n";
	print "1074 = hom(2)+multi(16)+lowQual(32)+pathogenic(1024)\n";
	print "1077 = het(1)+del(4)+multi(16)+lowQual(32)+pathogenic(1024)\n";
	print "1078 = hom(2)+del(4)+multi(16)+lowQual(32)+pathogenic(1024)\n";
	print "1081 = het(1)+ins(8)+multi(16)+lowQual(32)+pathogenic(1024)\n";
	print "1082 = hom(2)+ins(8)+multi(16)+lowQual(32)+pathogenic(1024)\n";
	print "1089 = het(1)+segDel(64)+pathogenic(1024)\n";
	print "1090 = hom(2)+segDel(64)+pathogenic(1024)\n";
	print "1093 = het(1)+del(4)+segDel(64)+pathogenic(1024)\n";
	print "1094 = hom(2)+del(4)+segDel(64)+pathogenic(1024)\n";
	print "1097 = het(1)+ins(8)+segDel(64)+pathogenic(1024)\n";
	print "1098 = hom(2)+ins(8)+segDel(64)+pathogenic(1024)\n";
	print "1105 = het(1)+multi(16)+segDel(64)+pathogenic(1024)\n";
	print "1106 = hom(2)+multi(16)+segDel(64)+pathogenic(1024)\n";
	print "1109 = het(1)+del(4)+multi(16)+segDel(64)+pathogenic(1024)\n";
	print "1110 = hom(2)+del(4)+multi(16)+segDel(64)+pathogenic(1024)\n";
	print "1113 = het(1)+ins(8)+multi(16)+segDel(64)+pathogenic(1024)\n";
	print "1114 = hom(2)+ins(8)+multi(16)+segDel(64)+pathogenic(1024)\n";
	print "1121 = het(1)+lowQual(32)+segDel(64)+pathogenic(1024)\n";
	print "1122 = hom(2)+lowQual(32)+segDel(64)+pathogenic(1024)\n";
	print "1125 = het(1)+del(4)+lowQual(32)+segDel(64)+pathogenic(1024)\n";
	print "1126 = hom(2)+del(4)+lowQual(32)+segDel(64)+pathogenic(1024)\n";
	print "1129 = het(1)+ins(8)+lowQual(32)+segDel(64)+pathogenic(1024)\n";
	print "1130 = hom(2)+ins(8)+lowQual(32)+segDel(64)+pathogenic(1024)\n";
	print "1137 = het(1)+multi(16)+lowQual(32)+segDel(64)+pathogenic(1024)\n";
	print "1138 = hom(2)+multi(16)+lowQual(32)+segDel(64)+pathogenic(1024)\n";
	print "1141 = het(1)+del(4)+multi(16)+lowQual(32)+segDel(64)+pathogenic(1024)\n";
	print "1142 = hom(2)+del(4)+multi(16)+lowQual(32)+segDel(64)+pathogenic(1024)\n";
	print "1145 = het(1)+ins(8)+multi(16)+lowQual(32)+segDel(64)+pathogenic(1024)\n";
	print "1146 = hom(2)+ins(8)+multi(16)+lowQual(32)+segDel(64)+pathogenic(1024)\n";
	print "1153 = het(1)+segDup(128)+pathogenic(1024)\n";
	print "1154 = hom(2)+segDup(128)+pathogenic(1024)\n";
	print "1157 = het(1)+del(4)+segDup(128)+pathogenic(1024)\n";
	print "1158 = hom(2)+del(4)+segDup(128)+pathogenic(1024)\n";
	print "1161 = het(1)+ins(8)+segDup(128)+pathogenic(1024)\n";
	print "1162 = hom(2)+ins(8)+segDup(128)+pathogenic(1024)\n";
	print "1169 = het(1)+multi(16)+segDup(128)+pathogenic(1024)\n";
	print "1170 = hom(2)+multi(16)+segDup(128)+pathogenic(1024)\n";
	print "1173 = het(1)+del(4)+multi(16)+segDup(128)+pathogenic(1024)\n";
	print "1174 = hom(2)+del(4)+multi(16)+segDup(128)+pathogenic(1024)\n";
	print "1177 = het(1)+ins(8)+multi(16)+segDup(128)+pathogenic(1024)\n";
	print "1178 = hom(2)+ins(8)+multi(16)+segDup(128)+pathogenic(1024)\n";
	print "1185 = het(1)+lowQual(32)+segDup(128)+pathogenic(1024)\n";
	print "1186 = hom(2)+lowQual(32)+segDup(128)+pathogenic(1024)\n";
	print "1189 = het(1)+del(4)+lowQual(32)+segDup(128)+pathogenic(1024)\n";
	print "1190 = hom(2)+del(4)+lowQual(32)+segDup(128)+pathogenic(1024)\n";
	print "1193 = het(1)+ins(8)+lowQual(32)+segDup(128)+pathogenic(1024)\n";
	print "1194 = hom(2)+ins(8)+lowQual(32)+segDup(128)+pathogenic(1024)\n";
	print "1201 = het(1)+multi(16)+lowQual(32)+segDup(128)+pathogenic(1024)\n";
	print "1202 = hom(2)+multi(16)+lowQual(32)+segDup(128)+pathogenic(1024)\n";
	print "1205 = het(1)+del(4)+multi(16)+lowQual(32)+segDup(128)+pathogenic(1024)\n";
	print "1206 = hom(2)+del(4)+multi(16)+lowQual(32)+segDup(128)+pathogenic(1024)\n";
	print "1209 = het(1)+ins(8)+multi(16)+lowQual(32)+segDup(128)+pathogenic(1024)\n";
	print "1210 = hom(2)+ins(8)+multi(16)+lowQual(32)+segDup(128)+pathogenic(1024)\n";
	print "\n";
	print "Condition of multi allelic:\n";
	print "0/1 het\n";
	print "1/1 hom\n";
	print "0/2 het multi-allelic\n";
	print "1/2 het multi-allelic\n";
	print "2/2 hom multi-allelic\n";
	print "x/y hom multi-allelic (y>1 is multi allelic) (y==x is hom)\n";
	print "\n";
	print "Note:\n";
	print "    If you are using BCF files, please install bcftools\n";
	print "    http://samtools.github.io/bcftools/bcftools.html\n";
	print "\n";
	print "File Suffix:\n";
	print "      .vcf  VCF file\n";
	print "    .g.vcf  GVCF file\n";
	print "  .avinput  AVINPUT file\n";
	print "\n";
	print "Author: akira.hasegawa\@riken.jp\n";
	print "Update: $program_version\n";
	print "\n";
	print "Note:\n";
	print "  If only one VCF file is specified, vcftable.pl assumes VCF file has multiple entries and parse\n";
	print "  Either dup or cov can be calculated in heatmap mode\n";
	print "\n";
}
############################## MAIN ##############################
if($ARGV[0]eq"sortsubs"){sortSubs();exit();}
elsif($ARGV[0]eq"test"){test();exit();}
elsif($ARGV[0]eq"combinations"){printCombinations();exit();}
elsif(defined($opt_h)||scalar(@ARGV)<1){help();exit(1);}
my $lowCovThreshold=0.5;
my $highCovThreshold=2.0;
my $alphaMissenses={};
if(defined($opt_a)){$alphaMissenses=handleAlphaMissense($opt_a);}
if(defined($opt_C)&&!defined($opt_c)){$opt_c="$lowCovThreshold,$highCovThreshold";}
if(defined($opt_c)&&$opt_c=~/,/){
  my @array=split(/,/,$opt_c);
  $lowCovThreshold=$array[0];
  $highCovThreshold=$array[1];
}
my $lowSegThreshold=0.5;
my $highSegThreshold=1.5;
if(defined($opt_S)&&!defined($opt_s)){$opt_s="$lowSegThreshold,$highSegThreshold";}
if(defined($opt_s)&&$opt_s=~/,/){
  my @array=split(/,/,$opt_s);
  $lowSegThreshold=$array[0];
  $highSegThreshold=$array[1];
}
my $threshold=defined($opt_t)?$opt_t:50;
my $tableFile=defined($opt_o)?$opt_o:'-';
my $tmpDir=(-e "/tmp")?"/tmp":"tmp";
mkdir($tmpDir);
my @vcfFiles=listFiles("\\.([vb]cf|avinput)\$",@ARGV);
my ($chromosomes,$basenames,$splitFiles,$covFiles,$refAltFiles)=parseInputs($tmpDir,$threshold,@vcfFiles);
createTable($tableFile,$splitFiles,$chromosomes,$basenames,$alphaMissenses);
if(defined($opt_C)){createTable($opt_C,$covFiles,$chromosomes,$basenames);}
if(defined($opt_S)){createTable($opt_S,$refAltFiles,$chromosomes,$basenames);}
foreach my $chr(keys(%{$splitFiles})){foreach my $splitFile(values(%{$splitFiles->{$chr}})){unlink($splitFile);}}
if($opt_a){foreach my $chr(keys(%{$alphaMissenses})){unlink($alphaMissenses->{$chr});}}
if($tmpDir ne "/tmp"){rmdir($tmpDir);}
############################## absolutePath ##############################
sub absolutePath{
	my $path=shift();
	my $directory=dirname($path);
	my $filename=basename($path);
	return Cwd::abs_path($directory)."/$filename";
}
############################## createFile ##############################
sub createFile{
	my @lines=@_;
	my $path=shift(@lines);
	mkdirs(dirname($path));
	open(OUT,">$path");
	foreach my $line(@lines){print OUT "$line\n";}
	close(OUT);
}
############################## createTable ##############################
sub createTable{
	my $tableFile=shift();
	my $splitFiles=shift();
	my $chromosomes=shift();
	my $basenames=shift();
	my $alphaMissenses=shift();
	if(!defined($opt_q)){print STDERR "#Writing: $tableFile\n";}
	my $writer=IO::File->new(">$tableFile");
	my @labels=();
	print $writer "#chromosome\tposition\tref\talt\t".join("\t",@{$basenames})."\n";
	foreach my $chromosome(@{$chromosomes}){
		my $alphaMissenseReader;
		if(exists($alphaMissenses->{$chromosome})){$alphaMissenseReader=openTable($alphaMissenses->{$chromosome});}
		my @readers=();
		foreach my $basename(@{$basenames}){push(@readers,openTable($splitFiles->{$chromosome}->{$basename}));}
		my $count=0;
		while(nextTable($writer,$alphaMissenseReader,@readers)){$count++;}
		if(!defined($opt_q)){print STDERR "#$chromosome count:\t$count\n";}
	}
	close($writer);
}
############################## equals ##############################
sub equals{
	my $obj1=shift();
	my $obj2=shift();
	my $ref1=ref($obj1);
	my $ref2=ref($obj2);
	if($ref1 ne $ref2){return;}
	if($ref1 eq "ARRAY"){
		my $len1=scalar(@{$obj1});
		my $len2=scalar(@{$obj2});
		if($len1!=$len2){return;}
		for(my $i=0;$i<$len1;$i++){if(!equals($obj1->[$i],$obj2->[$i])){return;}}
		return 1;
	}elsif($ref1 eq "HASH"){
		my @keys1=keys(%{$obj1});
		my @keys2=keys(%{$obj2});
		my $len1=scalar(@keys1);
		my $len2=scalar(@keys2);
		if($len1!=$len2){return;}
		foreach my $key(@keys1){
			if(!exists($obj2->{$key})){return;}
			my $val1=$obj1->{$key};
			my $val2=$obj2->{$key};
			if(!equals($val1,$val2)){return;}
		}
		return 1;
	}
	if($obj1 eq $obj2){return 1;}
}
############################## getAverageDepthOfVCF ##############################
sub getAverageDepthOfVCF{
	my $file=shift();
	my $reader=openFile($file);
	my $sums=[];
	my $count=0;
	while(<$reader>){
			chomp;s/r//g;
			if(/^#/){next;}
			my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@nas)=split(/\t/);
			for(my $i=0;$i<scalar(@nas);$i++){
				my $na=$nas[$i];
				my @formats=split(/:/,$format);
				my @ns=split(/:/,$na);
				for(my $j=0;$j<scalar(@formats);$j++){
					my $f=$formats[$j];
					if($f eq "COV"){
						my $n=$ns[$j];
						my @array=split(/,/,$n);
						my $rc=$array[0];
						my $ac=$array[1];
						$sums->[$i]+=$rc+$ac;
						$count+=1;
					}elsif($f eq "AD"){
						my $n=$ns[$j];
						my @array=split(/,/,$n);
						my $rc=$array[0];
						my $ac=$array[1];
						$sums->[$i]+=$rc+$ac;
						$count+=1;
					}
				}
			}
			if($count==0){last;}
			if($count>=1000000){last;}
	}
	close($reader);
	if($count>0){
		for(my $i=0;$i<scalar(@{$sums});$i++){$sums->[$i]=$sums->[$i]/$count;}
		return $sums;
	}else{
		print STDERR "No depth information\n";exit(1);
	}
}
############################## getDate ##############################
sub getDate{
	my $delim=shift();
	my $time=shift();
	if(!defined($delim)){$delim="";}
	if(!defined($time)||$time eq ""){$time=localtime();}
	else{$time=localtime($time);}
	my $year=$time->year+1900;
	my $month=$time->mon+1;
	if($month<10){$month="0".$month;}
	my $day=$time->mday;
	if($day<10){$day="0".$day;}
	return $year.$delim.$month.$delim.$day;
}
############################## handleAlphaMissense ##############################
sub handleAlphaMissense{
	my $file=shift();
	my $files={};
	my $currentChr;
	my $writer;
	my $reader=openFile($file);
	while(<$reader>){
		chomp;
		if(/^#/){next;}
		my ($chr,$pos,$ref,$alt,$type)=split(/\t/);
		if($currentChr	ne $chr){
			my ($fh,$tmpfile)=tempfile(DIR=>$tmpDir,TEMPLATE=>"alphaMissense.$chr.XXXXXX",SUFFIX=>".txt");
			$files->{$chr}=$tmpfile;
			$writer=$fh;
			$currentChr=$chr;
		}
		if($type ne "likely_pathogenic"){next;}
		print	$writer "$chr\t$pos\t$ref\t$alt\n";
	}
	if(defined($writer)){close($writer);}
	close($reader);
	foreach my $chr(keys(%{$files})){
		my $tmpfile=$files->{$chr};
		my ($fh,$tmpfile2)=tempfile(DIR=>$tmpDir,TEMPLATE=>"alphaMissense.$chr.sort.XXXXXX",SUFFIX=>".txt");
		close($fh);
		if(!defined($opt_q)){print STDERR "#Sorting $chr alphaMissense file: $tmpfile2\n";}
		system("sort -k1,1 -k2,2n -k3,3 -k4,4 $tmpfile>$tmpfile2");
		$files->{$chr}=$tmpfile2;
		unlink($tmpfile);
	}
	return $files;
}
############################## handleNextChromosomes ##############################
sub handleNextChromosomes{
	my $basename=shift();
	my $writers=shift();
	my $files=shift();
	my $preChr=shift();
	my $nextChr=shift();
	if(exists($writers->{$preChr})){close($writers->{$preChr});}
	if(exists($files->{$nextChr}->{$basename})){
		my $tmpfile=$files->{$nextChr}->{$basename};
		$writers->{$nextChr}=IO::File->new(">>$tmpfile");
	}else{
		my ($fh,$tmpfile)=tempfile(DIR=>$tmpDir,TEMPLATE=>"$basename.$nextChr.XXXXXX",SUFFIX=>".txt");
		$writers->{$nextChr}=$fh;
		if(!exists($files->{$nextChr})){$files->{$nextChr}={};}
		$files->{$nextChr}->{$basename}=$tmpfile;
	}
}
############################## listFiles ##############################
sub listFiles{
	my @input_directories=@_;
	my $file_suffix=shift(@input_directories);
	my @input_files=();
	foreach my $input_directory (@input_directories){
		$input_directory=absolutePath($input_directory);
		if(-f $input_directory){push(@input_files,$input_directory);next;}
		elsif(-l $input_directory){push(@input_files,$input_directory);next;}
		opendir(DIR,$input_directory);
		foreach my $file(readdir(DIR)){
			if($file eq "."){next;}
			if($file eq "..") {next;}
			if($file eq ""){next;}
			$file="$input_directory/$file";
			if(-d $file){next;}
			elsif($file!~/$file_suffix$/){next;}
			push(@input_files,$file);
		}
		closedir(DIR);
	}
	return sort(@input_files);
}
############################## mkdirs ##############################
sub mkdirs{
	my @directories=@_;
	foreach my $directory(@directories){
		if(-d $directory){next;}
		my @tokens=split(/[\/\\]/,$directory);
		if(($tokens[0] eq "")&&(scalar(@tokens)>1)){
			shift(@tokens);
			my $token=shift(@tokens);
			unshift(@tokens,"/$token");
		}
		my $string="";
		foreach my $token(@tokens){
			$string.=(($string eq "")?"":"/").$token;
			if(-d $string){next;}
			if(!mkdir($string)){return 0;}
		}
	}
	return 1;
}
############################## nextLine ##############################
sub nextLine{
	my $reader=shift();
	my $in=$reader->[0];
	my $prev=$reader->[1];
	if(eof($in)){
		close($in);
		$reader->[0]=undef;
		$reader->[1]=[];
	}else{
		my $line=<$in>;
		chomp($line);
		$line=~s/\r//g;
		my @tokens=split(/\t/,$line);
		$reader->[1]=\@tokens;
	}
}
############################## nextTable ##############################
sub nextTable{
	my @readers=@_;
	my $writer=shift(@readers);
	my $alphaMissenseReader=shift(@readers);
	my $currentChr;
	my $currentPos;
	my $currentRefAlt;
	foreach my $reader(@readers){
		if(!defined($reader->[0])){next;}
		my ($chr,$pos,$ref,$alt,$flag)=@{$reader->[1]};
		my $refAlt="$ref\t$alt";
		my $update=0;
		if(!defined($currentChr)){$update=1;}
		elsif(($chr cmp $currentChr)<0){$update=1;}
		elsif(($chr cmp $currentChr)>0){$update=0;}
		elsif($pos<$currentPos){$update=1;}
		elsif($pos>$currentPos){$update=0;}
		elsif(($refAlt cmp $currentRefAlt)<0){$update=1;}
		if($update==1){
			$currentChr=$chr;
			$currentPos=$pos;
			$currentRefAlt=$refAlt;
		}
	}
	if(!defined($currentChr)){return 0;}
	my $alphaFlag=0;
	if(defined($alphaMissenseReader)&&defined($alphaMissenseReader->[0])){
		my ($chr,$pos,$ref,$alt)=@{$alphaMissenseReader->[1]};
		my $update=0;
		my $refAlt="$ref\t$alt";
		if(($chr cmp $currentChr)<0){$update=1;}
		elsif(($chr cmp $currentChr)>0){$update=0;}
		elsif($pos<$currentPos){$update=1;}
		elsif($pos>$currentPos){$update=0;}
		elsif(($refAlt cmp $currentRefAlt)==0){$update=1;}
		if($update==1){
			$alphaFlag=1024;
			nextLine($alphaMissenseReader);
		}
	}
	my $line="$currentChr\t$currentPos\t$currentRefAlt";
	my $fileCount=1;
	foreach my $reader(@readers){
		$fileCount++;
		if(!defined($reader->[0])){$line.="\t0";next;}
		my ($chr,$pos,$ref,$alt,$flag)=@{$reader->[1]};
		$flag+=$alphaFlag;
		my $refAlt="$ref\t$alt";
		if($chr ne $currentChr||$pos!=$currentPos||$refAlt ne $currentRefAlt){$line.="\t0";next;}
		$line.="\t$flag";
		nextLine($reader);
	}
	print $writer "$line\n";
	return 1;
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
############################## openTable ##############################
sub openTable{
	my $file=shift();
	my $reader=IO::File->new($file);
	my $line=<$reader>;
	chomp($line);
	$line=~s/\r//g;
	my @tokens=split(/\t/,$line);
	return [$reader,\@tokens];
}
############################## parseAVINPUT ##############################
sub parseAVINPUT{
		my $file=shift();
		my $basename=shift();
		my $threshold=shift();
		my $reader=shift();
		my $writers=shift();
		my $splitFiles=shift();
		my $previousChr;
		while(<$reader>){
		if(/^#/){next;}#skip comment
		chomp;s/\r//g;
		my ($chr,$start,$end,$ref,$alt,$type,$score,$qual,$chr2,$start2,$id,$ref2,$alt2,$score2,$qual2,$info,$format,$na)=split(/\t/);
		if($chr!~/^chr/){$chr="chr$chr";}#add chr
		elsif($chr=~/^CHR(.+)$/i){$chr="chr$1";}#lowercase
		if($chr eq "chrMT"){$chr="chrM";}#handle mitochondria
		my ($a,$b)=split(/[:\|\/]/,$na);#example: 0/1,1/1,1/2
		if($a eq "."){$a=0;}#./2 => 0/2
		if($b eq "."){$b=0;}#1/. => 1/0, just in case
		my $flag=0;
		if($a==0&&$b==0){next;}#wild
		elsif($a==$b&&$a>0){$flag=2;}#homo
		else{$flag=1;}#het
		if($alt eq "-"){$alt=".";}
		elsif($alt eq "0"){$alt="";for(my $i=$start;$i<=$end;$i++){$alt.=".";}}
		if($ref eq "-"){$ref=".";}
		elsif($ref eq "0"){$ref="";for(my $i=$start;$i<=$end;$i++){$ref.=".";}}
		my $multiCount=$alt=~tr/,//;#A,B,C=>3, but haven't seen this happen
		if($multiCount<$b){$multiCount=$b;}#1/2=>2, 3/3=>3, 3/4=>4
		if($multiCount>1){$flag+=16;}#multi allelic
		my $c1=length($alt);# 'TTT'=>3
		my $c2=length($ref);# 'AA'=>2 '.'=>1
		if($c1<$c2){$flag+=4;}#deletion
		elsif($c1>$c2){$flag+=8;}#insertion
		elsif($alt eq "." && $ref ne "."){$flag+=4;}
		elsif($ref eq "." && $alt ne "."){$flag+=8;}
		if($qual eq "PASS"||$qual2 eq "PASS"){}
		elsif($qual=~/^\d+$/){if($qual<$threshold){$flag+=32;}}
		if(!defined($previousChr)||$previousChr ne $chr){#new chromosome
			if(exists($writers->{$previousChr})){close($writers->{$previousChr});}#close
			if(exists($splitFiles->{$chr}->{$basename})){#append to already existing file
				my $tmpfile=$splitFiles->{$chr}->{$basename};
				$writers->{$chr}=IO::File->new(">>$tmpfile");
			}else{#crete new temporary writer
				my ($fh,$tmpfile)=tempfile(DIR=>$tmpDir,TEMPLATE=>"$basename.$chr.XXXXXX",SUFFIX=>".txt");
				$writers->{$chr}=$fh;
				if(!exists($splitFiles->{$chr})){$splitFiles->{$chr}={};}
				$splitFiles->{$chr}->{$basename}=$tmpfile;
			}
			$previousChr=$chr;#update chromosome
		}
		my $writer=$writers->{$chr};#current writer
		print $writer "$chr\t$start\t$ref\t$alt\t$flag\n";#write
	}
	if(exists($writers->{$previousChr})){close($writers->{$previousChr});}#close last writer
	return $basename;
}
############################## parseGVCF ##############################
sub parseGVCF{
  my $file=shift();
  my $basename=shift();
  my $threshold=shift();
  my $reader=shift();
  my $writers=shift();
  my $splitFiles=shift();
  my $previousChr;
  while(<$reader>){
			chomp;s/\r//g;
			if(/^#/){}
			#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Pt2402_BL1811.Bt07
			my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$na)=split(/\t/);
			if($chr!~/^chr/){$chr="chr$chr";}
			elsif($chr=~/^CHR(.+)$/i){$chr="chr$1";}
			if($chr eq "chrMT"){$chr="chrM";}
			my ($a,$b)=split(/[:\|\/]/,$na);
			if($a eq "."){$a=0;}
			if($b eq "."){$b=0;}
			my $flag=0;
			if($a==0&&$b==0){next;}#wild
			elsif($a==$b&&$a>0){$flag=2;}#homo
			else{$flag=1;}#het
			$alt=~s/\,?\<NON_REF\>//g;
			my $multiCount=$alt=~tr/,//;#A,B,C=>3, but haven't seen this happen
			if($multiCount<$b){$multiCount=$b;}#1/2=>2, 3/3=>3, 3/4=>4
			if($multiCount>1){$flag+=16;}#multi allelic
			my $c1=length($alt);
			my $c2=length($ref);# 'AA'=>2 '.'=>1
			if($c1<$c2){$flag+=4;}
			elsif($c1>$c2){$flag+=8;}
			elsif($alt eq "." && $ref ne "."){$flag+=4;}
			elsif($ref eq "." && $alt ne "."){$flag+=8;}
			if($qual eq "PASS"){}
			elsif($qual=~/^\d+$/){if($qual<$threshold){$flag+=32;}}
			if(!defined($previousChr)||$previousChr ne $chr){
				if(exists($writers->{$previousChr})){close($writers->{$previousChr});}
				if(exists($splitFiles->{$chr}->{$basename})){
					my $tmpfile=$splitFiles->{$chr}->{$basename};
					$writers->{$chr}=IO::File->new(">>$tmpfile");
				}else{
					my ($fh,$tmpfile)=tempfile(DIR=>$tmpDir,TEMPLATE=>"$basename.$chr.XXXXXX",SUFFIX=>".txt");
					$writers->{$chr}=$fh;
					if(!exists($splitFiles->{$chr})){$splitFiles->{$chr}={};}
					$splitFiles->{$chr}->{$basename}=$tmpfile;
				}
			}
			$previousChr=$chr;
			my $writer=$writers->{$chr};
			print $writer "$chr\t$pos\t$ref\t$alt\t$flag\n";
  }
  if(exists($writers->{$previousChr})){close($writers->{$previousChr});}
		return $basename;
}
############################## parseInputs ##############################
sub parseInputs{
  my @inFiles=@_;
  my $tmpDir=shift(@inFiles);
  my $threshold=shift(@inFiles);
  my @basenames=();
  my $splitFiles={};
  my $covFiles={};
  my $refAltFiles={};
  my $chrHash={};
  if(!defined($opt_q)){foreach my $inFile(@inFiles){print STDERR "#Input file: $inFile\n";}}
  foreach my $inFile(@inFiles){
    my $reader;
    my $basename;
    if($inFile=~/^(.+)\.g(ip)?z$/i){$basename=$1;$reader=IO::File->new("gzip -cd $inFile|");}
    elsif($inFile=~/^(.+)\.b(ip)?z2$/i){$basename=$1;$reader=IO::File->new("bzip2 -cd $inFile|");}
    elsif($inFile=~/^(.+)\.bcf$/i){$basename=$1;$reader=IO::File->new("bcftools view $inFile|");}
    else{$basename=$inFile;$reader=IO::File->new($inFile);}
    my $type;
    if($basename=~/^(.+)\.g\.vcf$/i){$type="gvcf";$basename=$1;}
    elsif($basename=~/^(.+)\.vcf$/i){$type="vcf";$basename=$1;}
    elsif($basename=~/^(.+)\.bcf$/i){$type="vcf";$basename=$1;}
    elsif($basename=~/^(.+)\.avinput$/i){$type="avinput";$basename=$1;}
    $basename=basename($basename);
    my $writers={};
    my $covWriters={};
    my $refAltWriters={};
				my @handledNames=();
    if($type eq "vcf"){@handledNames=parseVCF($inFile,$basename,$threshold,$reader,$writers,$splitFiles,$covWriters,$covFiles,$refAltWriters,$refAltFiles);}
    elsif($type eq "gvcf"){@handledNames=parseGVCF($inFile,$basename,$threshold,$reader,$writers,$splitFiles);}
    elsif($type eq "avinput"){@handledNames=parseAVINPUT($inFile,$basename,$threshold,$reader,$writers,$splitFiles);}
    close($reader);
				push(@basenames,@handledNames);
				foreach my $handledName(@handledNames){
					foreach my $chr(keys(%{$writers})){
						$chrHash->{$chr}=1;
						sortSplitFile($splitFiles,$chr,$handledName);
						sortSplitFile($covFiles,$chr,$handledName);
						sortSplitFile($refAltFiles,$chr,$handledName);
					}
				}
  }
  my @chromosomes=sort{$a cmp $b}keys(%{$chrHash});
  return (\@chromosomes,\@basenames,$splitFiles,$covFiles,$refAltFiles);
}
############################## parseVCF ##############################
sub parseVCF{
	my $file=shift();
	my $basename=shift();
	my $threshold=shift();
	my $reader=shift();
	my $writers=shift();
	my $splitFiles=shift();
	my $covWriters=shift();
	my $writer;
	my $covWriter;
	my $refAltWriter;
	my $covFiles=shift();
	my $refAltWriters=shift();
	my $refAltFiles=shift();
	my $previousChr;
	my @formats;
	my $averageDepth;
	if(defined($opt_c)){
		$averageDepth=getAverageDepthOfVCF($file);
		if(!defined($opt_q)){
			print STDERR "#Average depth of ".basename($file).":\t".join(",",@{$averageDepth})."\n";
		}
	}
	my @basenames=();
	while(<$reader>){
		if(/^##/){next;}
		if(/^#/){
			chomp;s/\r//g;
			my @tokens=split(/\t/);
			@basenames=splice(@tokens,9);
			last;
		}else{
			close($reader);
			$reader=openFile($file);
			last;
		}
	}
	if(scalar(@basenames)==0){push(@basenames,$basename);}
	while(<$reader>){
		if(/^#/){next;}
		chomp;s/\r//g;
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA18939
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@nas)=split(/\t/);
		if($chr!~/^chr/){$chr="chr$chr";}
		elsif($chr=~/^CHR(.+)$/i){$chr="chr$1";}
		if($chr eq "chrMT"){$chr="chrM";}
		my @flags=();
		my @covs=();
		my @refAlts=();
		my $printFlag=0;
		my @formats=split(/:/,$format);
		my $multiCount=$alt=~tr/,//;#A,B,C=>3, but haven't seen this happen
		my $c1=length($alt);
		my $c2=length($ref);#'AA'=>2 '.'=>1
		my $index=0;
		foreach my $na(@nas){
			my $basename=$basenames[$index];
			my @ns=split(/:/,$na);
			my $flag=0;
			my $covFlag=0;
			my $refAltFlag=0;
			if($c1<$c2){$flag+=4;}
			elsif($c1>$c2){$flag+=8;}
			elsif($alt eq "." && $ref ne "."){$flag+=4;}
			elsif($ref eq "." && $alt ne "."){$flag+=8;}
			if($qual eq "PASS"||$qual eq "."){}
			elsif($qual=~/^\d+$/){if($qual<$threshold){$flag+=32;}}
			for(my $i=0;$i<scalar(@formats);$i++){
				my $f=$formats[$i];
				my $n=$ns[$i];
				if($f eq "GT"){
					my ($a,$b)=split(/[\|\/]/,$n);
					if($a eq "."){$a=0;}
					if($b eq "."){$b=0;}
					if($a==0&&$b==0){next;}#wild
					elsif($a==$b&&$a>0){$flag+=2;}#homo
					else{$flag+=1;}#het
					my $count=$multiCount;
					if($count<$b){$count=$b;}#1/2=>2, 3/3=>3, 3/4=>4
					if($count>1){$flag+=16;}#multi allelic
				}elsif($f eq "COV"||$f eq "AD"){
					my @array=split(/,/,$n);
					my $rc=$array[0];
					my $ac=$array[1];
					if(defined($opt_s)){
						my $ratio=($rc>0)?$ac/$rc:1;
						if($ratio<=$lowSegThreshold){$flag+=64;}
						elsif($ratio>=$highSegThreshold){$flag+=128;}
						if(defined($opt_S)){$refAltFlag=int($ratio*1000+0.5)/1000;}
					}
					if(defined($opt_c)){
							my $ratio=($rc+$ac)/$averageDepth->[0];
							if($ratio<=$lowCovThreshold){$flag+=256;}
							elsif($ratio>=$highCovThreshold){$flag+=512;}
							if(defined($opt_C)){$covFlag=$rc+$ac;}
					}
				}
			}
			if(!defined($previousChr)||$previousChr ne $chr){
				handleNextChromosomes($basename,$writers,$splitFiles,$previousChr,$chr);
				if(defined($opt_S)){handleNextChromosomes($basename,$refAltWriters,$refAltFiles,$previousChr,$chr);}
				if(defined($opt_C)){handleNextChromosomes($basename,$covWriters,$covFiles,$previousChr,$chr);}
				$previousChr=$chr;
			}
			$writer=$writers->{$chr};
			print $writer "$chr\t$pos\t$ref\t$alt\t$flag\n";
			if(defined($opt_C)){
					$covWriter=$covWriters->{$chr};
					print $covWriter "$chr\t$pos\t$ref\t$alt\t$covFlag\n";
			}
			if(defined($opt_S)){
					$refAltWriter=$refAltWriters->{$chr};
					print $refAltWriter "$chr\t$pos\t$ref\t$alt\t$refAltFlag\n";
			}
			$index++;
		}
	}
	close($reader);
	close($writer);
	if(defined($covWriter)){close($covWriter);}
	if(defined($refAltWriter)){close($refAltWriter);}
	return @basenames;
}
############################## printCombinations ##############################
sub printCombinations{
  print "print \"Combinations:\\n\";\n";
  my @labels=(["het","hom"],["del","ins"],["multi"],["lowQual"],["segDel","segDup"],["lowCov","highCov","pathogenic"]);
  my $numbers={"het"=>1,"hom"=>2,"del"=>4,"ins"=>8,"multi"=>16,"lowQual"=>32,"segDel"=>64,"segDup"=>128,"lowCov"=>256,"highCov"=>512,"pathogenic"=>1024};
  my $outputs={0=>""};
  for(my $i=0;$i<scalar(@labels);$i++){
    my $label=$labels[$i];
    my @keys=keys(%{$outputs});
    foreach my $key(@keys){
      my $output=$outputs->{$key};
      for(my $k=0;$k<scalar(@{$label});$k++){
        my $score=$key;
        my $line=$output;
        my $l=$label->[$k];
        my $n=$numbers->{$l};
        if($line eq ""){$line="$l($n)";}
        else{$line.="+$l($n)";}
        $score+=$n;
        $outputs->{$score}=$line;
      }
    }
    if($i==0){delete($outputs->{0});}
  }
  my @keys=sort{$a<=>$b}keys(%{$outputs});
  foreach my $key(@keys){
    my $val=$outputs->{$key};
    for(my $i=length($key);$i<4;$i++){$key=" $key";}
    print "\tprint \"$key = $val\\n\";\n";
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
############################## readText ##############################
sub readText{
	my $file=shift();
	my $text="";
	open(IN,$file);
	while(<IN>){s/\r//g;$text.=$_;}
	close(IN);
	return $text;
}
############################## sortSplitFile ##############################
sub sortSplitFile{
	my $hash=shift();
	my $chr=shift();
	my $name=shift();
	if(!exists($hash->{$chr})){return;}
	if(!exists($hash->{$chr}->{$name})){return;}
	my $tmpfile=$hash->{$chr}->{$name};
	my ($fh,$tmpfile2)=tempfile(DIR=>$tmpDir,TEMPLATE=>"$name.$chr.sort.XXXXXX",SUFFIX=>".txt");
	close($fh);
	if(!defined($opt_q)){print STDERR "#Sorting $chr $name file: $tmpfile2\n";}
	system("sort -k1,1 -k2,2n -k3,3 -k4,4 $tmpfile>$tmpfile2");
	unlink($tmpfile);
	$hash->{$chr}->{$name}=$tmpfile2;
}
############################## sortSubs ##############################
sub sortSubs{
	my $path="$program_directory/$program_name";
	my $reader=openFile($path);
	my @headers=();
	my $name;
	my $blocks={};
	my $block=[];
	my $date=getDate("/");
	my @orders=();
	while(<$reader>){
		chomp;s/\r//g;
		if(/^#{30}\s*(\S+)\s*#{30}$/){
			$name=$1;
			if($name!~/^[A-Z]+$/){push(@{$block},$_);last;}
		}elsif(/^my \$program_version=\"\S+\";/){$_="my \$program_version=\"$date\";";}
		push(@headers,$_);
	}
	while(<$reader>){
		chomp;s/\r//g;
		if(/^#{30}\s*(\S+)\s*#{30}$/){
			$blocks->{$name}=$block;
			push(@orders,$name);
			$name=$1;
			$block=[];
		}
		push(@{$block},$_);
	}
	close($reader);
	if(defined($name)){$blocks->{$name}=$block;push(@orders,$name);}
	my ($writer,$file)=tempfile(DIR=>"/tmp",SUFFIX=>".pl");
	foreach my $line(@headers){print $writer "$line\n";}
	foreach my $key(sort{$a cmp $b}@orders){foreach my $line(@{$blocks->{$key}}){print $writer "$line\n";}}
	close($writer);
	return system("mv $file $path");
}
############################## test ##############################
sub test{
  mkdirs("test");
  createFile(
    "test/input1.avinput",
    "chr1	10440	10440	C	-	hom	24.08	3	chr1	10439	rs112766696	AC	A	24.08	PASS	AC=2;AF=0.077;AN=26;BaseQRankSum=-4.310e-01;ClippingRankSum=0.00;DB;DP=72;ExcessHet=0.1703;FS=0.000;InbreedingCoeff=0.0754;MLEAC=2;MLEAF=0.077;MQ=15.37;MQRankSum=-9.670e-01;QD=8.03;ReadPosRankSum=0.967;SOR=0.223	GT:AD:DP:FT:GQ:PL	1/1:1,2:3:LowDP:4:64,4,0",
    "chr1	15903	15903	-	C	hom	4945.94	3	chr1	15903	rs557514207	G	GC	4945.94	PASS	AC=40;AF=0.870;AN=46;BaseQRankSum=-7.920e-01;ClippingRankSum=0.00;DB;DP=259;ExcessHet=0.0087;FS=24.539;InbreedingCoeff=0.6275;MLEAC=42;MLEAF=0.913;MQ=7.50;MQRankSum=0.00;QD=34.24;ReadPosRankSum=0.674;SOR=3.585	GT:AD:DP:FT:GQ:PL	1/1:0,3:3:LowDP:9:107,9,0",
    "chr1	20317	20317	A	-	het	3958.34	38	chr1	20316	.	GA	G	3958.34	PASS	AC=16;AF=0.348;AN=46;BaseQRankSum=-1.090e-01;ClippingRankSum=0.00;DP=984;ExcessHet=20.9140;FS=0.000;InbreedingCoeff=-0.5274;MLEAC=16;MLEAF=0.348;MQ=7.48;MQRankSum=-9.310e-01;QD=5.54;ReadPosRankSum=-4.720e-01;SOR=0.703	GT:AD:DP:GQ:PL	0/1:23,15:38:99:320,0,587",
    "chr1	81447	81447	A	G	het	492.85	41	chr1	81447	.	A	G	492.85	PASS	AC=1;AF=0.022;AN=46;BaseQRankSum=1.42;ClippingRankSum=0.00;DP=810;ExcessHet=3.0103;FS=0.000;InbreedingCoeff=-0.0223;MLEAC=1;MLEAF=0.022;MQ=44.01;MQRankSum=-1.032e+00;QD=12.02;ReadPosRankSum=0.450;SOR=0.616	GT:AD:DP:FT:GQ:PL	0/1:21,20:41:PASS:99:530,0,571",
    "chr1	82136	82137	AA	-	het	13691.92	22	chr1	82133	.	CAAAA	CAA	13691.92	PASS	AC=11,27,6;AF=0.239,0.587,0.130;AN=46;BaseQRankSum=1.70;ClippingRankSum=0.00;DP=635;ExcessHet=3.1079;FS=1.124;InbreedingCoeff=-0.0403;MLEAC=11,25,6;MLEAF=0.239,0.543,0.130;MQ=10.37;MQRankSum=0.260;QD=26.69;ReadPosRankSum=0.805;SOR=0.878	GT:AD:DP:FT:GQ:PL	2/3:0,0,10,12:22:PASS:99:766,694,676,325,322,261,268,283,0,255",
    "chr1	95011	95011	-	G	het	1320.14	40	chr1	95011	.	T	TG	1320.14	PASS	AC=7,2;AF=0.152,0.043;AN=46;BaseQRankSum=1.10;ClippingRankSum=0.00;DP=777;ExcessHet=7.2151;FS=1.793;InbreedingCoeff=-0.2425;MLEAC=6,2;MLEAF=0.130,0.043;MQ=19.55;MQRankSum=-4.795e+00;QD=4.26;ReadPosRankSum=0.999;SOR=0.866	GT:AD:DP:GQ:PGT:PID:PL	0/1:34,6,0:40:99:0|1:94986_C_T:150,0,2303,252,2321,2573",
    "chr1	104160	104160	-	ACAC	het	11577.69	19	chr1	104160	rs372078516	A	AACAC	11577.69	PASS	AC=26,16;AF=0.565,0.348;AN=46;BaseQRankSum=-1.490e-01;ClippingRankSum=0.00;DB;DP=584;ExcessHet=3.1125;FS=50.214;InbreedingCoeff=0.1776;MLEAC=24,13;MLEAF=0.522,0.283;MQ=7.79;MQRankSum=-1.549e+00;QD=29.00;ReadPosRankSum=0.942;SOR=5.670	GT:AD:DP:FT:GQ:PL	1/2:3,8,7:19:PASS:99:581,166,257,225,0,358",
    "chr1	790136	790136	A	G	het	5687.30	10	chr1	790136	rs6696240	A	G	5687.30	PASS	AC=5,7,2,5;AF=0.109,0.152,0.043,0.109;AN=46;BaseQRankSum=0.269;ClippingRankSum=0.00;DB;DP=590;ExcessHet=51.2979;FS=7.965;InbreedingCoeff=-0.7242;MLEAC=5,7,2,5;MLEAF=0.109,0.152,0.043,0.109;MQ=13.39;MQRankSum=0.00;QD=20.10;ReadPosRankSum=0.693;SOR=0.799	GT:AD:DP:FT:GQ:PGT:PID:PL	0/2:5,0,5,0,0:10:PASS:99:.:.:307,262,605,0,327,310,262,605,327,605,262,605,327,605,605",
    "chr1	791101	791101	T	G	hom	22194.47	23	chr1	791101	rs3131980	T	G	22194.47	PASS	AC=41,5;AF=0.891,0.109;AN=46;BaseQRankSum=1.02;ClippingRankSum=0.00;DB;DP=679;ExcessHet=3.0103;FS=0.000;InbreedingCoeff=0.0000;MLEAC=41,5;MLEAF=0.891,0.109;MQ=12.26;MQRankSum=0.00;QD=34.90;ReadPosRankSum=2.02;SOR=0.526	GT:AD:DP:GQ:PGT:PID:PL	1/1:0,23,0:23:69:.:.:741,69,0,741,69,741",
    "chr1	822428	822498	CCTGGCCAGCAGATCCACCCTGTCTATACTACCTGCCTGGCCAGCAGATCCACCCTGTCTATACTACCTGA	0	hom	18262.34	25	chr1	822428	.	CCTGGCCAGCAGATCCACCCTGTCTATACTACCTGCCTGGCCAGCAGATCCACCCTGTCTATACTACCTGA	0	18262.34	PASS	AC=1,37;AF=0.022,0.804;AN=46;BaseQRankSum=0.00;ClippingRankSum=0.00;DP=838;ExcessHet=6.1884;FS=2.276;InbreedingCoeff=-0.2105;MLEAC=1,37;MLEAF=0.022,0.804;MQ=8.79;MQRankSum=0.00;QD=29.74;ReadPosRankSum=-7.330e-01;SOR=0.870	GT:AD:DP:GQ:PL	2/2:0,0,25:25:78:1022,1022,1022,78,78,0",
    "chr1	840411	840411	A	-	hom	14651.44	22	chr1	840409	rs755461528	TAA	TA	14651.44	PASS	AC=9,8,29;AF=0.196,0.174,0.630;AN=46;BaseQRankSum=-4.920e-01;ClippingRankSum=0.00;DB;DP=775;ExcessHet=3.0103;FS=0.000;InbreedingCoeff=0.0007;MLEAC=8,7,29;MLEAF=0.174,0.152,0.630;MQ=11.74;MQRankSum=-3.780e-01;QD=21.33;ReadPosRankSum=0.406;SOR=0.722	GT:AD:DP:GQ:PL	3/3:0,0,0,22:22:65:463,463,463,463,463,463,65,65,65,0",
    "chr1	984611	984611	-	T	hom	19540.64	12	chr1	984611	.	CTTAT	CTTTAT	19540.64	PASS	AC=3,32,1,2,3,2;AF=0.065,0.696,0.022,0.043,0.065,0.043;AN=46;BaseQRankSum=1.37;ClippingRankSum=0.00;DP=526;ExcessHet=3.3099;FS=1.943;InbreedingCoeff=-0.0688;MLEAC=3,32,1,2,3,2;MLEAF=0.065,0.696,0.022,0.043,0.065,0.043;MQ=12.55;MQRankSum=0.00;QD=26.25;ReadPosRankSum=0.199;SOR=0.647	GT:AD:DP:FT:GQ:PGT:PID:PL	2/2:0,0,12,0,0,0,0:12:PASS:39:1|1:984611_C_CT:584,585,585,39,39,0,585,585,39,585,585,585,39,585,585,585,585,39,585,585,585,585,585,39,585,585,585,585",
    "chr1	1010440	1010440	C	-	hom	24.08	3	chr1	1010439	rs112766696	AC	A	24.08	10	AC=2;AF=0.077;AN=26;BaseQRankSum=-4.310e-01;ClippingRankSum=0.00;DB;DP=72;ExcessHet=0.1703;FS=0.000;InbreedingCoeff=0.0754;MLEAC=2;MLEAF=0.077;MQ=15.37;MQRankSum=-9.670e-01;QD=8.03;ReadPosRankSum=0.967;SOR=0.223	GT:AD:DP:FT:GQ:PL	1/1:1,2:3:LowDP:4:64,4,0",
    "chr1	1015903	1015903	-	C	hom	4945.94	3	chr1	1015903	rs557514207	G	GC	4945.94	11	AC=40;AF=0.870;AN=46;BaseQRankSum=-7.920e-01;ClippingRankSum=0.00;DB;DP=259;ExcessHet=0.0087;FS=24.539;InbreedingCoeff=0.6275;MLEAC=42;MLEAF=0.913;MQ=7.50;MQRankSum=0.00;QD=34.24;ReadPosRankSum=0.674;SOR=3.585	GT:AD:DP:FT:GQ:PL	1/1:0,3:3:LowDP:9:107,9,0",
    "chr1	1020317	1020317	A	-	het	3958.34	38	chr1	1020316	.	GA	G	3958.34	12	AC=16;AF=0.348;AN=46;BaseQRankSum=-1.090e-01;ClippingRankSum=0.00;DP=984;ExcessHet=20.9140;FS=0.000;InbreedingCoeff=-0.5274;MLEAC=16;MLEAF=0.348;MQ=7.48;MQRankSum=-9.310e-01;QD=5.54;ReadPosRankSum=-4.720e-01;SOR=0.703	GT:AD:DP:GQ:PL	0/1:23,15:38:99:320,0,587",
    "chr1	1081447	1081447	A	G	het	492.85	41	chr1	1081447	.	A	G	492.85	13	AC=1;AF=0.022;AN=46;BaseQRankSum=1.42;ClippingRankSum=0.00;DP=810;ExcessHet=3.0103;FS=0.000;InbreedingCoeff=-0.0223;MLEAC=1;MLEAF=0.022;MQ=44.01;MQRankSum=-1.032e+00;QD=12.02;ReadPosRankSum=0.450;SOR=0.616	GT:AD:DP:FT:GQ:PL	0/1:21,20:41:PASS:99:530,0,571",
    "chr1	1082136	1082137	AA	-	het	13691.92	22	chr1	1082133	.	CAAAA	CAA	13691.92	14	AC=11,27,6;AF=0.239,0.587,0.130;AN=46;BaseQRankSum=1.70;ClippingRankSum=0.00;DP=635;ExcessHet=3.1079;FS=1.124;InbreedingCoeff=-0.0403;MLEAC=11,25,6;MLEAF=0.239,0.543,0.130;MQ=10.37;MQRankSum=0.260;QD=26.69;ReadPosRankSum=0.805;SOR=0.878	GT:AD:DP:FT:GQ:PL	2/3:0,0,10,12:22:PASS:99:766,694,676,325,322,261,268,283,0,255",
    "chr1	1095011	1095011	-	G	het	1320.14	40	chr1	1095011	.	T	TG	1320.14	15	AC=7,2;AF=0.152,0.043;AN=46;BaseQRankSum=1.10;ClippingRankSum=0.00;DP=777;ExcessHet=7.2151;FS=1.793;InbreedingCoeff=-0.2425;MLEAC=6,2;MLEAF=0.130,0.043;MQ=19.55;MQRankSum=-4.795e+00;QD=4.26;ReadPosRankSum=0.999;SOR=0.866	GT:AD:DP:GQ:PGT:PID:PL	0/1:34,6,0:40:99:0|1:94986_C_T:150,0,2303,252,2321,2573",
    "chr1	10104160	10104160	-	ACAC	het	11577.69	19	chr1	10104160	rs372078516	A	AACAC	11577.69	16	AC=26,16;AF=0.565,0.348;AN=46;BaseQRankSum=-1.490e-01;ClippingRankSum=0.00;DB;DP=584;ExcessHet=3.1125;FS=50.214;InbreedingCoeff=0.1776;MLEAC=24,13;MLEAF=0.522,0.283;MQ=7.79;MQRankSum=-1.549e+00;QD=29.00;ReadPosRankSum=0.942;SOR=5.670	GT:AD:DP:FT:GQ:PL	1/2:3,8,7:19:PASS:99:581,166,257,225,0,358",
    "chr1	10790136	10790136	A	G	het	5687.30	10	chr1	10790136	rs6696240	A	G	5687.30	17	AC=5,7,2,5;AF=0.109,0.152,0.043,0.109;AN=46;BaseQRankSum=0.269;ClippingRankSum=0.00;DB;DP=590;ExcessHet=51.2979;FS=7.965;InbreedingCoeff=-0.7242;MLEAC=5,7,2,5;MLEAF=0.109,0.152,0.043,0.109;MQ=13.39;MQRankSum=0.00;QD=20.10;ReadPosRankSum=0.693;SOR=0.799	GT:AD:DP:FT:GQ:PGT:PID:PL	0/2:5,0,5,0,0:10:PASS:99:.:.:307,262,605,0,327,310,262,605,327,605,262,605,327,605,605",
    "chr1	10791101	10791101	T	G	hom	22194.47	23	chr1	10791101	rs3131980	T	G	22194.47	18	AC=41,5;AF=0.891,0.109;AN=46;BaseQRankSum=1.02;ClippingRankSum=0.00;DB;DP=679;ExcessHet=3.0103;FS=0.000;InbreedingCoeff=0.0000;MLEAC=41,5;MLEAF=0.891,0.109;MQ=12.26;MQRankSum=0.00;QD=34.90;ReadPosRankSum=2.02;SOR=0.526	GT:AD:DP:GQ:PGT:PID:PL	1/1:0,23,0:23:69:.:.:741,69,0,741,69,741",
    "chr1	10822428	10822498	CCTGGCCAGCAGATCCACCCTGTCTATACTACCTGCCTGGCCAGCAGATCCACCCTGTCTATACTACCTGA	0	hom	18262.34	25	chr1	10822428	.	CCTGGCCAGCAGATCCACCCTGTCTATACTACCTGCCTGGCCAGCAGATCCACCCTGTCTATACTACCTGA	0	18262.34	19	AC=1,37;AF=0.022,0.804;AN=46;BaseQRankSum=0.00;ClippingRankSum=0.00;DP=838;ExcessHet=6.1884;FS=2.276;InbreedingCoeff=-0.2105;MLEAC=1,37;MLEAF=0.022,0.804;MQ=8.79;MQRankSum=0.00;QD=29.74;ReadPosRankSum=-7.330e-01;SOR=0.870	GT:AD:DP:GQ:PL	2/2:0,0,25:25:78:1022,1022,1022,78,78,0",
    "chr1	10840411	10840411	A	-	hom	14651.44	22	chr1	10840409	rs755461528	TAA	TA	14651.44	20	AC=9,8,29;AF=0.196,0.174,0.630;AN=46;BaseQRankSum=-4.920e-01;ClippingRankSum=0.00;DB;DP=775;ExcessHet=3.0103;FS=0.000;InbreedingCoeff=0.0007;MLEAC=8,7,29;MLEAF=0.174,0.152,0.630;MQ=11.74;MQRankSum=-3.780e-01;QD=21.33;ReadPosRankSum=0.406;SOR=0.722	GT:AD:DP:GQ:PL	3/3:0,0,0,22:22:65:463,463,463,463,463,463,65,65,65,0",
    "chr1	10984611	10984611	-	T	hom	19540.64	12	chr1	10984611	.	CTTAT	CTTTAT	19540.64	21	AC=3,32,1,2,3,2;AF=0.065,0.696,0.022,0.043,0.065,0.043;AN=46;BaseQRankSum=1.37;ClippingRankSum=0.00;DP=526;ExcessHet=3.3099;FS=1.943;InbreedingCoeff=-0.0688;MLEAC=3,32,1,2,3,2;MLEAF=0.065,0.696,0.022,0.043,0.065,0.043;MQ=12.55;MQRankSum=0.00;QD=26.25;ReadPosRankSum=0.199;SOR=0.647	GT:AD:DP:FT:GQ:PGT:PID:PL	2/2:0,0,12,0,0,0,0:12:PASS:39:1|1:984611_C_CT:584,585,585,39,39,0,585,585,39,585,585,585,39,585,585,585,585,39,585,585,585,585,585,39,585,585,585,585",
  );
  testCommand(
    "perl $program_directory/vcftable.pl -q test/input1.avinput",
    "#chromosome	position	ref	alt	input1",
    "chr1	10440	C	.	6",
    "chr1	15903	.	C	10",  
    "chr1	20317	A	.	5",
    "chr1	81447	A	G	1",
    "chr1	82136	AA	.	21",
    "chr1	95011	.	G	9",
    "chr1	104160	.	ACAC	25",
    "chr1	790136	A	G	17",
    "chr1	791101	T	G	2",
    "chr1	822428	CCTGGCCAGCAGATCCACCCTGTCTATACTACCTGCCTGGCCAGCAGATCCACCCTGTCTATACTACCTGA	.......................................................................	18",
    "chr1	840411	A	.	22",
    "chr1	984611	.	T	26",
    "chr1	1010440	C	.	38",
    "chr1	1015903	.	C	42",
    "chr1	1020317	A	.	37",
    "chr1	1081447	A	G	33",
    "chr1	1082136	AA	.	53",
    "chr1	1095011	.	G	41",
    "chr1	10104160	.	ACAC	57",
    "chr1	10790136	A	G	49",
    "chr1	10791101	T	G	34",
    "chr1	10822428	CCTGGCCAGCAGATCCACCCTGTCTATACTACCTGCCTGGCCAGCAGATCCACCCTGTCTATACTACCTGA	.......................................................................	50",
    "chr1	10840411	A	.	54",
    "chr1	10984611	.	T	58",
  );
  unlink("test/input1.avinput");
  createFile(
    "test/input1.vcf",
    "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA18939",
    "1	887560	var_452	A	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	0/1:3,9:8.37:-43.31,-8.99,-17.36",
    "1	948921	var_803	T	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=134.03	GT:COV:GT_CONF:GL	1/1:0,13:9.01:-64.60,-16.09,-7.08",
    "1	10887560	var_452	AT	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	0/1:3,9:8.37:-43.31,-8.99,-17.36",
    "1	10948921	var_803	T	CA	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=134.03	GT:COV:GT_CONF:GL	1/1:0,13:9.01:-64.60,-16.09,-7.08",
    "1	20887560	var_452	AT	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	1/3:3,9:8.37:-43.31,-8.99,-17.36",
    "1	20948921	var_803	T	CA	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=134.03	GT:COV:GT_CONF:GL	2/2:0,13:9.01:-64.60,-16.09,-7.08",
  );
  testCommand(
    "perl $program_directory/vcftable.pl -q test/input1.vcf",
    "#chromosome	position	ref	alt	NA18939",
    "chr1	887560	A	C	1",
    "chr1	948921	T	C	2",
    "chr1	10887560	AT	C	5",
    "chr1	10948921	T	CA	10",
    "chr1	20887560	AT	C	21",
    "chr1	20948921	T	CA	26",
  );
  unlink("test/input1.vcf");
		createFile("test/NA18939.vcf",
		"1	887560	var_452	A	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	0/1:3,9:8.37:-43.31,-8.99,-17.36",
		"1	948921	var_803	T	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=134.03	GT:COV:GT_CONF:GL	1/1:0,13:9.01:-64.60,-16.09,-7.08",
		"3	892745	var_489	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=74.92	GT:COV:GT_CONF:GL	0/1:6,3:10.45:-15.80,-5.35,-29.49"
		);
		createFile("test/NA18940.vcf",
		"1	887560	var_452	A	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	0/1:3,9:8.37:-43.31,-8.99,-17.36",
		"2	889238	var_462	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=90.09	GT:COV:GT_CONF:GL	0/1:2,2:8.86:-11.53,-2.67,-11.53",
		"3	892745	var_489	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=74.92	GT:COV:GT_CONF:GL	0/1:6,3:10.45:-15.80,-5.35,-29.49"
		);
		createFile("test/NA18941.vcf",
		"1	1263144	var_2686	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=156.88	GT:COV:GT_CONF:GL	1/1:0,1:0.69:-9.34,-3.87,-3.18",
		"2	889238	var_462	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=90.09	GT:COV:GT_CONF:GL	0/1:2,2:8.86:-11.53,-2.67,-11.53",
		"2	892745	var_489	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=74.92	GT:COV:GT_CONF:GL	0/1:6,3:10.45:-15.80,-5.35,-29.49"
		);
		createFile("test/NA18942.vcf",
		"3	887560	var_452	A	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	0/1:3,9:8.37:-43.31,-8.99,-17.36",
		"3	889238	var_462	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=90.09	GT:COV:GT_CONF:GL	0/1:2,2:8.86:-11.53,-2.67,-11.53",
		"3	892745	var_489	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=74.92	GT:COV:GT_CONF:GL	0/1:6,3:10.45:-15.80,-5.35,-29.49"
		);
  testCommand("perl $program_directory/vcftable.pl -q test/NA18939.vcf test/NA18940.vcf test/NA18941.vcf test/NA18942.vcf","#chromosome	position	ref	alt	NA18939	NA18940	NA18941	NA18942",
	"chr1	887560	A	C	1	1	0	0",
	"chr1	948921	T	C	2	0	0	0",
	"chr1	1263144	G	A	0	0	2	0",
	"chr2	889238	G	A	0	1	1	0",
	"chr2	892745	G	A	0	0	1	0",
	"chr3	887560	A	C	0	0	0	1",
	"chr3	889238	G	A	0	0	0	1",
	"chr3	892745	G	A	1	1	0	1");
 unlink("test/NA18939.vcf");
 unlink("test/NA18940.vcf");
 unlink("test/NA18941.vcf");
 unlink("test/NA18942.vcf");
	createFile("test/NA18943.vcf","#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA18943","1	887560	var_452	A	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	0/1:3,9:8.37:-43.31,-8.99,-17.36","1	887560	var_453	G	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	0/1:3,9:8.37:-43.31,-8.99,-17.36","1	889238	var_462	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=90.09	GT:COV:GT_CONF:GL	0/1:2,2:8.86:-11.53,-2.67,-11.53","1	889238	var_463	T	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=90.09	GT:COV:GT_CONF:GL	0/1:2,2:8.86:-11.53,-2.67,-11.53","1	892745	var_489	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=74.92	GT:COV:GT_CONF:GL	0/1:6,3:10.45:-15.80,-5.35,-29.49","1	892745	var_490	C	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=74.92	GT:COV:GT_CONF:GL	0/1:6,3:10.45:-15.80,-5.35,-29.49","1	909073	var_590	C	T	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=143.24	GT:COV:GT_CONF:GL	0/1:5,4:15.22:-20.17,-4.95,-24.72","1	909073	var_591	A	T	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=143.24	GT:COV:GT_CONF:GL	0/1:5,4:15.22:-20.17,-4.95,-24.72","1	948921	var_803	T	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=134.03	GT:COV:GT_CONF:GL	1/1:0,13:9.01:-64.60,-16.09,-7.08","1	949654	var_804	C	G	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=141.94	GT:COV:GT_CONF:GL	1/1:0,1:0.69:-9.34,-3.87,-3.18","1	949654	var_810	A	G	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=141.94	GT:COV:GT_CONF:GL	1/1:0,1:0.69:-9.34,-3.87,-3.18");
 testCommand("perl $program_directory/vcftable.pl -q test/NA18943.vcf","#chromosome	position	ref	alt	NA18943","chr1	887560	A	C	1","chr1	887560	G	C	1","chr1	889238	G	A	1","chr1	889238	T	A	1","chr1	892745	C	A	1","chr1	892745	G	A	1","chr1	909073	A	T	1","chr1	909073	C	T	1","chr1	948921	T	C	2","chr1	949654	A	G	2","chr1	949654	C	G	2");
	createFile("test/NA18944.vcf","#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA18944","1	887560	var_452	A	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	0/1:3,9:8.37:-43.31,-8.99,-17.36","1	887560	var_453	G	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	0/1:3,9:8.37:-43.31,-8.99,-17.36","1	889238	var_462	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=90.09	GT:COV:GT_CONF:GL	0/1:2,2:8.86:-11.53,-2.67,-11.53","1	889238	var_463	T	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=90.09	GT:COV:GT_CONF:GL	0/1:2,2:8.86:-11.53,-2.67,-11.53","1	892745	var_489	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=74.92	GT:COV:GT_CONF:GL	0/1:6,3:10.45:-15.80,-5.35,-29.49","1	892745	var_490	C	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=74.92	GT:COV:GT_CONF:GL	0/1:6,3:10.45:-15.80,-5.35,-29.49","1	909073	var_590	C	T	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=143.24	GT:COV:GT_CONF:GL	0/1:5,4:15.22:-20.17,-4.95,-24.72","1	909073	var_591	A	T	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=143.24	GT:COV:GT_CONF:GL	0/1:5,4:15.22:-20.17,-4.95,-24.72","1	948921	var_803	T	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=134.03	GT:COV:GT_CONF:GL	1/1:0,13:9.01:-64.60,-16.09,-7.08","1	949654	var_804	C	G	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=141.94	GT:COV:GT_CONF:GL	1/1:0,1:0.69:-9.34,-3.87,-3.18","1	949654	var_810	A	G	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=141.94	GT:COV:GT_CONF:GL	1/1:0,1:0.69:-9.34,-3.87,-3.18");
 testCommand("perl $program_directory/vcftable.pl -q test/NA18943.vcf test/NA18944.vcf","#chromosome	position	ref	alt	NA18943	NA18944","chr1	887560	A	C	1	1","chr1	887560	G	C	1	1","chr1	889238	G	A	1	1","chr1	889238	T	A	1	1","chr1	892745	C	A	1	1","chr1	892745	G	A	1	1","chr1	909073	A	T	1	1","chr1	909073	C	T	1	1","chr1	948921	T	C	2	2","chr1	949654	A	G	2	2","chr1	949654	C	G	2	2");
	createFile("test/NA18945.vcf","#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA18945","1	887560	var_453	T	CG	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=110.65	GT:COV:GT_CONF:GL	0/1:3,9:8.37:-43.31,-8.99,-17.36","1	889238	var_463	T	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=90.09	GT:COV:GT_CONF:GL	0/1:2,2:8.86:-11.53,-2.67,-11.53","1	892745	var_489	G	A	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=74.92	GT:COV:GT_CONF:GL	0/1:6,3:10.45:-15.80,-5.35,-29.49","1	909073	var_590	C	T	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=143.24	GT:COV:GT_CONF:GL	0/1:5,4:15.22:-20.17,-4.95,-24.72","1	948921	var_803	T	C	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=134.03	GT:COV:GT_CONF:GL	1/1:0,13:9.01:-64.60,-16.09,-7.08","1	949654	var_804	C	G	.	PASS	SVLEN=0;SVTYPE=SNP:SITE_CONF=141.94	GT:COV:GT_CONF:GL	1/1:0,1:0.69:-9.34,-3.87,-3.18");
	testCommand("perl $program_directory/vcftable.pl -q test/NA18943.vcf test/NA18944.vcf test/NA18945.vcf","#chromosome	position	ref	alt	NA18943	NA18944	NA18945","chr1	887560	A	C	1	1	0","chr1	887560	G	C	1	1	0","chr1	887560	T	CG	0	0	9","chr1	889238	G	A	1	1	0","chr1	889238	T	A	1	1	1","chr1	892745	C	A	1	1	0","chr1	892745	G	A	1	1	1","chr1	909073	A	T	1	1	0","chr1	909073	C	T	1	1	1","chr1	948921	T	C	2	2	2","chr1	949654	A	G	2	2	0","chr1	949654	C	G	2	2	2",);
	createFile("test/alphaMissense.tsv","#CHROM  POS     REF     ALT     am_class",
	"chr1	887560	A	C	likely_pathogenic",
	"chr1	949654	T	C	likely_pathogenic",
	"chr1	887560	T	CG	likely_pathogenic");
	testCommand("perl $program_directory/vcftable.pl -q -a test/alphaMissense.tsv test/NA18943.vcf test/NA18944.vcf test/NA18945.vcf","#chromosome	position	ref	alt	NA18943	NA18944	NA18945","chr1	887560	A	C	1025	1025	0","chr1	887560	G	C	1	1	0","chr1	887560	T	CG	0	0	1033","chr1	889238	G	A	1	1	0","chr1	889238	T	A	1	1	1","chr1	892745	C	A	1	1	0","chr1	892745	G	A	1	1	1","chr1	909073	A	T	1	1	0","chr1	909073	C	T	1	1	1","chr1	948921	T	C	2	2	2","chr1	949654	A	G	2	2	0","chr1	949654	C	G	2	2	2",);
	unlink("test/NA18943.vcf");
	unlink("test/NA18944.vcf");
	unlink("test/NA18945.vcf");
	unlink("test/alphaMissense.tsv");
	rmdir("test");
}
############################## testCommand ##############################
sub testCommand{
	my @values=@_;
	my $command=shift(@values);
	my $value2=join("\n",@values);
	my ($writer,$file)=tempfile();
	close($writer);
	if(system("$command > $file")){
		print STDERR ">$command\n";
		print STDERR "Command failed...\n";
		return 1;
	}
	my $value1=readText($file);
	chomp($value1);
  $value1=~s/\r//g;
	if($value2 eq""){if($value1 eq""){return 0;}}
	if($value1 eq $value2){return 0;}
	print STDERR ">$command\n";
	print STDERR "$value1\n";
	print STDERR "$value2\n";
}
############################## testSub ##############################
sub testSub{
	my $command=shift();
	my $value2=shift();
	my $value1=eval($command);
	if(equals($value1,$value2)){return 0;}
  print STDERR ">$command\n";
	printTable($value1);
	printTable($value2);
}
