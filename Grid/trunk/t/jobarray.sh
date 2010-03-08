#!/bin/bash 
#$ -cwd
#$ -S /bin/bash
#$ -o /tmp/$USER.$JOB_ID.$TASK_ID

# qsub -P 810001 -N SGE_Submit -t 1-3:1 -o /usr/local/devel/ANNOTATION/naxelrod -e /usr/local/devel/ANNOTATION/naxelrod /home/naxelrod/lib/perl/Schedule/examples/jobarray.sh /usr/local/devel/ANNOTATION/naxelrod/examples

# qsub -t 1-2:1 -m abe -M naxelrod@jcvi.org -P 810001 -N SGE_Submit examples/jobarray.sh /home/naxelrod/lib/perl/Schedule/examples /home/naxelrod/lib/perl/Schedule/examples/out
INPUT_DIR=$1
OUTPUT_DIR=$2
TMP="/tmp/$USER.$JOB_ID.$SGE_TASK_ID"

# get an array of fasta files
files=( `ls $INPUT_DIR/data` )

# get the file to work on for this job
task_file=${files[$SGE_TASK_ID - 1]}

# Select which file based on our $SGE_TASK_ID
$INPUT_DIR/hello.pl $task_file > $TMP

# copy the output to our output directory
cp $TMP $OUTPUT_DIR/

