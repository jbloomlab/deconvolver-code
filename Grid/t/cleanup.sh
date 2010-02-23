#!/bin/bash 
#$ -cwd
#$ -S /bin/bash
#$ -o /tmp/$USER.$JOB_ID.$TASK_ID
# qsub -P 810001 -N SGE_Submit -t 1-3:1  /usr/local/devel/ANNOTATION/naxelrod/lib/Grid/t/cleanup.sh 
rm -f /tmp/*naxelrod*

