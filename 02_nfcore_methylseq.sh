#!/bin/bash

##  Set username
USER=glamka

## Set dir name
DIR=kratepi_02_nfcore_methylseq

## Set project names
PROJ1=bismark  #bwameth

## Set error location
#ERR = work/


## --------------------------------
## Activate conda environment and load modules

module load python/anaconda/3.10.9

source activate nextflow
module load java/15.0.1
module load singularity/3.8.4


## --------------------------------
## Download test data

mkdir /scratch/${USER}/${DIR}/
chmod -R 777 /scratch/${USER}/${DIR}/
mkdir /scratch/${USER}/${DIR}/${PROJ1}/
cd /scratch/${USER}/${DIR}/${PROJ1}/

mkdir data
cp -r ~/kratepi/ref_fasta .
cp ~/kratepi/kratepi_sample_names_byhand.csv ./data/
cp ~/kratepi/easleyconfig_gfl.conf .

## --------------------------------
## Run nf-score methylseq pipeline

mkdir results


## aligning with bismark -- requires a manual fix after this step
nextflow run nf-core/methylseq \
--input data/kratepi_sample_names_byhand.csv \
--outdir results \
--fasta ref_fasta/dspec_genbank_assem.fa.gz \
--save_reference \
-resume \
-bg \
-c easleyconfig_gfl.conf \
-profile singularity 



## for some reason, ref .fa can't be used when it's just a symbolic link? idk. this fixes
## it tho, alongside specif.
#cd /scratch/${USER}/${DIR}/${PROJ1}/${ERR}
#rm dspec_genbank_assem.fa
#cp /scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa

#cd /scratch/glamka/ms_kratepi_02_nfcore_methylseq/bismark/work/d5/082d525c5435b5851311157f54b344/BismarkIndex # /scratch/glamka/ms_kratepi_02_nfcore_methylseq_rerun/bismark/work/04/f271e46329a2d092e200bc2c09d84c/BismarkIndex
#rm dspec_genbank_assem.fa.gz
#cp /home/gfl0003/kratepi/ref_fasta/dspec_genbank_assem.fa /scratch/glamka/ms_kratepi_02_nfcore_methylseq/bismark/work/d5/082d525c5435b5851311157f54b344/BismarkIndex/dspec_genbank_assem.fa  # /scratch/glamka/ms_kratepi_02_nfcore_methylseq_rerun/bismark/work/04/f271e46329a2d092e200bc2c09d84c/BismarkIndex/dspec_genbank_assem.fa.gz 
###cd /scratch/avrilh/ms_kratroh_02_nfcore_methylseq/bismark/work/fc/e4a108a5ce95a58106631b57f1dd5b/BismarkIndex/
###rm dspec_genbank_assem.fa
###cp /scratch/glamka/kratroh_02_nfcore_methylseq/bismark_TEST/ref_fasta/dspec_genbank_assem.fa .             #/scratch/avrilh/kratroh_01_assembindex/dspec_genbank_assem.fa

# cd /scratch/${USER}/${DIR}/${PROJ1}/ 

#nextflow run nf-core/methylseq \
#--input data/kratepi_sample_names_byhand.csv \
#--outdir results \
#--fasta ref_fasta/dspec_genbank_assem.fa \
#--save_reference \
#-resume \
#-bg \
#-c easleyconfig_gfl.conf \
#-profile singularity \
#--aligner bwameth 

