#!/bin/bash

#SBATCH --job-name=genomesize
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12000 # mb
#SBATCH --time=100:00:00
#SBATCH --output=%x.%j.stdout
#SBATCH --error=%x.%j.stderr

#SBATCH --cpus-per-task=8

##USAGE  sbatch --array=1-14 run_genome_size.sh

module load   GCC/11.3.0  
module load OpenMPI/4.1.4
module load R-bundle-Bioconductor/3.16-R-4.2.2

m=$SLURM_ARRAY_TASK_ID
if [ ! $m ]; then
 m=$1
fi


genome=$(head -n $m genomes.txt | tail -n 1) 

len=5000
path=/data/scratch/projects/punim2009/genome_size/
Rscript ${path}/genome_size1.R  ${genome}.${len}.fq ${genome}.${len}.paf ${path}/${genome}.gsize ${genome}.${len}.fq ${len} false > ${genome}.${len}.out.txt

