#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --time=72:00:00
#SBATCH --mem=50GB
#SBATCH --job-name=DEEBOH
#SBATCH --output=./slurmouts/Correlation_DiffExp_EB_OneHop_%j.out

module purge


singularity exec --nv \
	    --overlay /scratch/jjb509/ORCBiomarkers/src/ORC-Overlay.ext3:ro \
	    /scratch/work/public/singularity/cuda11.6.124-cudnn8.4.0.27-devel-ubuntu20.04.4.sif\
	    /bin/bash -c "source /ext3/env.sh; python correlation_biomarkers.py -config ./configs/DE_OH.yaml"
