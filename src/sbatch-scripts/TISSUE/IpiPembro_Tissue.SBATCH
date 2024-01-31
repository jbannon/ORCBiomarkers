#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --time=96:00:00
#SBATCH --mem=50GB
#SBATCH --job-name=IpiPembro_Tissue
#SBATCH --output=./slurmouts/IpiPembro_Tissue_%j.out
#SBATCH --mail-type=TIME_LIMIT,END,FAIL
#SBATCH --mail-user=jjb509@nyu.edu

module purge


singularity exec --nv \
	    --overlay /scratch/jjb509/ORCBiomarkers/src/ORC-Overlay.ext3:ro \
	    /scratch/work/public/singularity/cuda11.6.124-cudnn8.4.0.27-devel-ubuntu20.04.4.sif\
	    /bin/bash -c "source /ext3/env.sh; python tissue_wide_biomarkers.py -config ./configs/CRI/IpiPembro.yaml"
