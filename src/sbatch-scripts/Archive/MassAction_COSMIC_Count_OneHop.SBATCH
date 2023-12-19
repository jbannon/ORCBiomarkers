#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --time=72:00:00
#SBATCH --mem=50GB
#SBATCH --job-name=CCCO
#SBATCH --output=./slurmouts/MassAction_COSMIC_Count_OneHop_%j.out

module purge


singularity exec --nv \
	    --overlay /scratch/jjb509/ORCBiomarkers/src/ORC-Overlay.ext3:ro \
	    /scratch/work/public/singularity/cuda11.6.124-cudnn8.4.0.27-devel-ubuntu20.04.4.sif\
	    /bin/bash -c "source /ext3/env.sh; python mass_action_biomarkers.py -config ./configs/COSMIC_COUNT_OH.yaml"
