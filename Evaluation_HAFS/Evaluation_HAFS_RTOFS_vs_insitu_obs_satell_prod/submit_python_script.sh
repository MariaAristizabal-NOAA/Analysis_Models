#!/bin/bash

#SBATCH --job-name=Glider_vs_HAFS_RTOFS_NESDIS_OHC_stats
#SBATCH --account=hurricane
#SBATCH --qos=batch
##SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --tasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH -t 10:00:00
##SBATCH -t 00:30:00
##SBATCH --partition=xjet
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

source /home/Maria.Aristizabal/load_miniconda_scratch2.sh
conda activate Analysis_env

cd /home/Maria.Aristizabal/Analysis/Evaluation_HAFS/Evaluation_HAFS_RTOFS_vs_satellite_prod
python Glider_vs_HAFS_RTOFS_NESDIS_OHC_stats.py
