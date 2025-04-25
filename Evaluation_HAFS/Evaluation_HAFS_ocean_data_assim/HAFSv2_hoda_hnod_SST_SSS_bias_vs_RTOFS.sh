#!/bin/sh
#SBATCH --job-name=hoda_hnod_SST_SSS_bias
#SBATCH --account=hurricane
#SBATCH --qos=batch
##SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --tasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH -t 06:00:00
##SBATCH -t 00:30:00
##SBATCH --partition=xjet
##SBATCH --partition=orion
#SBATCH -o hoda_hnod_SST_SSS_bias.log.%j
#SBATCH -e hoda_hnod_SST_SSS_bias.log.%j
##SBATCH --mem=0
##SBATCH --exclusive
#SBATCH -D.

source /home/Maria.Aristizabal/load_miniconda.sh
python /home/Maria.Aristizabal/Analysis/Evaluation_HAFS/Evaluation_HAFS_ocean_data_assim/HAFSv2_hoda_hnod_SST_SSS_bias_vs_RTOFS.py
