#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --job-name=metadata
#SBATCH --output=demux_to_metadata.out
#SBATCH --error=demux_to_metadata.err

module purge
module load mambaforge

# Should have everything needed for the script
conda activate HoMi

python demux_to_metadata.py --infile /pl/active/LozuponeLab/Zimbabwe_shotgun/230703_A00405_0710_AHCY7MDSX7_L1_Lozupone_demux.csv \
--outfile ../HoMi_metadata.csv \
--directory /pl/active/LozuponeLab/Zimbabwe_shotgun/
--skipnrows 5

echo "demux_to_metadata completed successfully"
