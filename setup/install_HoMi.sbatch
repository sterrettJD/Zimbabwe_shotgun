#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --partition=amilan
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --job-name=install_HoMi
#SBATCH --output=install_HoMi.out
#SBATCH --error=install_HoMi.err

module purge
module load mambaforge

conda create -n HoMi python
conda activate HoMi
pip install -e ../../HoMi/

echo "HoMi has been successfully installed!"
