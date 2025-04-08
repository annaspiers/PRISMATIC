#!/bin/bash

#SBATCH --qos=regular
#SBATCH --job-name=main_py        # Job name
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --constraint=cpu          # architecture
#SBATCH --time=0:30:00            # Time limit (hh:mm:ss)
#SBATCH --output=output_%j.log    # Standard output and error log
#SBATCH --account=m4460

. $HOME/.bashrc
conda activate prismatic
srun python $HOME/dev/PRISMATIC/main.py