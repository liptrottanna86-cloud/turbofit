#!/bin/bash
#SBATCH --job-name=synth_c13
#SBATCH --output=synth_c13_%j.log
#SBATCH --partition=nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00

source /mnt/aridata1/users/arialipt/venvs/turbofit/bin/activate
cd /users/arialipt/turbofit/examples
python -u synth_c12c13.py