#!/bin/bash
#SBATCH --job-name=gridM_setup
#SBATCH --output=gridM_setup_%j.log
#SBATCH --partition=nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00

source /mnt/aridata1/users/arialipt/venvs/turbofit/bin/activate
cd /users/arialipt/turbofit/examples
python -u build_gc_grid.py count M
python -u build_gc_grid.py setup M
python -u build_gc_grid.py time M 3