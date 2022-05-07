#!/bin/csh

#SBATCH --nodes=1
#SBATCH --time=02:25:00
#SBATCH --qos=regular
#SBATCH --constraint=knl
#SBATCH --account=m3525
#SBATCH --output=a.out
#SBATCH --error=a.err

#module load python
#conda activate esmac_diags

#python prep_HISCALE_allobs.py
#python prep_ACEENA_allobs.py
#python prep_MAGIC_allobs.py
#python prep_MARCUS_allobs.py
#python prep_CSET_allobs.py
#python prep_SOCRATES_allobs.py

python prep_HISCALE_E3SM.py

exit
