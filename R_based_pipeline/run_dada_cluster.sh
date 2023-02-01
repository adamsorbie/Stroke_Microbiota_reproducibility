## SLURM batch job example

# Adjust SBATCH commands below according to your cluster setup

#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH -J "16S_analysis"
#SBATCH --get-user-env
#SBATCH --clusters=cluster_name
#SBATCH --partition=cluster_partition
#SBATCH --cpus-per-task=16
#SBATCH --mem=50gb
#SBATCH --mail-type=end
#SBATCH --mail-user=user@email.com
#SBATCH --export=NONE
#SBATCH --time=08:00:00

# Pipeline parameters
FILEPATH=
OUTPUT=
TRIM_F=
TRIM_R
N_ERROR=
THREADS=
MODE="PE" # change to SE if working with single-end reads

if [ "$MODE" == "PE" ]; then
    Rscript run_dada2_PE.R -p $FILEPATH -o $OUTPUT -f $TRIM_F -r $TRIM_R -n $N_ERROR  -t $THREADS  
elif [ "$MODE" == "SE" ]
    Rscript run_dada2_SE.R -p $FILEPATH -o $OUTPUT -f $TRIM_F -n $N_ERROR  -t $THREADS
else 
   exit 
fi

