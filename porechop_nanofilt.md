'''shell

#!/bin/bash
### SLURM specific parameters, for e.g:

#SBATCH --output=/home/haicu/ttreska57/minimap2.log
#SBATCH -e /home/haicu/ttreska57/minimap2.log
#SBATCH -J run_minimap2nt
#SBATCH -p gpu_p
#SBATCH -q gpu
#SBATCH --mem=32G
#SBATCH -t 12:00:00
#SBATCH --qos gpu
#SBATCH --partition=gpu_p
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=ALL
#SBATCH --gres=gpu:1

### User's commands, apps, parameters, etc. for e.g:

### First Task:
#!/bin/bash
for file in *.fastq
do
    base=$(basename $file .fastq)
    porechop -i $file -o ${base}_trimmed.fastq
    NanoFilt -l 100 -q 7 ${base}_trimmed.fastq > ${base}_filtered.fastq
done
