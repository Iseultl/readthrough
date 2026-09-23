#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --mem 4G
#SBATCH -p genoa64
#SBATCH --qos=pipelines
#SBATCH --mail-type=ALL
#SBATCH --mail-user=iseult.leahy@crg.eu
#SBATCH --output=/no_backup/rg/ileahy/logs/nf_Cyanidium_caldarium_%A.out
#SBATCH --error=/no_backup/rg/ileahy/logs/nf_Cyanidium_caldarium_%A.err
set -e
module load Java
export NXF_JVM_ARGS="-Xms2g -Xmx5g"
cd /users/rg/ileahy/git/gitlab/readthrough
nextflow run main_protists.nf -params-file runs/params_Cyanidium_caldarium.yaml -profile cluster \
    --max_cpus 4 --max_memory 16GB -w /nfs/scratch01/rg/ileahy/nf_work/Cyanidium_caldarium
