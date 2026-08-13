#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --mem 8G
#SBATCH -p genoa64
#SBATCH --time 00:10:00
#SBATCH --output=/no_backup/rg/ileahy/logs/nf_Eimeria_necatrix_%A.out   # log per job
#SBATCH --error=/no_backup/rg/ileahy/logs/nf_Eimeria_necatrix_%A.err   # error log
 
# Configure bash
set -e          # exit immediately on error
set -u          # exit immediately if using undefined variables
set -o pipefail # ensure bash pipelines return non-zero status if any of their command fails
 
# Setup trap function to be run when canceling the pipeline job. It will propagate the SIGTERM signal
# to Nextflow so that all jobs launched by the pipeline will be cancelled too.
_term() {
        echo "Caught SIGTERM signal!"
        kill -s SIGTERM $pid
        wait $pid
}
 
trap _term TERM
 
# load Java module
module load Java
 
# limit the RAM that can be used by nextflow
export NXF_JVM_ARGS="-Xms2g -Xmx5g"
 
# Run the pipeline. The command uses the arguments passed to this script, e.g:
#
# $ sbatch submit_nf.sh nextflow/rnatoy -with-singularity
#
# will use "nextflow/rnatoy -with-singularity" as arguments

nextflow run main_protists.nf --params-file params.yaml -profile cluster --max_cpus 4 --max_memory 16GB
 
# Wait for the pipeline to finish
echo "Waiting for ${pid}"
wait $pid
 
# Return 0 exit-status if everything went well
exit 0
