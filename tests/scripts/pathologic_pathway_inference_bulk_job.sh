#!/bin/bash
#SBATCH --account=sortion
#SBATCH --partition=normal
#SBATCH --qos=default
#SBATCH --nodes=1
#SBATCH --time=00:00:30
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10M
#SBATCH --job-name="pathologic_sbatch_launcher"
#SBATCH --comment="Pathologic sbatch bulk launcher"
#SBATCH --output=log/tuto.%x-%j.out
#SBATCH --error=log/tuto.%x-%j.err

set -u
set -e
set -o pipefail

for folder in ./tmp/generated/*/; do
	if [[ -f "${folder}pathologic.pathways.list" ]]; then
		echo "$folder already treated."
	else
		sbatch ./scripts/pathologic_pathway_inference_job.sh "${folder}" "${folder}pathologic.pathways.list"
	fi
done
