#!/bin/bash
#SBATCH --account=sortion
#SBATCH --partition=normal
#SBATCH --qos=default
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name="pathologic_sbatch"
#SBATCH --comment="Pathologic sbatch"
#SBATCH --output=log/tuto.%x-%j.out
#SBATCH --error=log/tuto.%x-%j.err

set -u
set -e
set -o pipefail

printf "Script de soumission executé sur %s avec %d CPUs\n" $(hostname) $(nproc)

if [[ $# -eq 2 ]]; then
	input="${1}"
	output="${2}"
	echo "reading $input, writing to $output"
else
	echo "Error: check the command line arguments" > /dev/stderr
	echo "Usage: $0 <pathologic_input> <pathway_list>" > /dev/stderr
	exit 1
fi

module purge
module load extenv/ibfj
module load extenv/labgem

module load pathway-tools
module load pythoncyc

TMPDIR="$(mktemp -d "/tmp/tmp.$USER.ptools_local.XXXXXXXXXX")"
create_ptools_local "${TMPDIR}"
echo "Temporary ptools-local created in ${TMPDIR}"
eval $(use_ptools_local "${TMPDIR}")

pathologic_batch "$input" & PATHOLOGIC_PID=$!

wait $PATHOLOGIC_PID

org_id="$(basename "${input}")"
# pathway_list_folder="./tmp/generated.pathologic.out/${org_id}"
# mkdir -p "${pathway_list_folder}"
# pathway_list_file="${pathway_list_folder}/pathway.list"
bash scripts/export-pathways.sh "${org_id}" "${output}"
