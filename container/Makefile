build:
	apptainer build pan2met.sif pan2met.def

test:
	apptainer run --bind ${WORK},${TMPDIR} pan2met.sif python -m pan2met --help
