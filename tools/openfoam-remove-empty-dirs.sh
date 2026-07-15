#!/usr/bin/env bash

# Cleaning up stray functionObjectProperties files to not polute the post-processed time steps,
# see https://github.com/precice/openfoam-adapter/issues/26
openfoam_remove_empty_dirs() {
	(
		set -e -u
		echo "Cleaning up any time directories without results to make post-processing easier..."

		DIRECTORIES_TO_CHECK=("." "processor"*)
		OPENFOAM_RESULT_FILES=("p" "U" "T" "D" "pointD" "DD" "pointDD")

		# Search the current and every processor* directory for common results files.
		# If there are none, remove the directory.
		for pd in "${DIRECTORIES_TO_CHECK[@]}"; do
				for d in "${pd}"/[0-9]* "${pd}"/[0-9]*.[0-9]*; do
					KEEP_DIRECTORY=false
					for r in "${OPENFOAM_RESULT_FILES[@]}"; do
						# OpenFOAM can be configured to store files either as compressed or uncompressed
						if [[ -f "${d}/${r}" ]] || [[ -f "${d}/${r}.gz" ]]; then
							KEEP_DIRECTORY=true
							break;
						fi
					done
					if ! [[ ${KEEP_DIRECTORY} == true ]]; then
						# None of the expected result files found - delete the directory
						rm -rfv "${d}"
					fi
				done
		done
		echo "Done."
	)
}
