#!/usr/bin/env bash

# Cleaning up stray functionObjectProperties files, see https://github.com/precice/openfoam-adapter/issues/26
openfoam_remove_empty_dirs() {
	(
		set -e -u -pipefail
		echo "Cleaning up any time directories without results"

		for d in [0-9]* [0-9]*.[0-9]*; do
                        if ! [ -f "${d}/U" ] && ! [ -f "${d}/T" ] && ! [ -f "${d}/U.gz" ] && ! [ -f "${d}/T.gz" ] && ! [ -f "${d}/D" ] && ! [ -f "${d}/pointD" ] && ! [ -f "${d}/DD" ] && ! [ -f "${d}/pointDD" ] && ! [ -f "${d}/D.gz" ] && ! [ -f "${d}/pointD.gz" ] && ! [ -f "${d}/DD.gz" ] && ! [ -f "${d}/pointDD.gz" ]; then
				rm -rf "${d}"
			fi
		done
		if [[ $(("ls | grep processor | wc -l")) -gt 0 ]]; then
			for pd in processor*; do
				cd "${pd}"
				for d in [0-9]* [0-9]*.[0-9]*; do
                                        if ! [ -f "${d}/U" ] && ! [ -f "${d}/T" ] && ! [ -f "${d}/U.gz" ] && ! [ -f "${d}/T.gz" ] && ! [ -f "${d}/D" ] && ! [ -f "${d}/pointD" ] && ! [ -f "${d}/DD" ] && ! [ -f "${d}/pointDD" ] && ! [ -f "${d}/D.gz" ] && ! [ -f "${d}/pointD.gz" ] && ! [ -f "${d}/DD.gz" ] && ! [ -f "${d}/pointDD.gz" ]; then
						rm -rf "${d}"
					fi
				done
				cd ..
			done
		fi
                echo "Done."
	)
}
