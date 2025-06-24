#! /bin/sh

# Cleaning up stray functionObjectProperties files, see https://github.com/precice/openfoam-adapter/issues/26
openfoam_remove_empty_dirs() {
	(
		set -e -u
		echo "Cleaning up any time directories without results"

		for f in [0-9]* [0-9]*.[0-9]*; do
                        if ! [ -f "${f}/U" ] && ! [ -f "${f}/T" ] && ! [ -f "${f}/U.gz" ] && ! [ -f "${f}/T.gz" ] && ! [ -f "${f}/D" ] && ! [ -f "${f}/pointD" ] && ! [ -f "${f}/DD" ] && ! [ -f "${f}/pointDD" ] && ! [ -f "${f}/D.gz" ] && ! [ -f "${f}/pointD.gz" ] && ! [ -f "${f}/DD.gz" ] && ! [ -f "${f}/pointDD.gz" ]; then
				rm -rf "${f}"
			fi
		done
		if [ -d processor0 ]; then
			for d in processor*; do
				cd "${d}"
				for f in [0-9]* [0-9]*.[0-9]*; do
                                        if ! [ -f "${f}/U" ] && ! [ -f "${f}/T" ] && ! [ -f "${f}/U.gz" ] && ! [ -f "${f}/T.gz" ] && ! [ -f "${f}/D" ] && ! [ -f "${f}/pointD" ] && ! [ -f "${f}/DD" ] && ! [ -f "${f}/pointDD" ] && ! [ -f "${f}/D.gz" ] && ! [ -f "${f}/pointD.gz" ] && ! [ -f "${f}/DD.gz" ] && ! [ -f "${f}/pointDD.gz" ]; then
						rm -rf "${f}"
					fi
				done
				cd ..
			done
		fi
                echo "Done."
	)
}
