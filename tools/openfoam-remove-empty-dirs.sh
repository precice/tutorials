#! /bin/sh
set -e -u

openfoam_remove_empty_dirs() {
	echo "Looking for any time directories without results (e.g. stray functionObjectProperties files, see openfoam-adapter issue #26 on GitHub)..."

	for f in [0-9]* [0-9]*.[0-9]*; do
		if ! [ -f "${f}/U" ] && ! [ -f "${f}/T" ]; then
			rm -rfv "${f}"
		fi
	done
	if [ -d processor0 ]; then
		for d in processor*; do
			cd "${d}"
			for f in [0-9]* [0-9]*.[0-9]*; do
				if ! [ -f "${f}"/U ] && ! [ -f "${f}"/T ]; then
					rm -rfv "${f}"
				fi
			done
			cd ..
		done
	fi
}
