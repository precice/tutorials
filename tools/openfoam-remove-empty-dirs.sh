#! /bin/sh

openfoam_remove_empty_dirs() {
	(
		set -e -u
		echo "Looking for any time directories without results (e.g. stray functionObjectProperties files, see openfoam-adapter issue #26 on GitHub)..."

		for f in [0-9]* [0-9]*.[0-9]*; do
			if ! [ -f "${f}/U" ] && ! [ -f "${f}/T" ] && ! [ -f "${f}/U.gz" ] && ! [ -f "${f}/T.gz" ]; then
				rm -rfv "${f}"
			fi
		done
		if [ -d processor0 ]; then
			for d in processor*; do
				cd "${d}"
				for f in [0-9]* [0-9]*.[0-9]*; do
					if ! [ -f "${f}/U" ] && ! [ -f "${f}/T" ] && ! [ -f "${f}/U.gz" ] && ! [ -f "${f}/T.gz" ]; then
						rm -rfv "${f}"
					fi
				done
				cd ..
			done
		fi
	)
}
