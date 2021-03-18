#!/bin/bash
# Run this script at the root of the repository to check images for correct prefixes

tutorials=$(find . -mindepth 2 -maxdepth 2 -type d -name images | cut -d/ -f2)

CODE=0
for tutorial in $tutorials; do
  images=$(find ./$tutorial/images -type f)
  prefix="tutorials-$tutorial-"
  for img in $images; do
    if ! [[ $img =~ ^./$tutorial/images/$prefix ]]; then
      echo "$img: error: wrong filename"
      echo "$img: note: expected prefix \"$prefix\""
      CODE=1
    else
      echo "$img: info: correct filename"
    fi
    echo
  done
done

[ ! "$CODE" -eq "0" ] && echo "There have been errors"
exit $CODE
