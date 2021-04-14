#!/bin/sh
# Run this script at the root of the repository to check the size of images

CODE=0
MAXIMUMSIZE=250

# Check tutorials
IGNORE="tools"
tutorials=$(find . -maxdepth 1 -type d -not -name ".*" | grep -vE $IGNORE | sed "s/^.\///")

for tutorial in $tutorials; do
  images=$(find ./"${tutorial}"/images -type f 2> /dev/null | sed "s/^.\///")
  for img in $images; do
    actualsize=$(du -k "$img" | cut -f 1)
    if [ "${actualsize}" -ge "${MAXIMUMSIZE}" ]; then
        echo "Size of ""${img}"" is with ""${actualsize}"" kilobytes larger than the  of ""${MAXIMUMSIZE}"" kilobytes"
        CODE=1
    else
      echo "$img: info: correct file size"
    fi
  done
done

[ ! "$CODE" -eq "0" ] && echo "There have been errors"
exit $CODE
