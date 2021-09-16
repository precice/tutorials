#!/bin/bash
# Run this script at the root of the repository to check the size of images

CODE=0
MAXIMUMSIZE=250
MAXIMUMGIFSIZE=2200

RED='\033[0;31m'
NOCOLOR='\033[0m'

# Check tutorials
IGNORE="tools"
tutorials=$(find . -maxdepth 1 -type d -not -name ".*" | grep -vE $IGNORE | sed "s/^.\///")

echo "Limit for regular images: ${MAXIMUMSIZE} kb"
echo "Limit for gifs: ${MAXIMUMGIFSIZE} kb"
# For all tutorials do
for tutorial in $tutorials; do
  images=$(find ./"${tutorial}"/images -type f 2> /dev/null | sed "s/^.\///")
  for img in $images; do
    actualsize=$(du -k "$img" | cut -f 1)
    # Check gifs
    if [[ "${img}" == *.gif ]]; then
      if [ "${actualsize}" -ge "${MAXIMUMGIFSIZE}" ]; then
          echo -e "$img:$RED $actualsize kb exceeds the limit of $MAXIMUMGIFSIZE kb. $NOCOLOR"
          CODE=1
      else
          echo -e "$img: $actualsize kb (Ok)."
      fi
    else
      if [ "${actualsize}" -ge "${MAXIMUMSIZE}" ]; then
         echo -e "$img:$RED $actualsize kb exceeds the limit of $MAXIMUMSIZE kb. $NOCOLOR"
         CODE=1
      else
        echo -e "$img: $actualsize kb (Ok)."
      fi
    fi
  done
done

[ ! "$CODE" -eq "0" ] && echo "There have been errors"
exit $CODE
