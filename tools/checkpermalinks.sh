#!/bin/bash
# Run this script at the root of the repository to check permalinks

FILES=$(find . -mindepth 2 -maxdepth 2 -type f -name "*.md" -not -path "./quickstart/*")
tutorials=$(echo "$FILES" | cut -d/ -f2 | sort -u)

CODE=0
for tutorial in $tutorials; do
  files=$(find ./$tutorial -maxdepth 1 -type f -name "*.md")
  for file in $files; do
    link=$(grep "permalink:" $file | sed "s/permalink: \+//")
    prefix="tutorials-$tutorial"

    if ! [[ $link =~ ^$prefix ]]; then
      echo "$file: error: wrong permalink"
      echo "$file: note: permalink \"$link\" does not start with \"$prefix\""
      CODE=1
    else
      echo "$file: info: correct permalink"
      echo "$file: note: permalink is \"$link\""
    fi
    echo
  done
done

[ ! "$CODE" -eq "0" ] && echo "There have been errors"
exit $CODE
