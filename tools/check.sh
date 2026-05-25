#!/usr/bin/env bash
# Run this script at the root of the repository to check images and permalinks

CODE=0

# Check tutorials
IGNORE="tools|quickstart"
tutorials=$(find . -maxdepth 1 -type d -not -name ".*" | grep -vE $IGNORE | sed "s/^.\///")

for tutorial in $tutorials; do
  # Check permalinks
  docs=$(find "./$tutorial" -maxdepth 1 -type f -name "*.md" -print0 | xargs -0 grep -l "permalink:" | sed "s/^.\///")
  for doc in $docs; do
    # Only check files that are tracked by git
    git ls-files --error-unmatch "$doc" >/dev/null 2>&1 || continue

    link=$(grep "permalink:" "$doc" | sed "s/permalink: \+//")
    prefix="tutorials-$tutorial"

    if ! [[ $link =~ ^$prefix ]]; then
      echo "$doc: error: wrong permalink"
      echo "$doc: note: permalink \"$link\" does not start with \"$prefix\""
      echo
      CODE=1
    else
      if [ "${1:-}" = "-v" ]; then
        echo "$doc: info: correct permalink"
        echo "$doc: note: permalink is \"$link\""
        echo
      fi
    fi
  done

  images=$(find "./$tutorial/images" -type f 2> /dev/null | sed "s/^.\///")
  prefix="tutorials-$tutorial-"
  for img in $images; do
    # Only check files that are tracked by git
    git ls-files --error-unmatch "$img" >/dev/null 2>&1 || continue
    if ! [[ $img =~ ^$tutorial/images/$prefix ]]; then
      echo "$img: error: wrong filename"
      echo "$img: note: expected prefix \"$prefix\""
      echo
      CODE=1
    else
      if [ "${1:-}" = "-v" ]; then
        echo "$img: info: correct filename"
        echo
      fi
    fi
  done
done

# Check quickstart
docs=$(find ./quickstart -maxdepth 1 -type f -name "*.md" -print0 | xargs -0 grep -l "permalink:" | sed "s/^.\///")
for doc in $docs; do
  # Only check files that are tracked by git
  git ls-files --error-unmatch "$doc" >/dev/null 2>&1 || continue
  link=$(grep "permalink:" "$doc" | sed "s/permalink: \+//")
  prefix="quickstart"

  if ! [[ $link =~ ^$prefix ]]; then
    echo "$doc: error: wrong permalink"
    echo "$doc: note: permalink \"$link\" does not start with \"$prefix\""
    echo
    CODE=1
  else
    if [ "${1:-}" = "-v" ]; then
      echo "$doc: info: correct permalink"
      echo "$doc: note: permalink is \"$link\""
      echo
    fi
  fi
done

images=$(find ./quickstart/images -type f 2> /dev/null | sed "s/^.\///")
prefix="quickstart-"
for img in $images; do
  # Only check files that are tracked by git
  git ls-files --error-unmatch "$img" >/dev/null 2>&1 || continue
  if ! [[ $img =~ ^quickstart/images/$prefix ]]; then
    echo "$img: error: wrong filename"
    echo "$img: note: expected prefix \"$prefix\""
    echo
    CODE=1
  else
    if [ "${1:-}" = "-v" ]; then
      echo "$img: info: correct filename"
      echo
    fi
  fi
done


[ ! "$CODE" -eq "0" ] && echo "There have been errors"
exit $CODE
