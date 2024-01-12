#!bash

rm -fr precice-run runs

for solid in solid-{cpp,python,rust}; do
  for fluid in fluid-{cpp,python,rust}; do
    echo "Running $solid - $fluid"
    set -m
    (
      cd $solid && ./clean.sh && ./run.sh &
      cd $fluid && ./clean.sh && ./run.sh &
      wait
    )
    mkdir -p "runs/${solid}-${fluid}"
    cp $fluid/precice-Fluid-watchpoint-Middle.log runs/${solid}-${fluid}/middle.log
    echo "Done ($solid - $fluid)"
  done
done
