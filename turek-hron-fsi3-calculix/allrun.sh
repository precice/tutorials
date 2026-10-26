echo "START FLUID"

cd fluid-openfoam/
./run.sh -parallel > log-fluid &
cd ..

echo "START SOLID"

cd solid-calculix/
./run.sh > log-solid &
cd ..

wait


echo "FINISHED"
