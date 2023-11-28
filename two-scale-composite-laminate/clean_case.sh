rm -fv slurm*
rm -rfv precice-profiling
rm -rfv precice-run

echo "Cleaning meso Abaqus participant"
cd meso_laminate_abaqus
./clean.sh

echo "Cleaning micro Abaqus participant"
cd ../micro_ruc_abaqus
./clean.sh

echo "Cleaning micro NASMAT participant"
cd ../micro_ruc_nasmat
./clean.sh
