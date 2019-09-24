python3 create_all_configs.py -mth ie -t s -stol 100 -wri quadratic -qntol 0.00001
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py -wrr 3 5 -wrl 3 5 
mv runexperiments.sh ..
