python3 create_all_configs.py -mth tr -pp passive -t s -stol 100 -qntol 0.00001
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py
mv runexperiments.sh ..
