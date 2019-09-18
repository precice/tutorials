python3 create_all_configs.py -mth tr -t s -exec heat_subcycling.py -stol 100 -qntol 0.00001
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py
mv runexperiments.sh ..
