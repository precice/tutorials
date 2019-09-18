python3 create_all_configs.py -mth tr -t l -stol 100
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py
mv runexperiments.sh ..
