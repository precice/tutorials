python3 create_all_configs.py -mth ie -t l -stol 100 -pp passive
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py
mv runexperiments.sh ..
