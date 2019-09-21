python3 create_all_configs.py -mth ie -t s -stol 100 -ctol 0.00001 -pp none
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py
mv runexperiments.sh ..
