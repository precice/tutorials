python3 create_all_configs.py -dd DN -g 1.0 
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py
mv runexperiments.sh ..
