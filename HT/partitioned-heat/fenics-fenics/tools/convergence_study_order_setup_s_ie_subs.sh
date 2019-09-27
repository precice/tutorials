python3 create_all_configs.py -mth ie -t s -stol 100 -wri linear -exec heat_subcycling.py -T 1 -Dts 1.0 0.5 0.2 0.1 0.05 0.025 0.0125 -wrr 3 -wrl 5
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py -Dts 1.0 0.5 0.2 0.1 0.05 0.025 0.0125 -wrr 3 -wrl 5 
mv runexperiments.sh ..
