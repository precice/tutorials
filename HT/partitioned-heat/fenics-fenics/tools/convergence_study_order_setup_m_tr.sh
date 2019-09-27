python3 create_all_configs.py --monolithic -mth tr -t s -stol 100 -wrl 1 -wrr 1 -Dts 1.0 0.5 0.2 0.1 0.05 0.025 0.0125 -T 1
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py -wrl 1 -wrr 1 -Dts 1.0 0.5 0.2 0.1 0.05 0.025 0.0125
mv runexperiments.sh ..
