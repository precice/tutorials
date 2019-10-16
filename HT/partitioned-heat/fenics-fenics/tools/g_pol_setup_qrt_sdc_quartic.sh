python3 create_all_configs.py -mth sdc -wri quartic -t qrt -Dts 1.0 0.5 0.25 0.1 0.05 0.025 0.0125 -T 1 --sdc-K 16 -stol 100
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py -wrr 5 -wrl 5 -Dts 1.0 0.5 0.25 0.1 0.05 0.025 0.0125
mv runexperiments.sh ..
