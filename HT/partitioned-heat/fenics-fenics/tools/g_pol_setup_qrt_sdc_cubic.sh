python3 create_all_configs.py -mth sdc -wri cubic -t qrt -Dts 1.0 0.5 -T 1 --sdc-K 16 -stol 100
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py -wrr 5 3 -wrl 5 3 -Dts 1.0 0.5
mv runexperiments.sh ..
