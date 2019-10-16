python3 create_all_configs.py -mth tr -wri quadratic -t q -T 1 -Dts 1.0 0.5 0.25 0.1 0.05 0.025 0.0125
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py -wrr 2 3 5 -wrl 2 3 5 -Dts 1.0 0.5 0.25 0.1 0.05 0.025 0.0125
mv runexperiments.sh ..
