python3 create_all_configs.py -mth tr -t s -stol 100 -wri quadratic -qntol 0.00001 -T 1 -Dts 1.0 0.5 0.2 0.1 0.05 0.025 0.0125 0.00625 -wrr 2 5 -wrl 2 5
./config_creation.sh
mv experiments ..
python3 create_runexperiments.py -Dts 1.0 0.5 0.2 0.1 0.05 0.025 0.0125 0.00625 -wrr 2 5 -wrl 2 5 
mv runexperiments.sh ..
