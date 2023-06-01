# Systemtests
## How to use: 
The intented use of this python script hodgepodge is to ease the execution of systemtests. 

Therefore the workflow for the user is to just execute the `systemtests.py` file. Depending on the options given to the file i will then read in all the metada files and generate the appropriate `docker-compose.yaml` files. 

To test the current state, which only supports openfoam run:
` python systemtests.py --components=openfoam-adapter`

