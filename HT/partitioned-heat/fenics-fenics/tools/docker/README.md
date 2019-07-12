# Use docker to run everything

* Build precice baseimage: `docker build -f Dockerfile.Ubuntu1804.home -t experiments_base .`
* Build fenics-adapter image on top: `docker build -f Dockerfile.ExperimentRunner --build-arg from=experiments_base:latest -t experiments_fenics .`
* Start an interactive session to run experiments: `docker run -it -v $(pwd)/dockerout:/tutorials/HT/partitioned-heat/fenics-fenics/experiments experiments_fenics:latest`. Output from the experiments inside the docker container can be written to `/tutorials/HT/partitioned-heat/fenics-fenics/experiments`. This is mapped to `$(pwd)/dockerout` on the hostmachine. **DON'T DELETE THE EXPERIMENTS FOLDER ON THE CONTAINER!**
