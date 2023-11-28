"""
Class MicroSimulation, which solves a RUC using Abaqus
"""
import os
import subprocess


class MicroSimulation:

    def __init__(self, sim_id):
        """
        Constructor of MicroSimulation class.
        sim_id : int
            global ID of the micro simulation.
        """
        self._dims = 3
        self._sim_id = sim_id
        self._n = 0  # Time counter
        self._first_step = True  # Toggle for first solve step

        # Add 0 before a single digit ID to have a consistent two digit number which can be used in file and folder naming
        # Only works if number of micro simulations is < 100
        if self._sim_id < 10:
            self._id_as_string = '0' + str(self._sim_id)
        else:
            self._id_as_string = str(self._sim_id)

        # File and folder names
        self._foldername = 'ruc_' + self._id_as_string

        # Set the working directory to the micro_ruc_abaqus folder
        os.chdir(os.getcwd() + '/micro_ruc_abaqus')

        # Create a new directory for this micro simulation
        subprocess.call('mkdir ' + self._foldername, shell=True)

        # Copy the input files into the RUC folder
        subprocess.call('cp RUC_*.inp ' + self._foldername, shell=True)

        # Copy the postprocessing script to get stresses from Abaqus
        subprocess.call('cp get_stresses.py ' + self._foldername, shell=True)

        # Copy the cleaning script as we clean the working directory in every time step
        subprocess.call('cp clean_ruc.sh ' + self._foldername, shell=True)

        # Change the working directory to the ruc_ folder
        os.chdir(self._foldername)

        self._jobname = 'RUC_' + self._id_as_string

        log_filename = 'log_ruc_' + self._id_as_string + '_' + str(self._n)

        # Run the initial Abaqus simulation
        subprocess.call('abaqus job=' + self._jobname + ' input=RUC_initial \
                  scratch=' + os.getcwd() + ' interactive double=both \
                  &> ' + log_filename + '.log', shell=True)

    def solve(self, strains, dt):
        assert dt != 0

        # Set the working directory to the micro_ruc/ folder
        os.chdir(os.getcwd() + self._foldername)

        # Clean the working directory
        subprocess.call('sh clean_ruc.sh', shell=True)

        if self._first_step:
            subprocess.call('mv RUC_initial.inp RUC_nm1.inp', shell=True)
            self._first_step = False
        else:
            subprocess.call('mv RUC_iterate.inp RUC_nm1.inp', shell=True)

        # Open the current input file and read all the lines
        input_file = open('RUC_nm1.inp', "r")
        line_list = input_file.readlines()
        input_file.close()

        # Create a new input file to write the modified lines in to
        new_file = open('RUC_iterate.inp', "w")

        # Modify the lines with the strain values
        line_list[4] = 'deps11   = {}\n'.format(strains["strains1to3"][0])
        line_list[5] = 'deps22   = {}\n'.format(strains["strains1to3"][1])
        line_list[6] = 'deps33   = {}\n'.format(strains["strains1to3"][2])
        line_list[7] = 'dgamma12 = {}\n'.format(strains["strains4to6"][0])
        line_list[8] = 'dgamma23 = {}\n'.format(strains["strains4to6"][1])
        line_list[9] = 'dgamma13 = {}\n'.format(strains["strains4to6"][2])

        # Write the modified lines to the input file and close it
        new_file.writelines(line_list)
        new_file.close()

        self._n += 1

        # Run the Abaqus simulation
        subprocess.call('abaqus job=' + self._jobname + ' input=RUC_iterate \
                  scratch=' + os.getcwd() + ' interactive double=both \
                  &> log_ruc_run.log', shell=True)

        # Make sure that .odb file has been created
        assert os.path.exists(os.getcwd() + '/RUC_' + self._id_as_string + '.odb')

        # Get the stresses
        subprocess.call('abaqus cae noGUI=get_stresses.py', shell=True)

        # Open output file and read stress values
        stresses_file = open('stresses.txt', 'r')
        stresses_as_strings = stresses_file.readlines()
        stresses_file.close()

        stresses = []
        for stress_as_string in stresses_as_strings:
            stresses.append(float(stress_as_string))

        return {"stresses1to3": stresses[0:3], "stresses4to6": stresses[3:6]}
