case_name = "default"
precice_config_file = "default"

# scenario name for log file
scenario_name = "muscle"

# Fixed units in cellMl models:
# These define the unit system.
# 1 cm = 1e-2 m
# 1 ms = 1e-3 s
# 1 uA = 1e-6 A
# 1 uF = 1e-6 F
# 
# derived units:
#   (F=s^4*A^2*m^-2*kg^-1) => 1 ms^4*uA^2*cm^-2*x*kg^-1 = (1e-3)^4 s^4 * (1e-6)^2 A^2 * (1e-2)^-2 m^-2 * (x)^-1 kg^-1 = 1e-12 * 1e-12 * 1e4 F = 1e-20 * x^-1 F := 1e-6 F => x = 1e-14
# 1e-14 kg = 10e-15 kg = 10e-12 g = 10 pg

# (N=kg*m*s^-2) => 1 10pg*cm*ms^2 = 1e-14 kg * 1e-2 m * (1e-3)^-2 s^-2 = 1e-14 * 1e-2 * 1e6 N = 1e-10 N = 10 nN
# (S=kg^-1*m^-2*s^3*A^2, Siemens not Sievert!) => (1e-14*kg)^-1*cm^-2*ms^3*uA^2 = (1e-14)^-1 kg^-1 * (1e-2)^-2 m^-2 * (1e-3)^3 s^3 * (1e-6)^2 A^2 = 1e14 * 1e4 * 1e-9 * 1e-12 S = 1e-3 S = 1 mS
# (V=kg*m^2*s^-3*A^-1) => 1 10pg*cm^2*ms^-3*uA^-1 = (1e-14) kg * (1e-2)^2 m^2 * (1e-3)^-3 s^-3 * (1e-6)^-1 A^-1 = 1e-14 * 1e-4 * 1e6 * 1e6 V = 1e-6 V = 1mV
# (Hz=s^-1) => 1 ms^-1 = (1e-3)^-1 s^-1 = 1e3 Hz
# (kg/m^3) => 1 10 pg/cm^3 = 1e-14 kg / (1e-2 m)^3 = 1e-14 * 1e6 kg/m^3 = 1e-8 kg/m^3
# (Pa=kg/(m*s^2)) => 1e-14 kg / (1e-2 m * 1e-3^2 s^2) = 1e-14 / (1e-8) Pa = 1e-6 Pa

# Hodgkin-Huxley
# t: ms
# STATES[0], Vm: mV
# CONSTANTS[1], Cm: uF*cm^-2
# CONSTANTS[2], I_Stim: uA*cm^-2
# -> all units are consistent

# Shorten
# t: ms
# CONSTANTS[0], Cm: uF*cm^-2
# STATES[0], Vm: mV
# ALGEBRAIC[32], I_Stim: uA*cm^-2
# -> all units are consistent

# Fixed units in mechanics system
# 1 cm = 1e-2 m
# 1 ms = 1e-3 s
# 1 N
# 1 N/cm^2 = (kg*m*s^-2) / (1e-2 m)^2 = 1e4 kg*m^-1*s^-2 = 10 kPa
# (kg = N*s^2*m^-1) => N*ms^2*cm^-1 = N*(1e-3 s)^2 * (1e-2 m)^-1 = 1e-4 N*s^2*m^-1 = 1e-4 kg
# (kg/m^3) => 1 * 1e-4 kg * (1e-2 m)^-3 = 1e2 kg/m^3
# (m/s^2) => 1 cm/ms^2 = 1e-2 m * (1e-3 s)^-2 = 1e4 m*s^-2

# material parameters
# --------------------
# quantities in mechanics unit system
rho = 10                    # [1e-4 kg/cm^3] density of the muscle (density of water)

# Mooney-Rivlin parameters [c1,c2,b,d] of c1*(Ibar1 - 3) + c2*(Ibar2 - 3) + b/d (λ - 1) - b*ln(λ)
# Heidlauf13: [6.352e-10 kPa, 3.627 kPa, 2.756e-5 kPa, 43.373] = [6.352e-11 N/cm^2, 3.627e-1 N/cm^2, 2.756e-6 N/cm^2, 43.373], pmax = 73 kPa = 7.3 N/cm^2
# Heidlauf16: [3.176e-10 N/cm^2, 1.813 N/cm^2, 1.075e-2 N/cm^2, 9.1733], pmax = 7.3 N/cm^2

c1 = 3.176e-10              # [N/cm^2]
c2 = 1.813                  # [N/cm^2]
b  = 1.075e-2               # [N/cm^2] anisotropy parameter
d  = 9.1733                 # [-] anisotropy parameter

material_parameters = [c1, c2, b, d]   # material parameters
pmax = 7.3                  # [N/cm^2] maximum isometric active stress (30-40)
#pmax = 0.73

# load
constant_body_force = (0,0,-9.81e-4)   # [cm/ms^2], gravity constant for the body force
bottom_traction = [0.0,0.0,0.0]        # [N]

# Monodomain parameters
# --------------------
# quantities in CellML unit system
Conductivity = 3.828      # [mS/cm] sigma, conductivity
Am = 500.0                  # [cm^-1] surface area to volume ratio (this is not used, instead values of motor_units are used)
Cm = 0.58                   # [uF/cm^2] membrane capacitance, (1 = fast twitch, 0.58 = slow twitch)
# diffusion prefactor = Conductivity/(Am*Cm)

# timing and activation parameters
# -----------------
# motor units from paper Klotz2019 "Modelling the electrical activity of skeletal muscle tissue using a multi‐domain approach"
import random
random.seed(0)  # ensure that random numbers are the same on every rank
# radius: [μm], stimulation frequency [Hz], jitter [-]
motor_units = [
  {"radius": 40.00, "activation_start_time": 0.0, "stimulation_frequency": 23.92, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},    # low number of fibers
  {"radius": 42.35, "activation_start_time": 0.2, "stimulation_frequency": 23.36, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 45.00, "activation_start_time": 0.4, "stimulation_frequency": 23.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 48.00, "activation_start_time": 0.6, "stimulation_frequency": 22.46, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 51.42, "activation_start_time": 0.8, "stimulation_frequency": 20.28, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 55.38, "activation_start_time": 1.0, "stimulation_frequency": 16.32, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 60.00, "activation_start_time": 1.2, "stimulation_frequency": 12.05, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 65.45, "activation_start_time": 1.4, "stimulation_frequency": 10.03, "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 72.00, "activation_start_time": 1.6, "stimulation_frequency": 8.32,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},
  {"radius": 80.00, "activation_start_time": 1.8, "stimulation_frequency": 7.66,  "jitter": [0.1*random.uniform(-1,1) for i in range(100)]},    # high number of fibers
]

# timing parameters
# -----------------
end_time = 20000.0                      # [ms] end time of the simulation
stimulation_frequency = 100*1e-3    # [ms^-1] sampling frequency of stimuli in firing_times_file, in stimulations per ms, number before 1e-3 factor is in Hertz.
stimulation_frequency_jitter = 0    # [-] jitter in percent of the frequency, added and substracted to the stimulation_frequency after each stimulation
dt_0D = 2e-4                        # [ms] timestep width of ODEs (1e-3)
dt_1D = 2e-4                        # [ms] timestep width of diffusion (1e-3)
dt_splitting = 2e-4                 # [ms] overall timestep width of strang splitting (1e-3)
dt_3D = 1                           # [ms] time step width of coupling, when 3D should be performed, also sampling time of monopolar EMG
output_timestep_fibers = 4e0       # [ms] timestep for fiber output, 0.5
output_timestep_3D = dt_3D              # [ms] timestep for output of fibers and mechanics, should be a multiple of dt_3D


# input files
fiber_file = "../../../../input/left_biceps_brachii_9x9fibers.bin"
#fiber_file = "../../../../input/left_biceps_brachii_31x31fibers.bin"
fat_mesh_file = fiber_file + "_fat.bin"
firing_times_file = "../../../../input/MU_firing_times_always.txt"    # use setSpecificStatesCallEnableBegin and setSpecificStatesCallFrequency
fiber_distribution_file = "../../../../input/MU_fibre_distribution_10MUs.txt"
cellml_file             = "../../../../input/2020_06_03_hodgkin-huxley_shorten_ocallaghan_davidson_soboleva_2007.cellml"

# stride for sampling the 3D elements from the fiber data
# a higher number leads to less 3D elements
sampling_stride_x = 2
sampling_stride_y = 2
sampling_stride_z = 74

# Tolerance value in the element coordinate system of the 3D elements, [0,1]^3
# when a fiber point is still considered part of the element.
# Try to increase this such that all mappings have all points.
mapping_tolerance = 0.5

# other options
paraview_output = True
adios_output = False
exfile_output = False
python_output = False
disable_firing_output = False

# functions, here, Am, Cm and Conductivity are constant for all fibers and MU's
def get_am(fiber_no, mu_no):
  # get radius in cm, 1 μm = 1e-6 m = 1e-4*1e-2 m = 1e-4 cm
  r = motor_units[mu_no]["radius"]*1e-4
  # cylinder surface: A = 2*π*r*l, V = cylinder volume: π*r^2*l, Am = A/V = 2*π*r*l / (π*r^2*l) = 2/r
  return 2./r
  #return Am

def get_cm(fiber_no, mu_no):
  return Cm
  
def get_conductivity(fiber_no, mu_no):
  return Conductivity

def get_specific_states_call_frequency(fiber_no, mu_no):
  stimulation_frequency = motor_units[mu_no % len(motor_units)]["stimulation_frequency"]
  return stimulation_frequency*1e-3

def get_specific_states_frequency_jitter(fiber_no, mu_no):
  #return 0
  return motor_units[mu_no % len(motor_units)]["jitter"]

def get_specific_states_call_enable_begin(fiber_no, mu_no):
  return motor_units[mu_no % len(motor_units)]["activation_start_time"]*1e3