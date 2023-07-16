# This script declares to SCons how to compile the example.
# It has to be called from a SConstruct file.
# The 'env' object is passed from there and contains further specification like directory and debug/release flags.
#
# Note: If you're creating a new example and copied this file, adjust the desired name of the executable in the 'target' parameter of env.Program.


Import('env')     # import Environment object from calling SConstruct

# if the option no_tests was given, quit the script
if not env['no_examples']:
    
  # create the main executable
  env.Program(target = 'muscle_contraction', source = "src/muscle_contraction.cpp")
