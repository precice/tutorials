import glob
import os
import yaml
import itertools

# Set the directory containing the tutorial folders
tutorial_dir = '../'



# Use globbing to find all the test.yaml files in the directory
yaml_files = glob.glob(f'{tutorial_dir}/*/tests.yaml')

print(f'found {len(yaml_files)} test files to be loaded')
print("\n")
all_combinations=0
# Loop over all the YAML files found and extract the information
for yaml_path in yaml_files:
    # Read the YAML file and extract the information
    with open(yaml_path, 'r') as yaml_file:
        data = yaml.safe_load(yaml_file)

    # Get the tutorial name and its corresponding directory
    tutorial_name = data['name']
    tutorial_path = os.path.join(tutorial_dir, tutorial_name)

    # Get the available solvers for the tutorial
    domains = list(data['couples'].keys())
    solvers = {}
    for domain in domains:
        solvers[domain] = data['couples'][domain]
    
    # Find all possible solver combinations
    solver_combinations = []
    for combination in itertools.product(*[solvers[domain] for domain in domains]):
        solver_combinations.append('::'.join(combination))  

    all_combinations += len(solver_combinations)

    # Print the tutorial name and its available solver combinations
    print(f'Tutorial: {tutorial_name}')
    print(f'Couples the domains: {"::".join(domains)}')
    print(f'Available solver combinations: {", ".join(solver_combinations)}')
    print('-' * 40)

print(f'Found {all_combinations} possible solver combinations in  {len(yaml_files)} test specification files')