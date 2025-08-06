input_file = "probes/0/p"
output_file = "probes.txt"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.strip() and not line.startswith('#'):  # skip empty lines and comments
            parts = line.split()
            outfile.write(';'.join(parts) + '\n')
