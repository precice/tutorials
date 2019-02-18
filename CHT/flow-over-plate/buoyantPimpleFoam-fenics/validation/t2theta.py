import csv
import numpy as np

T_inf = 300
dT = 10

def convert(filename):
    reader = csv.DictReader(open(filename))

    result = {}
    for row in reader:
        for key in row:
            print(key)
            try:
                result[key].append(float(row[key]))
            except KeyError:
                result[key] = []
                result[key].append(float(row[key]))

    x = np.array(result['arc_length'])
    y = (np.array(result['T']) - T_inf) / dT

    with open(filename+'.theta', 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',',
			                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['x','theta'])
        for i in range(x.__len__()):
            writer.writerow([x[i], y[i]])

convert('out_FE_OF.csv')
convert('out_OF_OF.csv')
