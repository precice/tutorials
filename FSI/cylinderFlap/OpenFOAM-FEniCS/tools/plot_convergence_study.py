import csv
import numpy as np
import matplotlib.pyplot as plt


def outfile_to_data(filename):
    csv_file = open(filename, "r")

    csv_reader = csv.reader(csv_file, delimiter=',')

    uu = []
    dtdt = []
    for row in csv_reader:
        print(row)
        u, dt = row[0].split(";")

        uu.append(u)
        dtdt.append(dt)

    u_values = np.array(uu[1:], dtype=float)
    dt_values = np.array(dtdt[1:], dtype=float)


    def unzip(iterable):
        return zip(*iterable)


    paired = list(zip(u_values, dt_values))
    paired.sort(key=lambda x: x[1])
    print(paired)
    return [u for u, _ in paired], [dt for _, dt in paired]


u_values, dt_values = outfile_to_data("../out.txt")

reference_value = u_values[0]

errors = np.abs(u_values - reference_value)

u_values = np.array(u_values)
dt_values = np.array(dt_values)

plt.loglog(dt_values[1:], errors[1:], '.', label='num')

u_values, dt_values = outfile_to_data("../subiteration_out.txt")

errors = np.abs(u_values - reference_value)

u_values = np.array(u_values)
dt_values = np.array(dt_values)

plt.loglog(dt_values, errors, '^', label='num w subiterations')

u_values, dt_values = outfile_to_data("../wr_out.txt")

errors = np.abs(u_values - reference_value)

u_values = np.array(u_values)
dt_values = np.array(dt_values)

plt.loglog(dt_values, errors, 'x', label='num w waveform')

u_values, dt_values = outfile_to_data("../wr5_out.txt")

errors = np.abs(u_values - reference_value)

u_values = np.array(u_values)
dt_values = np.array(dt_values)

plt.loglog(dt_values, errors, '2', label='num w5 waveform5 (dt excitation = 5 * dt solid solver)')

plt.loglog(dt_values[1:], dt_values[1:], '--', label='O(h^1)')
plt.loglog(dt_values[1:], dt_values[1:]**2, ':', label='O(h^2)')
plt.xlabel("dt solid solver")
plt.ylabel("tip error")
plt.legend()
plt.show()


