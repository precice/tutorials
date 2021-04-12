#!/usr/bin/python

"""
Use this script to join the visual results (.frd files) of beam1 and beam2 and visualize them.
"""

import subprocess

#### params ####

# size of complete node sets
nsize1 = 201
nsize2 = 81
nsizem = 261  # nsize1 + nsize2 - nsize_coupling_surface

# time steps. set in precice-config.xml and .inp files
nsteps = 50


##################################################################
def join_frd(frd1, frd2):
    """
    Append the nodes and elements from frd2 to those in frd1 and write the result in a new frd.
    """
    with open(frd1, "r") as f1, open(frd2, "r") as f2, open("beam_full.frd", "w") as fp:

        # copy frd header in new file
        for i in xrange(11):
            fp.write(f1.readline())
            f2.readline()

        # node header (change number of nodes)
        line_f1 = f1.readline()
        line_f1 = line_f1[:33] + str(nsizem) + line_f1[36:]
        fp.write(line_f1)
        f2.readline()

        # merging node lines. each iteration in the for uses a new line from frd2. lines in frd1 are advanced manually
        line_f1 = f1.readline()
        for line_f2 in iter(f2.readline, " -3\n"):  # -3 indicates end of line block
            # same node in both files (interface): write any (assuming their values are correct!!)
            if line_f1[:13] == line_f2[:13]:
                if line_f1[2] == "3":
                    continue
                else:
                    fp.write(line_f1)
                    line_f1 = f1.readline()
            # sorting lines according to node index
            elif line_f2 < line_f1:
                fp.write(line_f2)
            else:
                while line_f2 > line_f1:
                    fp.write(line_f1)
                    line_f1 = f1.readline()

                if line_f1[:13] != line_f2[:13]:
                    fp.write(line_f2)

        fp.write(" -3\n")

        # element header (change number of elements)
        line_f1 = f1.readline()
        line_f1 = line_f1[:34] + "32" + line_f1[36:]
        fp.write(line_f1)
        f2.readline()

        # merge element lines. assuming they are sorted and non-overlapping in frd1 and frd2
        for line_f1 in iter(f1.readline, " -3\n"):
            fp.write(line_f1)
        for line_f2 in iter(f2.readline, " -3\n"):
            fp.write(line_f2)
        fp.write(" -3\n")

        # merging blocks of lines for each step
        for i in xrange(nsteps):
            print "step", i + 1
            # step header
            fp.write(f1.readline())
            f2.readline()
            line_f1 = f1.readline()
            line_f1 = line_f1[:33] + str(nsizem) + line_f1[36:]
            fp.write(line_f1)
            f2.readline()
            for j in xrange(5):
                fp.write(f1.readline())
                f2.readline()

            line_f1 = f1.readline()
            for line_f2 in iter(f2.readline, " -3\n"):  # -3 indicates end of line block
                # same node in both files (interface): write any (assuming their values are correct!!)
                if line_f1[:13] == line_f2[:13]:
                    if line_f1[2] == "3":
                        continue
                    else:
                        # this is an interface node for both beams. write the mean of the values in beam1 and beam2
                        mean_vals = [(float(x) + float(y)) / 2. for x, y in zip([line_f1[13:25],
                                                                                 line_f1[25:37], line_f1[37:49]], [line_f2[13:25], line_f2[25:37], line_f2[37:49]])]
                        fp.writelines(line_f1[:13] +
                                      '{:12.5E}'.format(mean_vals[0]) +
                                      '{:12.5E}'.format(mean_vals[1]) +
                                      '{:12.5E}'.format(mean_vals[2]) +
                                      '\n')
                        line_f1 = f1.readline()
                # sorting lines according to node index
                elif line_f2 < line_f1:
                    fp.write(line_f2)
                else:
                    while line_f2[:13] > line_f1[:13]:
                        fp.write(line_f1)
                        line_f1 = f1.readline()

                    if line_f1[:13] != line_f2[:13]:
                        fp.write(line_f2)
                    else:
                        mean_vals = [(float(x) + float(y)) / 2. for x, y in zip([line_f1[13:25],
                                                                                 line_f1[25:37], line_f1[37:49]], [line_f2[13:25], line_f2[25:37], line_f2[37:49]])]
                        fp.writelines(line_f1[:13] +
                                      '{:12.5E}'.format(mean_vals[0]) +
                                      '{:12.5E}'.format(mean_vals[1]) +
                                      '{:12.5E}'.format(mean_vals[2]) +
                                      '\n')
                        line_f1 = f1.readline()

            fp.write(" -3\n")

        fp.write("9999\n")  # EOF


########################## MAIN #################################

join_frd("dirichlet-calculix/beam1.frd", "neumann-calculix/beam2.frd")
subprocess.call(["cgx", "-b", "visualize.fbd"])
