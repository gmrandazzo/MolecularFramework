#/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) 2012-2019 gmrandazzo@gmail.com
# This file is part of MolecularFramework.
# You can use,modify, and distribute it under
# the terms of the GNU General Public Licenze, version 3.
# See the file LICENSE for details

import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count
from ADescGenerator import DescGen_


def PrepareMatrix(mollst, fother_desc, fparm, grid_size, step_size, gsmin, gsmax, esmin, esmax, nrot):
    # Read the user descriptors
    other_desc = {}
    other_desc_header = []

    if Path(fother_desc).exists():
        f = open(fother_desc, "r")
        for line in f:
            if "Molecule" in line:
                other_desc_header.append(str.split(line.strip(), ",")[1:])
            else:
                v = str.split(line.strip(), ",")
                other_desc[v[0]] = v[1:]
        f.close()

    # Compute the automatic descriptors

    p = Pool(processes=cpu_count())

    #def DescGen_((mol, pfile, g, s, gsmin, gsmax, esmin, esmax)):
    res = p.map(DescGen_, [(mol, str(fparm), int(grid_size), int(step_size), float(gsmin), float(gsmax), float(esmin), float(esmax), int(nrot)) for mol in mollst])
    p.close()
    descs = {}
    for row in res:
        descs[row[0]] = row[1:]

    # Write the final data matrix for prediction
    m = {}
    if Path(fother_desc).exists():
        for key in other_desc.keys():
            if key in descs.keys():
                m[key] = []
                d = list(descs[key])
                for j in range(len(d)):
                    m[key].append(float(d[j]))
                for j in range(len(other_desc[key])):
                    m[key].append(float(other_desc[key][j]))
            else:
                continue
    else:
        for key in descs.keys():
            m[key] = []
            d = list(descs[key])
            for j in range(len(d)):
                m[key].append(float(d[j]))

    header = []
    ndesc = int(len(descs.values()[0])/2.)
    for i in range(ndesc):
        header.append("EPOTDESC%d" % (i))
    for i in range(ndesc):
        header.append("GPOTDESC%d" % (i))
    if Path(fother_desc).exists():
        header.extend(other_desc_header)

    return m, header


def main():
    if len(sys.argv) > 10:
        # PrepareMatrix(mollst, fother_desc, fparm, grid_size, step_size, gsmin, gsmax, esmin, esmax, nrot):
        m, h =  PrepareMatrix(sys.argv[1:-10],
                              sys.argv[-10],
                              sys.argv[-9],
                              sys.argv[-8],
                              sys.argv[-7],
                              sys.argv[-6],
                              sys.argv[-5],
                              sys.argv[-4],
                              sys.argv[-3],
                              sys.argv[-2],)
        f = open(sys.argv[-1], "w")
        f.write("Molecule,")
        for i in range(len(h)-1):
            f.write("%s," % (h[i]))
        f.write("%s\n" % (h[-1]))
        for key in m.keys():
            f.write("%s," % (key))
            for i in range(len(m[key])-1):
                f.write("%.8f," % (m[key][i]))
            f.write("%.8f\n" % (m[key][-1]))
        f.close()
    else:
        print("\nUsage: %s [mol2] [descriptors] [generic pot parm] [grid size] [step size] [generic scale min] [generic scale max] [epot scale min] [epot scale max] [output file] [n rotation] [output file]" % (sys.argv[0]))

if __name__ in "__main__":
    main()
