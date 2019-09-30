#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) 2012-2019 gmrandazzo@gmail.com
# This file is part of MolecularFramework.
# You can use,modify, and distribute it under
# the terms of the GNU General Public Licenze, version 3.
# See the file LICENSE for details

import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


def nsplit(s, delim=None):
    return [x for x in s.split(delim) if x]


def SavePC2Mol2File(chgs, fmol2, fout):
    fi = open(fmol2, "r")
    fo = open(fout, "w")
    i = 0
    next_atom = False
    for line in fi:
        if "ATOM" in line:
            next_atom = True
            fo.write(line)
        elif "BOND" in line:
            next_atom = False
            fo.write(line)
        else:
            if next_atom is True:
                v = nsplit(line.strip(), " ")
                fo.write("%7d %-4s       % 3.4f   % 3.4f   % 3.4f %-5s %4d  %4s       % 2.4f\n" %
                         (int(v[0]),
                          v[1],
                          float(v[2]),
                          float(v[3]),
                          float(v[4]),
                          v[5],
                          int(v[6]),
                          v[7],
                          chgs[i]))
                i += 1
            else:
                fo.write(line)
    fo.close()
    fi.close()
    return 0


def CalcEEMCharge(fsdf):
    m = Chem.SDMolSupplier(fsdf, removeHs=False)[0]
    chgs = AllChem.CalcEEMcharges(m)
    # chgs = [x.GetDoubleProp('_GasteigerCharge') for x in m.GetAtoms()]
    return chgs


def CalcGasteigerCharge(fsdf):
    m = Chem.SDMolSupplier(fsdf, removeHs=False)[0]
    AllChem.ComputeGasteigerCharges(m)
    chgs = [x.GetDoubleProp('_GasteigerCharge') for x in m.GetAtoms()]
    return chgs


def main():
    # Calculate charge from sdf and tranfer them to MOL2 file....
    if len(sys.argv) == 4:
        chgs = CalcEEMCharge(sys.argv[1])
        # chgs = CalcGasteigerCharge(sys.argv[1])
        SavePC2Mol2File(chgs, sys.argv[2], sys.argv[3])
    else:
        print("\nUsage %s [SDF input] [MOL2 input] [MOL2 output]\n" % (sys.argv[0]))
    return


if __name__ == '__main__':
    main()
