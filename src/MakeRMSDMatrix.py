#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) 2012-2019 gmrandazzo@gmail.com
# This file is part of MolecularFramework.
# You can use,modify, and distribute it under
# the terms of the GNU General Public Licenze, version 3.
# See the file LICENSE for details

import sys
import subprocess

mol3daligner="/Users/marco/projects/MolDesc/build/src/test3Dmol2aligner"
def WriteTable(data, fout):
  fo = open(fout, "w")
  for line in data:
    for item in line[:-1]:
      fo.write("%s\t" % (item))
    fo.write("%s\n" % (line[-1]))
  fo.close()

def alignmols(mol1, mol2):
    cmd = "%s \"%s\" \"%s\"" % (mol3daligner, mol1, mol2)
    output = ""
    output = subprocess.check_output(cmd, shell=True)
    #return RMSD
    return float(str.split(output.strip(), ":")[-1])

def CalculateRMSDMatrix(mol2lst, outfile):
    rmsdsimmx = []
    for i in range(len(mol2lst)+1):
        rmsdsimmx.append([0.0 for j in range(len(mol2lst)+1)])
    mol_count = len(mol2lst)
    rmsdsimmx[0][0] = "%dx%d matrix" % (mol_count,mol_count)
    for i in range(mol_count):
        rmsdsimmx[i+1][0] = rmsdsimmx[0][i+1] = mol2lst[i].replace(".mol2", "")
    for i in range(mol_count - 1):
        for j in range(i+1, mol_count):
            rmsdsimmx[i+1][j+1] = rmsdsimmx[j+1][i+1] = alignmols(mol2lst[i], mol2lst[j])
        rmsdsimmx[i+1][i+1] = float(0.0)
    WriteTable(rmsdsimmx, ("%s" % (outfile)))

def main():
  if len(sys.argv) > 1:
      CalculateRMSDMatrix(list(sys.argv[1:-1]), sys.argv[-1])
  else:
      print "Usage %s mol2_1 mol2_2 ... mol2_N output.txt" % (sys.argv[0])

if __name__ == '__main__':
    main()
