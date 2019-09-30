#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) 2012-2019 gmrandazzo@gmail.com
# This file is part of MolecularFramework.
# You can use,modify, and distribute it under
# the terms of the GNU General Public Licenze, version 3.
# See the file LICENSE for details

from rdkit.Chem.Draw import SimilarityMaps
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import string


def Get2DCoordFromMol2(mol2file):
    m = Chem.MolFromMol2File(mol2file)
    Chem.SanitizeMol(m)
    AllChem.Compute2DCoords(m)
    x = []
    y = []
    conf = m.GetConformer()
    for i in range(m.GetNumAtoms()):
        coord = list(conf.GetAtomPosition(i))
        x.append(coord[0])
        y.append(coord[1])
    return m

def main():
    if len(sys.argv) == 3:
        #mol = Chem.MolFromMolFile(sys.argv[1])
        mol = Chem.MolFromMol2File(sys.argv[1])
        if mol is not None:
            Chem.SanitizeMol(mol)
            AllChem.Compute2DCoords(mol)
            contribs = []
            fw = open(sys.argv[2], "r")
            for line in fw:
                v = str.split(line.strip(), " ")
                contribs.append(float(v[2]))
            fw.close()
            fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap='jet', contourLines=10)
            fig.savefig("%s.png" % sys.argv[1].replace(".mol2",""),  bbox_inches='tight')
    else:
        print("Usage: %s  file.mol2 file_atomcontrib.txt")

if __name__ == "__main__":
  #demo()
  main()
