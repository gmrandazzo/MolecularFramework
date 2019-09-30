#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) 2012-2019 gmrandazzo@gmail.com
# This file is part of MolecularFramework.
# You can use,modify, and distribute it under
# the terms of the GNU General Public Licenze, version 3.
# See the file LICENSE for details

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalForceFields

#m = Chem.MolFromSmiles('N')
m = Chem.MolFromSmiles("O")
#m = Chem.MolFromSmiles("[CH4]")
m2=Chem.AddHs(m)
#m2 = m
AllChem.EmbedMolecule(m2)
AllChem.MMFFOptimizeMolecule(m2)
mp = ChemicalForceFields.MMFFGetMoleculeProperties(m2)

for i in range(m2.GetNumAtoms()):
  print(m2.GetAtoms()[i].GetAtomicNum())
  print(mp.GetMMFFPartialCharge(i))
  print(mp.GetMMFFVdWParams(i, i))
  print("-"*8)
