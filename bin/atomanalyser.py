#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) 2012-2019 gmrandazzo@gmail.com
# This file is part of MolecularFramework.
# You can use,modify, and distribute it under
# the terms of the GNU General Public Licenze, version 3.
# See the file LICENSE for details


import sys
from rdkit import Chem


class ATYPE(object):
    """ Data structure to define the atom type table """

    def __init__(self, aprops, patt):
        """ init construct which define:
        - adjacence matrix wich contain atoms hybridisation (adj)
        - atom properties (aprops) which contain the atom charge and
          the van der waals alpha etc..
        """
        self.aprops = aprops
        self.patt = patt


def ReadATypes(fatypes):
    """Read atom type table"""
    atab = []
    f = open(fatypes, "r")
    for line in f:
        # 1.00 2.00 ; {'C:SP2': ['O:SP2']} || {'C:SP2': ['S:SP2']} # sp2 carbon: C=O and C=S
        v = str.split(line.strip(), ";")
        if len(v) == 2:
            aprops = str.strip(v[0], " ")
            frgs = str.split(v[1], "||")
            for frg in frgs:
                patt = Chem.MolFromSmarts(frg)
                atab.append(ATYPE(aprops, patt))
        else:
            continue
    f.close()
    return atab


def AtomAnalyser(m, atab):
    m_graph = nx.Graph()
    for i in range(m.GetNumBonds()):
        from_ = m.GetBondWithIdx(i).GetBeginAtomIdx()
        to_ = m.GetBondWithIdx(i).GetEndAtomIdx()
        atm = ""
        # print m.GetAtomWithIdx(to_).GetSymbol(), m.GetAtomWithIdx(to_).GetHybridization()
        if m.GetAtomWithIdx(to_).GetIsAromatic() == True:
            atm = "%s_Ar" % (m.GetAtomWithIdx(to_).GetSymbol())
        else:
            atm = "%s_%s" % (m.GetAtomWithIdx(to_).GetSymbol(),
                             m.GetAtomWithIdx(to_).GetHybridization())
        m_graph.add_edge(from_, to_, weight=atypes[atm.strip()])

    print nx.adjacency_matrix(m_graph).todense()
    for frg in atab:
        GM = iso.GraphMatcher(m_graph, frg.G, edge_match=my_edge_match)
        if GM.subgraph_is_isomorphic():
            print frg.atype, GM.mapping
        # for mapping in GM.subgraph_isomorphisms_iter():
        #    print frg.atype, mapping

    # m.GetAtomWithIdx(i).GetTotalNumHs()
    # m.GetAtomWithIdx(i).GetIsAromatic()
    # m.GetAtomWithIdx(i).GetHybridization()
    # m.GetAtomWithIdx(i).GetNeighbors()
    return 0


def main():
    atab = ReadATypes(sys.argv[2])
    m = Chem.MolFromMol2File(sys.argv[1], removeHs=False)
    AtomAnalyser(m, atab)


if __name__ == '__main__':
    main()
