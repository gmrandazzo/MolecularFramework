#!/usr/bin/env python3

import sys
from time import sleep

def get_relatom_mass(aname, e):
    r = 0.
    for row in e:
        mass = float(row[0].split("(")[0])
        abund = float(row[1].split("(")[0])
        r += mass*abund
    print('{"%s", %.12f},' % (aname, r))


def get_monoisotopic_mass(aname, e):
    monomass = float(e[0][0].split("(")[0])
    cabund = float(e[0][1].split("(")[0])
    for i in range(1, len(e)):
        abund = float(e[i][1].split("(")[0])
        if abund > cabund:
            cabund = abund
            monomass = float(e[i][0].split("(")[0])
    print('%s: %.12f},' %(aname, monomass))


def main():
    e = []
    aname = ""
    f = open(sys.argv[1], "r")
    for line in f:
        v = str.split(line.strip(), " ")
        v = list(filter(None, v))
        if len(v) >= 5:
            #first line
            if len(e) > 0:
                try:
                    get_monoisotopic_mass(aname, e)
                    #get_relatom_mass(aname, e)
                except:
                    print(e)
                    print("Problem with %s" % (aname))
            aname = v[1]
            e = [[v[3], v[4]]]
        else:
            if len(v) >= 3:
                e.append([v[1], v[2]])


if __name__ in "__main__":
    main()

