#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) 2012-2019 gmrandazzo@gmail.com
# This file is part of MolecularFramework.
# You can use,modify, and distribute it under
# the terms of the GNU General Public Licenze, version 3.
# See the file LICENSE for details

'''
@package afmodeltrain

afmodeltrain was writen by Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
Geneve Gen 2015

afmodeltrain will calculate a property by recalculating
each time trough simplex an atom fragment weight.

afmodeltrain need afmodel binary installed with MolDesc package.

'''

import sys
from scipy.optimize import fmin
from scipy.optimize import fmin_cg
from scipy.optimize import fmin_l_bfgs_b
from math import sqrt
import multiprocessing as mp
import commands
import os.path

#import subprocess
#if "check_output" not in dir( subprocess ): # duck punch it in!
    #def f(*popenargs, **kwargs):
        #if 'stdout' in kwargs:
            #raise ValueError('stdout argument not allowed, it will be overridden.')
        #process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        #output, unused_err = process.communicate()
        #retcode = process.poll()
        #if retcode:
            #cmd = kwargs.get("args")
            #if cmd is None:
                #cmd = popenargs[0]
            #raise subprocess.CalledProcessError(retcode, cmd)
        #return output
    #subprocess.check_output = f

def square(x):
    return x*x

afmodel_bin = "/Users/marco/projects/MolDesc/build/src/afmodel"

def GetNCPU():
    return mp.cpu_count()

def afmodel(cmd, output):
    status, outstr = commands.getstatusoutput("%s %s" % (afmodel_bin, cmd))
    #outstr = subprocess.check_output("%s %s" % (afmodel_bin, cmd), shell=True)
    #return float(str.split(output.strip(), ",")[-1])
    output.put(float(str.split(outstr.strip(), ",")[-1]))

def MultiAFModel(mollst, awfile, ncpu):
    """ Multiprocess AFModel """
    results = []
    output = mp.Queue()
    for i in range(0, len(mollst), ncpu):
        processes = []
        for obj in  mollst[i:i+ncpu]:
            if os.path.isfile(obj) == True:
                cmd = ("1 %s %s" % (awfile, obj))
                processes.append(mp.Process(target=afmodel, args=(cmd, output)))
            else:
                print("File %s not found! Verfy your input file!" % (obj))
        for p in processes:
            p.start()
        for p in processes:
            p.join()
        for p in processes:
            results.append(output.get())
        del processes
    return results

class AFModelTrainer(object):
    """ Class to train the model with afmodel """
    def __init__(self, awfile, fmollst):
        self.awfile = awfile
        self.mollst = []
        self.depvar = []
        fi = open(fmollst, "r")
        for line in fi:
            v = str.split(line.strip(), ",")
            if len(v) == 2:
                self.mollst.append(v[0])
                self.depvar.append(float(v[1]))
            else:
                continue
        fi.close()

    def LoadAFWeight(self):
        afw = []
        fi = open(self.awfile, "r")
        for line in fi:
            afw.append(str.split(line.strip(), ","))
        fi.close()
        return afw

    def WriteAFWeight(self, weights):
        afw = self.LoadAFWeight()
        fo = open(self.awfile, "w")
        for i in range(len(afw)):
            str_ = "%s,%s,%f,%s\n" % (afw[i][0], afw[i][1], weights[i], afw[i][3])
            fo.write(str_)
        fo.close()

    def rosen(self, weights):
        self.WriteAFWeight(weights)
        p_depval = MultiAFModel(self.mollst, self.awfile, GetNCPU())
        # Serialization OK
        #p_depval = []
        #for mol in self.mollst:
        #    p_depval.append(afmodel("1 %s %s" % (self.awfile, mol)))

        #sdep = 0.
        #for i in range(len(self.depvar)):
            #sdep += square(self.depvar[i] - p_depval[i])
        #sdep = sqrt(sdep/float(len(self.depvar)))
        #print("SDEP %f" % (sdep))
        #return sdep

        err = 0.
        for i in range(len(self.depvar)):
            err += sqrt(square((self.depvar[i] - p_depval[i])/self.depvar[i]))
        err /= float(len(self.depvar))

        print("Error: %f" % (err*100))
        return err


    def OptimiseAtomWeights(self):
        afw = self.LoadAFWeight()
        init_weight = [float(afw[i][2]) for i in range(len(afw))]

        # l_bfgs_b optimization with constraint
        #bounds = []
        #for i in range(len(afw)):
        #    bounds.append((-10., 10.))
        #final_weights = list(fmin_l_bfgs_b(self.rosen, init_weight, bounds, approx_grad=True, epsilon=0.01)[0])

        #simplex optimization
        final_weights = fmin(self.rosen, init_weight, xtol=1e-6, maxiter=99999, maxfun=99999)
        self.WriteAFWeight(final_weights)

def main():
    if len(sys.argv) != 3:
        print("\nUsage: %s [weight database] [mollst]" % (sys.argv[0]))
        print("\nweight database:")
        print("The weight database must be generated with afmodel.")
        print("\nmollst")
        print("The mollst must be in this format:")
        print("mol1.sdf,0.23")
        print("mol2.sdf,0.54")
        print("mol3.sdf,0.12")
        print("mol4.sdf,0.27")
        print("mol5.sdf,0.93")
        print("............\n")
    else:
        afmtrain = AFModelTrainer(sys.argv[1], sys.argv[2])
        afmtrain.OptimiseAtomWeights()

if __name__ == "__main__":
    main()
