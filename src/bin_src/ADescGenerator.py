#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# (c) 2012-2019 gmrandazzo@gmail.com
# This file is part of MolecularFramework.
# You can use,modify, and distribute it under
# the terms of the GNU General Public Licenze, version 3.
# See the file LICENSE for details

# ADescGen.py: Artificial Descriptor Generator
#              Learn from data the best parameters
import argparse
import sys
import os
from pathlib import Path
from multiprocessing import cpu_count, Pool
from scipy.optimize import minimize as fmin
from sklearn import preprocessing
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.model_selection import RepeatedKFold as SKRepeatedKFold
from sklearn.linear_model import BayesianRidge
from sklearn.ensemble import RandomForestRegressor as RF
from sklearn.svm import SVR
from pyearth import Earth
import numpy as np
import shutil
import subprocess
import random


dgen_bin = "/home/marco/MolDesc/build/src/MolGPot3DDesc"
depot_bin = "/home/marco/MolDesc/build/src/MolEPot3DDesc"
sphdgen_bin = "/home/marco/MolDesc/build/src/MolGPot3DSphericalDesc"
sphdepot_bin = "/home/marco/MolDesc/build/src/MolEPot3DSphericalDesc"
# ljgen_bin = "/home/marco/MolDesc/build/src/MolPot3DDesc"
pls = "/home/marco/QStudioMetrics/build/src/qsm-pls"


def RepeatedKFold(n_splits, n_repeats, dict_target):
    indexes = [i for i in range(len(dict_target))]
    keys = list(dict_target.keys())
    for i in range(n_repeats):
        random.shuffle(indexes)
        for j in range(n_splits):
            train_keys = []
            test_keys = []
            for k in range(len(dict_target)):
                if k % n_splits != j:
                    train_keys.append(keys[indexes[k]])
                else:
                    test_keys.append(keys[indexes[k]])
            yield (train_keys,
                   test_keys)


def CVGroupRead(fcvgroup):
    d = {}
    g = 0
    f = open(fcvgroup, "r")
    for line in f:
        d[g] = list(str.split(line.strip(), ","))
        g += 1
    f.close()
    return d


def StaticGroupCV(d):
    keys = list(d.keys())
    for i in keys:
        test_keys = d[i]
        train_keys = []
        for j in keys:
            if i == j:
                continue
            else:
                train_keys.extend(d[j])
        yield (train_keys, test_keys)


def PLSCrossValidation(descs, targets, lvs=3, modpath="pls_mod"):
    fx = open("x.txt", "w")
    fy = open("y.txt", "w")
    names = []
    for key in targets.keys():
        if key in descs.keys():
            names.append(key)
            d = list(descs[key])
            for j in range(len(d)-1):
                fx.write("%.10f," % (float(d[j])))
            fx.write("%.10f\n" % (float(d[-1])))
            fy.write("%.4f\n" % (float(targets[key])))
        else:
            continue
    fx.close()
    fy.close()

    if Path(modpath).is_dir() is True:
        shutil.rmtree(modpath, ignore_errors=True)

    cmd = "%s -rmodel -x x.txt -y y.txt -o pls_mod -c %d -xa 1 -s ','" % (pls, lvs)
    os.system(cmd)
    cmd = "%s -cv -x x.txt -y y.txt -xa 1 -dm pls_mod -g 5 -i 20 -c %d -s ','" % (pls, lvs)
    os.system(cmd)
    ycv = {}
    f = open("%s/Validated_Predicted_Y.txt" % (modpath), "r")
    i = 0
    for line in f:
        v = str.split(line.strip())
        ycv[names[i]] = []
        for item in v:
            if "nan" in item:
                ycv[names[i]].append(0.)
            else:
                ycv[names[i]].append(float(item))
        i = i + 1
    # Cleanup!
    shutil.rmtree(modpath, ignore_errors=True)
    return ycv


def BayesianRidgeRegressionCrossValidation(descs, targets):
    print("Bayesian")
    x = []
    y = []
    names = []
    for key in targets.keys():
        if key in descs.keys():
            names.append(key)
            d = list(descs[key])
            x.append(d)
            y.append(targets[key])
        else:
            continue
    x = np.array(x).astype(float)
    y = np.array(y).astype(float)
    ycv = {}
    for name in names:
        ycv[name] = []

    rkf = SKRepeatedKFold(n_splits=5, n_repeats=20, random_state=2652124)
    for train_index, test_index in rkf.split(x):
        # print("TRAIN:", train_index, "TEST:", test_index)
        X_train, X_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        scaler = preprocessing.StandardScaler().fit(np.array(X_train))
        X_train_scaled = scaler.transform(np.array(X_train))
        # clf = Earth()
        # clf = RF(max_depth=2, random_state=0, n_estimators=100)
        clf = BayesianRidge()
        # clf = SVR(kernel='rbf', C=100, gamma=0.1, epsilon=.1)
        # clf = SVR(kernel='rbf')
        clf.fit(X_train_scaled, y_train)
        X_test_scaled = scaler.transform(np.array(X_test))
        y_pred = clf.predict(X_test_scaled)
        for i in range(len(test_index)):
            ycv[names[test_index[i]]].append(y_pred[i])
    return ycv


def WriteParmFile(params, out):
    f = open(out, "w")
    for key in params.keys():
        f.write("%s %e\n" % (key, params[key]))
    f.close()


def DescGen_((mol, pfile, g, s, gsmin, gsmax, esmin, esmax, nrot)):
    cmd = "%s %s %s %d %d %f %f %d" % (dgen_bin, mol, pfile, g, s, gsmin, gsmax, nrot)
    gen_desc = subprocess.check_output(cmd, shell=True)
    cmd = "%s %s %d %d %f %f %d" % (depot_bin, mol, g, s, esmin, esmax, nrot)
    epot_desc = subprocess.check_output(cmd, shell=True)
    gdesc = str.split(str(gen_desc.decode('UTF-8').rstrip()), "\t")
    edesc = str.split(str(epot_desc.decode('UTF-8').rstrip()), "\t")[1:]
    return list(gdesc+edesc)


def DescGen(mollst,
            params,
            grid=25,
            size=25,
            gpot_scale_min=-10.,
            gpot_scale_max=10.,
            epot_scale_min=-0.5,
            epot_scale_max=0.5,
            nrot=10):
    WriteParmFile(params, "parm.txt")

    """
    for mol in mollst:
        # molname = str(Path(mol).resolve().stem.split(".")[0])
        cmd = "%s %s parm.txt %d %d" % (dgen_bin, mol, grid, size)
        out = subprocess.check_output(cmd, shell=True)
        row = str.split(str(out.decode('UTF-8').rstrip()), "\t")
        descs[row[0]] = row[1:]
    """

    p = Pool(processes=cpu_count())
    res = p.map(DescGen_, [(mol, "parm.txt", grid, size, gpot_scale_min, gpot_scale_max, epot_scale_min, epot_scale_max, nrot) for mol in mollst])
    p.close()
    descs = {}
    for row in res:
        descs[row[0]] = row[1:]
    return descs


def SphericalDescGen_((mol, pfile, npnt, gsmin, gsmax, esmin, esmax, nrot)):
    cmd = "%s %s %s %d %f %f %d" % (sphdgen_bin, mol, pfile, npnt, gsmin, gsmax, nrot)
    gen_desc = subprocess.check_output(cmd, shell=True)
    cmd = "%s %s %d %f %f %d" % (sphdepot_bin, mol, npnt, esmin, esmax, nrot)
    epot_desc = subprocess.check_output(cmd, shell=True)
    gdesc = str.split(str(gen_desc.decode('UTF-8').rstrip()), "\t")
    edesc = str.split(str(epot_desc.decode('UTF-8').rstrip()), "\t")[1:]
    return list(gdesc+edesc)


def SphericalDescGen(mollst,
                     params,
                     npnt=100,
                     gpot_scale_min=-10.,
                     gpot_scale_max=10.,
                     epot_scale_min=-0.5,
                     epot_scale_max=0.5,
                     nrot=10):
    WriteParmFile(params, "parm.txt")

    p = Pool(processes=cpu_count())
    res = p.map(SphericalDescGen_, [(mol, "parm.txt", npnt, gpot_scale_min, gpot_scale_max, epot_scale_min, epot_scale_max, nrot) for mol in mollst])
    p.close()
    descs = {}
    for row in res:
        descs[row[0]] = row[1:]
    return descs


def ErrorCalculator(descs, target, modpath, nlv=3, method="pls"):
    ycvp = None
    if method == "pls":
        ycvp = PLSCrossValidation(descs, target, nlv, modpath)
    else:
        ycvp = BayesianRidgeRegressionCrossValidation(descs, target)

    yt = []
    yp = []
    mse = []
    mae = []
    real_nlv = len(ycvp.values()[0])
    for lv in range(real_nlv):
        for key in ycvp.keys():
            yt.append(target[key])
            yp.append(ycvp[key][lv])
        mae.append(mean_absolute_error(yt, yp))
        mse.append(mean_squared_error(yt, yp))

    # Check if the mse
    # decrease by adding more latent variables
    for i in range(1, len(mse)):
        if mse[i-1] < mse[i]:
            return mse[i-1], mae[i-1]
        else:
            continue
    return mse[-1], mae[-1]


class ADescGen(object):
    def __init__(self,
                 mol2_dir,
                 csv_target,
                 other_descs=None,
                 uparams=None,
                 descmethod="voxel_based",
                 regmethod="pls"):
        self.regmethod = str(regmethod)
        self.descmethod = str(descmethod)
        self.mol2_dir = str(Path(mol2_dir).absolute())
        self.mol2lst = []
        for p in Path(self.mol2_dir).iterdir():
            if p.is_file() and ".mol2" in str(p):
                self.mol2lst.append(str(p))
            else:
                continue

        self.target = self.ReadTarget(csv_target)
        self.modpath = "pls_mod"
        self.nlv = 4

        # Electrostatic potential and Generic Descriptor variables
        # 50 molecular rotation
        self.nrot = 50
        # Used for the spherical method
        self.npnt = 50
        self.grid_size = 25
        self.step_size = 25

        # Generic potential descriptor variables
        #self.finit = 0.5
        self.finit = 1e-8
        #self.finit = random.uniform
        self.gsmin = -10
        self.gsmax = 10

        # Electrostatic potential descriptor variables
        self.esmin = -0.5
        self.esmax = 0.5

        self.params = None
        self.uparams = None
        if uparams is not None:
            self.uparams = self.readParams(uparams)
            self.params = self.uparams
        else:
            self.params = self.initParams()

        self.other_descs = None
        if other_descs is not None:
            self.other_descs = self.ReadOtherDescs(other_descs)

        self.orig_mse = None
        self.orig_mae = None
        self.train_test_molecules = None
        self.val_molecules = None

    def ReadTarget(self, csv_target):
        fi = open(csv_target, "r")
        d = {}
        for line in fi:
            v = str.split(line.strip(), ",")
            if "Molecule" in line:
                continue
            else:
                d[v[0]] = float(v[1])
        fi.close()
        return d

    def ReadOtherDescs(self, csv_target):
        fi = open(csv_target, "r")
        d = {}
        for line in fi:
            v = str.split(line.strip(), ",")
            if "Molecule" in line:
                continue
            else:
                d[v[0]] = []
                for item in v[1:]:
                    d[v[0]].append(float(item))
        fi.close()
        return d

    def readParams(self, fparams):
        d = {}
        f = open(fparams, "r")
        for line in f:
            v = str.split(line.strip(), " ")
            d[v[0]] = float(v[1])
        f.close()
        return d

    def initParams(self):
        params = {}
        random.seed(len(self.mol2lst))
        for mol in self.mol2lst:
            f = open(mol, "r")
            for line in f:
                v = list(filter(None, str.split(line.strip(), " ")))
                if len(v) != 9:
                    continue
                else:
                    if v[5] not in params.keys():
                        # params[v[5]] = random.uniform(-1, 1)
                        # params[v[5]] = 1e-4
                        if callable(self.finit) is True:
                            params[v[5]] = self.finit(-1., 1.)
                        else:
                            params[v[5]] = self.finit
                    else:
                        continue
            f.close()
        return params

    def f(self, x):
        """
        1) Create a database
        2) Train the network
        3) Return the score result
        """
        new_params = {}
        i = 0
        for key in self.params.keys():
            new_params[key] = x[i]
            i += 1
        """
        for key in self.params.keys():
            print("%s o: %e n: %e" % (key, self.params[key], new_params[key]))
        """
        # descs = DescGen(self.mol2lst, new_params, self.grid_size, self.step_size, self.gsmin, self.gsmax, self.esmin, self.esmax)
        descs = None
        if self.descmethod == "voxel_based":
            descs = DescGen(self.train_test_molecules, new_params, self.grid_size, self.step_size, self.gsmin, self.gsmax, self.esmin, self.esmax, self.nrot)
        else:
            descs = SphericalDescGen(self.train_test_molecules, new_params, self.npnt, self.gsmin, self.gsmax, self.esmin, self.esmax, self.nrot)

        print(len(descs.keys()))
        if self.other_descs is not None:
            for key in self.other_descs.keys():
                if key in descs.keys():
                    descs[key].extend(self.other_descs[key])
                else:
                    continue

        mse, mae = ErrorCalculator(descs, self.target, self.modpath, self.nlv, self.regmethod)

        if self.orig_mse is None:
            self.orig_mse = mse
            self.orig_mae = mae

        print("Original MSE: %.8f MAE: %.8f" % (self.orig_mse, self.orig_mae))
        print("Actual MSE: %.8f MAE: %.8f" % (mse, mae))
        return mse

    def cv(self, fout, n_split=5, n_repeats=20, fcvgroup=None):
        if self.other_descs is not None:
            self.orig_mse, self.orig_mae = ErrorCalculator(self.other_descs, self.target, self.modpath, self.nlv, self.regmethod)

        cvmethod = None
        if fcvgroup is not None:
            cvgroups = CVGroupRead(fcvgroup)
            cvmethod = StaticGroupCV(cvgroups)
        else:
            cvmethod = RepeatedKFold(n_split, n_repeats, self.target)

        # Final parameters
        fparams = {}
        for key in self.params.keys():
            fparams[key] = []

        nval = 0
        val_errors = []
        for dataset_keys, val_keys in cvmethod:
            # Reinitialize the parameters
            if self.uparams is not None:
                self.params = self.uparams
            else:
                self.params = self.initParams()

            self.train_test_molecules = []
            for mol in dataset_keys:
                self.train_test_molecules.append("%s/%s.mol2"% (self.mol2_dir, mol))

            #res = fmin(self.f, list(self.params.values()), method='Nelder-Mead', tol=1e-3)
            res = fmin(self.f, list(self.params.values()), method='Powell', tol=1e-3)
            #res = fmin(self.f, list(self.params.values()), method='L-BFGS-B', options = {'eps': 1e-3})

            i = 0
            tmpparms = {}
            for key in self.params.keys():
                fparams[key].append(res.x[i])
                tmpparms[key] = res.x[i]
                i += 1
            if len(val_keys) > 2:
                valmol = []
                for mol in val_keys:
                    valmol.append("%s/%s.mol2"% (self.mol2_dir, mol))
                descs = None
                if self.descmethod is "voxel_based":
                    descs = DescGen(valmol, tmpparms, self.grid_size, self.step_size, self.gsmin, self.gsmax, self.esmin, self.esmax, self.nrot)
                else:
                    descs = SphericalDescGen(valmol, tmpparms, self.npnt, self.gsmin, self.gsmax, self.esmin, self.esmax, self.nrot)
                if self.other_descs is not None:
                    for key in descs.keys():
                        if key in self.other_descs.keys():
                            descs[key].extend(self.other_descs[key])
                        else:
                            continue
                mse, mae = ErrorCalculator(descs, self.target, self.modpath, self.nlv, self.regmethod)
                print("Validation MSE: %.8f MAE: %.8f" % (mse, mae))
                val_errors.append([mse, mae, val_keys])
            else:
                val_errors.append([0.0, 0.0, val_keys])
            nval += 1

        for i in range(nval):
            p = {}
            for key in fparams.keys():
                p[key] = fparams[key][i]
            WriteParmFile(p, "cv%d_%s" % (i,fout))

        print("Final Results")
        for i in range(nval):
            print(val_errors[i][0], val_errors[i][1], val_errors[i][2])

    def train(self, fout):
        if self.other_descs is not None:
            self.orig_mse, self.orig_mae = ErrorCalculator(self.other_descs, self.target, self.modpath, self.nlv, self.regmethod)
        self.train_test_molecules = self.mol2lst
        #res = fmin(self.f, list(self.params.values()), method='Nelder-Mead', tol=1e-3)
        res = fmin(self.f, list(self.params.values()), method='Powell', tol=1e-3)
        #res = fmin(self.f, list(self.params.values()), method='L-BFGS-B', options = {'eps': 1e-3})

        fparams = {}
        i = 0
        for key in self.params.keys():
            fparams[key] = res.x[i]
            i += 1
        WriteParmFile(self.params, "original_params.txt")
        WriteParmFile(fparams, fout)

        descs = None
        if self.descmethod is "voxel_based":
            descs = DescGen(self.mol2lst, fparams, self.grid_size, self.step_size, self.gsmin, self.gsmax, self.esmin, self.esmax, self.nrot)
        else:
            descs = SphericalDescGen(self.mol2lst, fparams, self.npnt, self.gsmin, self.gsmax, self.esmin, self.esmax, self.nrot)

        if self.other_descs is not None:
            for key in self.other_descs.keys():
                descs[key].extend(self.other_descs[key])

        mse, mae = ErrorCalculator(descs, self.target, self.modpath, self.nlv, self.regmethod)
        print("Final MSE: %.8f MAE: %.8f" % (mse, mae))



def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mol2dir',
                   default=None,
                   type=str,
                   help='Mol2 file directory')
    p.add_argument('--desc_csv',
                   default=None,
                   type=str,
                   help='Molecule descriptors')
    p.add_argument('--target',
                   default=None,
                   type=str,
                   help='Target to predict')
    p.add_argument('--initparams',
                   default=None,
                   type=str,
                   help='Init parameterss')
    p.add_argument('--fout',
                   default=None,
                   type=str,
                   help='Out parameters')
    p.add_argument('--cv',
                   default=False,
                   type=bool,
                   help='Run Cross Validation')
    p.add_argument('--n_splits',
                   default=5,
                   type=int,
                   help='Number of Cross Validation splits')
    p.add_argument('--n_repeats',
                   default=20,
                   type=int,
                   help='Number of Cross Validation splits')
    p.add_argument('--cvgroupfile',
                   default=None,
                   type=str,
                   help='Static cross-validation group file')
    p.add_argument('--regmethod',
                   default="pls",
                   type=str,
                   help='pls or bayesian')
    p.add_argument('--descmethod',
                   default="voxel_based",
                   type=str,
                   help='voxel_based or spherical')
    args = p.parse_args(sys.argv[1:])

    if args.fout is None and args.mol2dir is None and args.fout is None:
        print("\nUsage %s --mol2dir [mol2_dir] --target [csv_target] --desc_csv [csv_other_descs (optional)] --initparams [init params (optional)] --cv [True/False (optional)] --cvgroupfile [file csv (optional)] --n_splits [default=5 (optional)] --n_repeats [default=20 (optional)]" % (sys.argv[0]))
    else:
        adgen = ADescGen(args.mol2dir,
                         args.target,
                         args.desc_csv,
                         args.initparams,
                         str(args.descmethod),
                         str(args.regmethod))
        if args.cv is True:
            adgen.cv(args.fout, args.n_splits, args.n_repeats, args.cvgroupfile)
        else:
            adgen.train(args.fout)



if __name__ in "__main__":
    main()
