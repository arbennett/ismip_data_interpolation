"""
Python utility to produce a master data set from the interpolated data sets 
produced by the interpolate_ismip.py script.  Two potential files are created,
one for each lmla and stokes models.

@author arbennett
"""

import os
import sys
import csv
import glob
import argparse
import fnmatch
import numpy as np
from itertools import product
from scipy import interpolate

def parse_args(args):
    parser = argparse.ArgumentParser(description="ISMIP data interpolation utility.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            fromfile_prefix_chars='@')
    parser.add_argument('-i', '--input', default=os.path.join('.', 'ismip_interp'),
                        help='The location of the ISMIP data directory.')
    parser.add_argument('--lmla', action='store_true',
                        help="Analyze LMLA models only.")
    parser.add_argument('--stokes', action='store_true',
                        help="Analyze Full Stokes models only.")
    return parser.parse_args(args)


def get_stats(var):
    mu = np.mean(var, axis=0)
    sigma = np.std(var, axis=0)
    mini = np.min(var, axis=0)
    maxi = np.max(var,axis=0)
    return (mu, sigma, mini, maxi)


def assimilate_exp_a(file_list, model, opts):
    print("Beginning analysis of experiment A for " + model + " models....")
    res_80 = fnmatch.filter(file_list, '?????080_interp.txt')
    res_20 = fnmatch.filter(file_list, '?????020_interp.txt')
    header = ['x', 'y', 'vs_x mean', 'vs_x stdev', 'vs_x min', 'vs_x max', 'vs_y mean',
              'vs_y stdev', 'vs_y min', 'vs_y max', 'tb_xy mean', 'tb_xy stdev', 
              'tb_xy min', 'tb_xy max', 'tb_yz mean', 'tb_yz stdev', 'tb_yz min',
              'tb_yz max']

    X, Y, VS_X, VS_Y, VS_Z, TB_XY, TB_YZ, DP = [],[],[],[],[],[],[],[]
    for fname in res_80:
        fpath = os.path.join(opts.input, model,fname)
        x, y, vs_x, vs_y, vs_z, tb_xy, tb_yz, dp = np.loadtxt(fpath, unpack=True)
        X.append(x)
        Y.append(y)
        VS_X.append(vs_x)
        VS_Y.append(vs_y)
        VS_Z.append(vs_z)
        TB_XY.append(tb_xy)
        TB_YZ.append(tb_yz)
        DP.append(dp)
    
    X = X[0] # These should all be the same so just grab the first
    Y = Y[0]
    VS_X = get_stats(VS_X)
    VS_Y = get_stats(VS_Y)
    VS_Z = get_stats(VS_Z)
    TB_XY = get_stats(TB_XY)
    TB_YZ = get_stats(TB_YZ)
    DP = get_stats(DP)

    data = np.transpose([X,Y,VS_X[0], VS_X[1], VS_X[2], VS_X[3],
        VS_Y[0], VS_Y[1], VS_Y[2], VS_Y[3], VS_Z[0], VS_Z[1], VS_Z[2], VS_Z[3],
        TB_XY[0], TB_XY[1], TB_XY[2], TB_XY[3], TB_YZ[0], TB_YZ[1], 
        TB_YZ[2], TB_YZ[3], DP[0], DP[1], DP[2], DP[3]])
    fmt = "%.6f"
    with open('ismip_a_80_' + model +'_combined.txt','wb') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        np.savetxt(f, data, fmt=fmt, delimiter=',')
    
    X, Y, VS_X, VS_Y, VS_Z, TB_XY, TB_YZ, DP = [],[],[],[],[],[],[],[]
    for fname in res_20:
        fpath = os.path.join(opts.input, model,fname)
        x, y, vs_x, vs_y, vs_z, tb_xy, tb_yz, dp = np.loadtxt(fpath, unpack=True)
        X.append(x)
        Y.append(y)
        VS_X.append(vs_x)
        VS_Y.append(vs_y)
        VS_Z.append(vs_z)
        TB_XY.append(tb_xy)
        TB_YZ.append(tb_yz)
        DP.append(dp)

    X = X[0] # These should all be the same so just grab the first
    Y = Y[0]
    VS_X = get_stats(VS_X)
    VS_Y = get_stats(VS_Y)
    VS_Z = get_stats(VS_Z)
    TB_XY = get_stats(TB_XY)
    TB_YZ = get_stats(TB_YZ)
    DP = get_stats(DP)

    data = np.transpose([X,Y,VS_X[0], VS_X[1], VS_X[2], VS_X[3],
        VS_Y[0], VS_Y[1], VS_Y[2], VS_Y[3], VS_Z[0], VS_Z[1], VS_Z[2], VS_Z[3],
        TB_XY[0], TB_XY[1], TB_XY[2], TB_XY[3], TB_YZ[0], TB_YZ[1], 
        TB_YZ[2], TB_YZ[3], DP[0], DP[1], DP[2], DP[3]])
    fmt = "%.6f"
    with open('ismip_a_20_' + model + '_combined.txt','wb') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        np.savetxt(f, data, fmt=fmt, delimiter=',')
 

def assimilate_exp_c(file_list, model, opts):
    print("Beginning analysis of experiment C for " + model + " models....")
    res_80 = fnmatch.filter(file_list, '?????080_interp.txt')
    res_20 = fnmatch.filter(file_list, '?????020_interp.txt')
    header = ['x', 'y', 'vs_x mean', 'vs_x stdev', 'vs_x min', 'vs_x max', 'vs_y mean',
              'vs_y stdev', 'vs_y min', 'vs_y max', 'vs_z mean', 'vs_z stdev', 'vs_z min', 'vs_z max', 
              'vb_x mean', 'vb_x stdev', 'vb_x min', 'vb_x max', 'vb_y mean', 
              'vb_y stdev', 'vb_y min', 'vb_y max', 'tb_xy mean', 'tb_xy stdev', 
              'tb_xy min', 'tb_xy max', 'tb_yz mean', 'tb_yz stdev', 'tb_yz min',
              'tb_yz max']

    X, Y, VS_X, VS_Y, VS_Z, VB_X, VB_Y, TB_XY, TB_YZ, DP = [], [], [], [], [], [], [], [], [], []
    for fname in res_20:
        fpath = os.path.join(opts.input, model,fname)
        x, y, vs_x, vs_y, vs_z, vb_x, vb_y, tb_xy, tb_yz, dp = np.loadtxt(fpath, unpack=True)
        X.append(x)
        Y.append(y)
        VS_X.append(vs_x)
        VS_Y.append(vs_y)
        VS_Z.append(vs_z)
        VB_X.append(vb_x)
        VB_Y.append(vb_y)
        TB_XY.append(tb_xy)
        TB_YZ.append(tb_yz)
        DP.append(dp)
    
    X = X[0] # These should all be the same so just grab the first
    Y = Y[0]
    VS_X = get_stats(VS_X)
    VS_Y = get_stats(VS_Y)
    VS_Z = get_stats(VS_Z)
    VB_X = get_stats(VB_X)
    VB_Y = get_stats(VB_Y)
    TB_XY = get_stats(TB_XY)
    TB_YZ = get_stats(TB_YZ)
    DP = get_stats(DP)

    data = np.transpose([X,Y,VS_X[0], VS_X[1], VS_X[2], VS_X[3],
        VS_Y[0], VS_Y[1], VS_Y[2], VS_Y[3], VS_Z[0], VS_Z[1], VS_Z[2], VS_Z[3],
        VB_X[0], VB_X[1], VB_X[2], VB_X[3], VB_Y[0], VB_Y[1], VB_Y[2], VB_Y[3],
        TB_XY[0], TB_XY[1], TB_XY[2], TB_XY[3], TB_YZ[0], TB_YZ[1], 
        TB_YZ[2], TB_YZ[3], DP[0], DP[1], DP[2], DP[3]])
    fmt = "%.6f"
    with open('ismip_c_80_' + model + '_combined .txt','wb') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        np.savetxt(f, data, fmt=fmt, delimiter=',')

    X, Y, VS_X, VS_Y, VS_Z, VB_X, VB_Y, TB_XY, TB_YZ, DP = [], [], [], [], [], [], [], [], [], []
    for fname in res_20:
        fpath = os.path.join(opts.input, model,fname)
        x, y, vs_x, vs_y, vs_z, vb_x, vb_y, tb_xy, tb_yz, dp = np.loadtxt(fpath, unpack=True)
        X.append(x)
        Y.append(y)
        VS_X.append(vs_x)
        VS_Y.append(vs_y)
        VS_Z.append(vs_z)
        VB_X.append(vb_x)
        VB_Y.append(vb_y)
        TB_XY.append(tb_xy)
        TB_YZ.append(tb_yz)
        DP.append(dp)

    X = X[0] # These should all be the same so just grab the first
    Y = Y[0]
    VS_X = get_stats(VS_X)
    VS_Y = get_stats(VS_Y)
    VS_Z = get_stats(VS_Z)
    VB_X = get_stats(VB_X)
    VB_Y = get_stats(VB_Y)
    TB_XY = get_stats(TB_XY)
    TB_YZ = get_stats(TB_YZ)
    DP = get_stats(DP)

    data = np.transpose([X,Y,VS_X[0], VS_X[1], VS_X[2], VS_X[3],
        VS_Y[0], VS_Y[1], VS_Y[2], VS_Y[3], VS_Z[0], VS_Z[1], VS_Z[2], VS_Z[3],
        VB_X[0], VB_X[1], VB_X[2], VB_X[3], VB_Y[0], VB_Y[1], VB_Y[2], VB_Y[3],
        TB_XY[0], TB_XY[1], TB_XY[2], TB_XY[3], TB_YZ[0], TB_YZ[1], 
        TB_YZ[2], TB_YZ[3], DP[0], DP[1], DP[2], DP[3]])
    fmt = "%.6f"
    with open('ismip_c_20_' + model + '_combined.txt','wb') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        np.savetxt(f, data, fmt=fmt, delimiter=',')
 

def assimilate_exp_f(file_list, model, opts):
    print("Beginning analysis of experiment F for " + model + " models....")
    res_100 = fnmatch.filter(file_list, '?????00?_interp.txt')
    header = ['x', 'y', 'vs_x mean', 'vs_x stdev', 'vs_x min', 'vs_x max', 'vs_y mean',
              'z mean', 'z stdev', 'z min', 'z max', 'vs_y stdev', 'vs_y min', 
              'vs_y max', 'vs_z mean', 'vs_z stdev', 'vs_z min', 'vs_z max']
   
    X, Y, Z, VS_X, VS_Y, VS_Z = [],[],[],[],[],[]
    for fname in res_100:
        fpath = os.path.join(opts.input, model,fname)
        x, y, z, vs_x, vs_y, vs_z = np.loadtxt(fpath, unpack=True)
        X.append(x)
        Y.append(y)
        Z.append(z)
        VS_X.append(vs_x)
        VS_Y.append(vs_y)
        VS_Z.append(vs_z)

    X = X[0] # These should all be the same so just grab the first
    Y = Y[0]
    Z = get_stats(Z)
    VS_X = get_stats(VS_X)
    VS_Y = get_stats(VS_Y)
    VS_Z = get_stats(VS_Z)

    data = np.transpose([X,Y,VS_X[0], VS_X[1], VS_X[2], VS_X[3],
        VS_Y[0], VS_Y[1], VS_Y[2], VS_Y[3], VS_Z[0], VS_Z[1], VS_Z[2], VS_Z[3]])
    fmt = "%.6f"
    with open('ismip_f_100_' + model + '_combined.txt','wb') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        np.savetxt(f, data, fmt=fmt, delimiter=',')


def analyze_lmla(opts):
    """
    Analyzes all of the LMLA based models 
    """
    data_a, data_c, data_f = [], [], []
    data_dir = os.path.join(opts.input, 'lmla')
    file_list = os.listdir(data_dir)
    data_a.append(fnmatch.filter(file_list, '????a???_interp.txt')) 
    data_c.append(fnmatch.filter(file_list, '????c???_interp.txt'))
    data_f.append(fnmatch.filter(file_list, '????f???_interp.txt'))

    data_a = [fname for sublist in data_a for fname in sublist]
    data_c = [fname for sublist in data_c for fname in sublist]
    data_f = [fname for sublist in data_f for fname in sublist]

    assimilate_exp_a(data_a, 'lmla', opts)
    assimilate_exp_c(data_c, 'lmla', opts)
    assimilate_exp_f(data_f, 'lmla', opts)


def analyze_stokes(opts):
    """
    Analyzes all of the full Stokes based models
    """
    data_a, data_c, data_f = [], [], []
    data_dir = os.path.join(opts.input, 'stokes')
    file_list = os.listdir(data_dir)
    data_a.append(fnmatch.filter(file_list, '????a???_interp.txt')) 
    data_c.append(fnmatch.filter(file_list, '????c???_interp.txt'))
    data_f.append(fnmatch.filter(file_list, '????f???_interp.txt'))
    
    data_a = [fname for sublist in data_a for fname in sublist]
    data_c = [fname for sublist in data_c for fname in sublist]
    data_f = [fname for sublist in data_f for fname in sublist]

    assimilate_exp_a(data_a, 'stokes', opts)
    assimilate_exp_c(data_c, 'stokes', opts)
    assimilate_exp_f(data_f, 'stokes', opts)


def assimilate_ismip(opts):
    """
    Delegates interpolation to helper functions depending on the 
    implementation that the model that produced the data uses.
    """
    if not opts.lmla and not opts.stokes:
        print "Collecting all model runs...."
        analyze_lmla(opts)
        analyze_stokes(opts)
    elif opts.lmla and not opts.stokes:
        print "Collecting LMLA model runs...."
        analyze_lmla(opts)
    elif opts.stokes and not opts.lmla:
        print "Collecting Stokes model runs...."
        analyze_stokes(opts)
    elif opts.stokes and opts.lmla:
        print "Collecting all model runs...."
        analyze_lmla(opts)
        analyze_stokes(opts)


if __name__ == '__main__':
    """ Do everything when called from command line """
    opts = parse_args(sys.argv[1:])
   
    # Make sure that paths are okay
    if not os.path.exists(opts.input) or not os.path.isdir(opts.input):
        print("ERROR: Could not find input directory for analysis.")
        exit()

    assimilate_ismip(opts)


