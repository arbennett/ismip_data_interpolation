"""
Python utility to interpolate ISMIP benchmark data onto regular grids.
Current grid arrangement is compatible with CISM output from ISMIP tests.
Only covers experiments A, C, and F.  To expand to further experiments use
the interp_exp_* functions and the column definitions found in the ISMIP 
document found at http://homepages.ulb.ac.be/~fpattyn/ismip/ismiphom.pdf
on pages 25-27.

@author arbennett
"""

import os
import sys
import csv
import glob
import argparse
import numpy as np
from itertools import product
from scipy import interpolate

def parse_args(args):
    parser = argparse.ArgumentParser(description="ISMIP data interpolation utility.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            fromfile_prefix_chars='@')
    parser.add_argument('-i', '--input', default=os.path.join('.', 'ismip_all'),
                        help='The location of the ISMIP data directory.')
    parser.add_argument('-o', '--output', default=os.path.join('.','ismip_interp'), 
                        help='The directory to write the output files to.')
    parser.add_argument('--lmla', action='store_true',
                        help="Analyze LMLA models only.")
    parser.add_argument('--stokes', action='store_true',
                        help="Analyze Full Stokes models only.")
    return parser.parse_args(args)


def interp_exp_a(fname, out_dir):
    """ Used to interpolate data from experiment A. """
    print("  Beginning interpolation of " + fname)
    # The variables from the data
    print("    Reading data....")
    x, y, vs_x, vs_y, vs_z, tb_xy, tb_yz, dp = np.loadtxt(fname, unpack=True)
     
    #v_norm = np.sqrt(v_x**2 + v_y**2)
    res = 40  #int(fname.split(os.sep)[-1][5:8])
    
    # The given points
    x_pts = np.asarray(sorted(set(x)))
    y_pts = np.asarray(sorted(set(y)))
    points = np.transpose((x,y))

    # The points we want
    x_out = np.arange(0,1.0001,1.0/res)
    y_out = np.arange(0,1.0001,1.0/res)
    out_points = [[i,j] for i in x_out for j in y_out]
    x_out = np.transpose(out_points)[0]
    y_out = np.transpose(out_points)[1]

    # Interpolate each list separately
    print("    Interpolating data....")
    vs_x_i  = interpolate.griddata(points, vs_x , out_points)
    vs_y_i  = interpolate.griddata(points, vs_y , out_points)
    vs_z_i  = interpolate.griddata(points, vs_z , out_points)
    tb_xy_i = interpolate.griddata(points, tb_xy, out_points)
    tb_yz_i = interpolate.griddata(points, tb_yz, out_points)
    dp_i    = interpolate.griddata(points, dp   , out_points)

    out_file = os.path.join(out_dir, os.path.basename(fname).replace('.txt','_interp.txt'))
    print("    Writing data....")
    np.savetxt(out_file,np.transpose([x_out,y_out,vs_x_i,vs_y_i,vs_z_i,tb_xy_i,tb_yz_i,dp_i]))


def interp_exp_c(fname, out_dir):
    """ Used to interpolate data from experiment C. """
    print("  Beginning interpolation of " + fname)
    # The variables from the data
    print("    Reading data....")
    x, y, vs_x, vs_y, vs_z, vb_x, vb_y, tb_xy, tb_yz, dp = np.loadtxt(fname, unpack=True)
    
    #v_norm = np.sqrt(v_x**2 + v_y**2)
    res = 40  #int(fname.split(os.sep)[-1][5:8])
    
    # The given points
    x_pts = np.asarray(sorted(set(x)))
    y_pts = np.asarray(sorted(set(y)))
    points = (x,y)

    # The points we want
    x_out = np.arange(0,1.0001,1.0/res)
    y_out = np.arange(0,1.0001,1.0/res)
    out_points = [[i,j] for i in x_out for j in y_out]
    x_out = np.transpose(out_points)[0]
    y_out = np.transpose(out_points)[1]


    # Interpolate each list separately
    print("    Interpolating data....")
    vs_x_i  = interpolate.griddata(points, vs_x , out_points)
    vs_y_i  = interpolate.griddata(points, vs_y , out_points)
    vs_z_i  = interpolate.griddata(points, vs_z , out_points)
    vb_x_i  = interpolate.griddata(points, vb_x , out_points) 
    vb_y_i  = interpolate.griddata(points, vb_y , out_points)
    tb_xy_i = interpolate.griddata(points, tb_xy, out_points)
    tb_yz_i = interpolate.griddata(points, tb_yz, out_points)
    dp_i    = interpolate.griddata(points, dp   , out_points)
    
    out_file = os.path.join(out_dir,os.path.basename(fname).replace('.txt','_interp.txt'))
    print("    Writing data....")
    np.savetxt(out_file,np.transpose([x_out,y_out,vs_x_i,vs_y_i,vs_z_i,vb_x_i,vb_y_i,tb_xy_i,tb_yz_i,dp_i]))


def interp_exp_f(fname, out_dir):
    """ Used to interpolate data from experiment F. """
    print("  Beginning interpolation of " + fname)
    # The variables from the data
    print("    Reading data....")
    x, y, z_s, v_x, v_y, v_z = np.loadtxt(fname, unpack=True)

    #v_norm = np.sqrt(v_x**2 + v_y**2)
    res = 40 #int(fname.split(os.sep)[-1][5:8])
    
    # The given points
    x_pts = np.asarray(sorted(set(x)))
    y_pts = np.asarray(sorted(set(y)))
    points = (x,y)

    # The points we want
    x_out = np.arange(-50,50.0001,100.0/res)
    y_out = np.arange(-50,50.0001,100.0/res)
    out_points = [[i,j] for i in x_out for j in y_out]
    x_out = np.transpose(out_points)[0]
    y_out = np.transpose(out_points)[1]

    # Interpolate each list separately
    print("    Interpolating data....")
    z_s_i = interpolate.griddata(points, z_s, out_points)
    v_x_i = interpolate.griddata(points, v_x, out_points)
    v_y_i = interpolate.griddata(points, v_y, out_points)
    v_z_i = interpolate.griddata(points, v_z, out_points)
 
    out_file = os.path.join(out_dir, os.path.basename(fname).replace('.txt','_interp.txt'))
    print("    Writing data....")
    np.savetxt(out_file,np.transpose([x_out,y_out,z_s_i,v_x_i,v_y_i,v_z_i]))


def analyze_lmla(opts):
    """
    Analyzes all of the LMLA based models 
    """
    lmla = ['ahu1', 'ahu2', 'bds', 'fpa1', 'mbr', 'rhi2', 'tpa']
    exp_functs = {'a': interp_exp_a,
                  'c': interp_exp_c,
                  'f': interp_exp_f}
    out_dir = os.path.join(opts.output, 'lmla')
    os.mkdir(out_dir)
    for run in lmla:
        print("Analyzing " + run + "....")
        if not os.path.exists(os.path.join(opts.input, run)):
            print("  WARNING: Could not find directory for " + run + ".  Continuing....")
            continue
        files = glob.glob(os.path.join(opts.input, run, '*.txt'))
        if len(files) == 0:
            print("  WARNING: Could not find data for " + run + ".  Continuing....")
            continue
        for f in files:
            fname = os.path.basename(f).replace('.txt','')
            exp_type = fname[4]
            if exp_type not in ['a','c','f']:
                continue
            exp_functs[exp_type](f, os.path.join(opts.output,'lmla'))


def analyze_stokes(opts):
    """
    Analyzes all of the full Stokes based models
    """
    stokes = ['aas2', 'fpa2', 'jvj', 'mmr', 'oga', 'rhi1', 'rhi3', 'spr', 'ssu']
    exp_functs = {'a': interp_exp_a,
                  'c': interp_exp_c,
                  'f': interp_exp_f}
    out_dir = os.path.join(opts.output, 'stokes')
    os.mkdir(out_dir)
    for run in stokes:
        print("Analyzing " + run + "....")
        if not os.path.exists(os.path.join(opts.input, run)):
            print("  WARNING: Could not find directory for " + run + ".  Continuing....")
            continue
        files = glob.glob(os.path.join(opts.input, run, '*.txt'))
        if len(files) == 0:
            print("  WARNING: Could not find data for " + run + ".  Continuing....")
            continue
        for f in files:
            fname = os.path.basename(f).replace('.txt','')
            exp_type = fname[4]
            if exp_type not in ['a','c','f']:
                continue
            exp_functs[exp_type](f, os.path.join(opts.output,'stokes'))


def interpolate_ismip(opts):
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
    if os.path.exists(opts.output) and not os.path.isdir(opts.output):
        print("ERROR: Output path exists and is not a directory!")
        exit()
    if not os.path.exists(opts.output):
        os.mkdir(opts.output)
    if os.path.exists(os.path.join(opts.output, "stokes")) or os.path.exists(os.path.join(opts.output, "lmla")):
        print("ERROR: Output subdirectories already exist!")
        print("       Back up your work and run again.")
        exit()

    interpolate_ismip(opts)


