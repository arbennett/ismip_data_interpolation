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
import argparse
import numpy as np

from scipy import interpolate

def parse_args(args):
    parser = argparse.ArgumentParser(description="ISMIP data interpolation utility.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            fromfile_prefix_chars='@')
    parser.add_argument('-i', '--input', default='.', 
                        help='The location of the ISMIP data directory.')
    parser.add_argument('-o', '--output', default=os.path.join('.','ismip_interp'), 
                        help='The directory to write the output files to.')
    return parser.parse_args(args)


def interp_exp_a(fname):
    """ Used to interpolate data from experiment A. """
    x, y, vs_x, vs_y, vs_z, tb_xy, tb_yz, dp = np.loadtxt(fname, unpack=True)
    #vs_norm = np.sqrt(vs_x**2 + vs_y**2)
    return


def interp_exp_c(fname):
    """ Used to interpolate data from experiment C. """
    x, y, vs_x, vs_y, vs_z, vb_x, vb_y, vb_z, tb_xy, tb_yz, dp = np.loadtxt(fname, unpack=True)
    #vs_norm = np.sqrt(vs_x**2 + vs_y**2)
    #vb_norm = np.sqrt(vb_x**2 + vb_y**2)
    return


def interp_exp_f(fname):
    """ Used to interpolate data from experiment F. """
    print("Beginning interpolation of " + fname)
    # The variables from the data
    print("  Reading data....")
    x, y, z_s, v_x, v_y, v_z = np.loadtxt(fname, unpack=True)

    #v_norm = np.sqrt(v_x**2 + v_y**2)
    res = 100 #int(fname.split(os.sep)[-1][5:8])
    
    # The given points
    x_pts = np.asarray(sorted(set(x)))
    y_pts = np.asarray(sorted(set(x)))

    # The points we want
    x_out = np.arange(-50,50.0001,100.0/res)
    y_out = np.arange(-50,50.0001,100.0/res)
    xx, yy = np.meshgrid(x_out,y_out)
    xx = np.transpose(xx)
    yy = np.transpose(yy)
    
    # The interpolation variables
    zi_s = np.zeros(len(x_out)*len(y_out)) 
    vi_x = np.zeros(len(x_out)*len(y_out)) 
    vi_y = np.zeros(len(x_out)*len(y_out)) 
    vi_z = np.zeros(len(x_out)*len(y_out)) 
    interp_vars = [zi_s, vi_x, vi_y, vi_z]
   
    # Interpolate each list separately
    print("  Interpolating data....")
    z_s_i = interpolate.interp2d(x_pts, y_pts, z_s)
    v_x_i = interpolate.interp2d(x_pts, y_pts, v_x)
    v_y_i = interpolate.interp2d(x_pts, y_pts, v_y)
    v_z_i = interpolate.interp2d(x_pts, y_pts, v_z)
    data = [[j, i, z_s_i(i,j)[0], v_x_i(i,j)[0], v_y_i(i,j)[0], v_z_i(i,j)[0]] for i in x_out for j in y_out]
    
    with open(fname.replace('.txt', '') + "_interp.txt", 'w') as f:
        print("  Writing data....")
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(data)


def interpolate_ismip():
    """
    Delegates interpolation to helper functions depending on the 
    implementation that the model that produced the data uses.
    """
    pass


if __name__ == '__main__':
    """ Do everything when called from command line """
    opts = parse_args(sys.argv[1:])
   
    # Make sure that paths are okay
    if not os.path.exists(opts.input):
        print("ERROR: Could not find input directory for analysis.")
    if not os.path.exists(os.path.exists(opts.input, "ismip_all")) or not opts.input.split(os.sep)[-1] == 'ismip_all':
        print("ERROR: Could not find ismip_all directory for analysis.")
    if os.path.exists(opts.output) and not os.path.isdir(opts.output):
        print("ERROR: Output path exists and is not a directory!")
    if not os.path.exists(opts.output):
        os.path.mkdir(opts.output)

    interpolate_ismip()



