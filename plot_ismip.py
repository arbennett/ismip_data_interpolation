"""
Python utility to plot a comparison of the original data to the interpolated
data that is the result of the interpolate_ismip.py script.  

@author arbennett
"""

from mpl_toolkits.mplot3d import axes3d
from itertools import product
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys 
import os

def parse_args(args):
    parser = argparse.ArgumentParser(description="ISMIP data comparison utility.",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            fromfile_prefix_chars='@')
    parser.add_argument('-o', '--original', help='The original ISMIP data file.')
    parser.add_argument('-i', '--interp', help='The interpolated ISMIP data file.')
    parser.add_argument('-v', '--var', help='The variable to plot')
    parser.add_argument('--list-vars', action='store_true', help='Show which variables are available to be plotted')

    return parser.parse_args(args)

def list_vars():
    print("Variable names for experiment a:")
    print("  vs_x, vs_y, vs_z, tb_xy, tb_yz, dp")
    print("")
    print("Variable names for experiment c:")
    print("  vs_x, vs_y, vs_z, vb_x, vb_y, tb_xy, tb_yz, dp")
    print("")
    print("Variable names for experiment f:")
    print("  zs, v_x, v_y, v_z")
    print("")
    exit()

def show_plot(orig_data, interp_data):
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.scatter(orig_data[0],orig_data[1],c=orig_data[2])
    ax2.scatter(interp_data[0],interp_data[1],c=interp_data[2])
    ax1.set_title("Original Data")
    ax2.set_title("Interpolated Data")
    plt.show()
    

def read_exp_a(orig, interp, var):
    """ Reads in experiment a data """
    x_o, y_o, vs_x_o, vs_y_o, vs_z_o, tb_xy_o, tb_yz_o, dp_o = np.loadtxt(orig, unpack=True)
    x_i, y_i, vs_x_i, vs_y_i, vs_z_i, tb_xy_i, tb_yz_i, dp_i = np.loadtxt(interp, unpack=True)
    vardict = {'vs_x' : (vs_x_o, vs_x_i), 'vs_y' : (vs_y_o, vs_y_i), 'vs_z' : (vs_z_o, vs_z_i), 
               'tb_xy' : (tb_xy_o, tb_xy_i), 'tb_yz' : (tb_yz_o, tb_yz_i), 'dp' : (dp_o, dp_i)}
    if var not in vardict.keys():
        print("ERROR: Could not find " + var + " in available data.  Use one of the following instead:")
        for v in vardict.keys(): print(v)
        print("Or use the --list-vars option to view all variables for each experiment.")
        exit()
    orig_data, interp_data = vardict[var]
    show_plot([x_o, y_o, orig_data], [x_i, y_i, interp_data])


def read_exp_c(orig, interp, var):
    """ Reads in experiment c data """
    x_o, y_o, vs_x_o, vs_y_o, vs_z_o, vb_x_o, vb_y_o, tb_xy_o, tb_yz_o, dp_o = np.loadtxt(orig, unpack=True)
    x_i, y_i, vs_x_i, vs_y_i, vs_z_i, vb_x_i, vb_y_i, tb_xy_i, tb_yz_i, dp_i = np.loadtxt(interp, unpack=True)
    vardict = {'vs_x' : (vs_x_o, vs_x_i), 'vs_y' : (vs_y_o, vs_y_i), 'vs_z' : (vs_z_o, vs_z_i),
               'vb_x' : (vb_x_o, vb_x_i), 'vb_y' : (vb_y_o, vb_y_i), 
               'tb_xy' : (tb_xy_o, tb_xy_i), 'tb_yz' : (tb_yz_o, tb_yz_i), 'dp' : (dp_o, dp_i)}
    if var not in vardict.keys():
        print("ERROR: Could not find " + var + " in available data.  Use one of the following instead:")
        for v in vardict.keys(): print(v)
        print("Or use the --list-vars option to view all variables for each experiment.")
        exit()
    orig_data, interp_data = vardict[var]
    show_plot([x_o, y_o, orig_data], [x_i, y_i, interp_data])


def read_exp_f(orig, interp, var):
    """ Reads in experiment f data """
    x_o, y_o, zs_o, v_x_o, v_y_o, v_z_o = np.loadtxt(orig, unpack=True)
    x_i, y_i, zs_i, v_x_i, v_y_i, v_z_i = np.loadtxt(interp, unpack=True)
    vardict = {'zs' : (zs_o, zs_i), 'v_x' : (v_x_o, v_x_i), 'v_y' : (v_y_o, v_y_i), 'v_z' : (v_z_o, v_z_i)}
    if var not in vardict.keys():
        print("ERROR: Could not find " + var + " in available data.  Use one of the following instead:")
        for v in vardict.keys(): print(v)
        print("Or use the --list-vars option to view all variables for each experiment.")
        exit()
    orig_data, interp_data = vardict[var]
    show_plot([x_o, y_o, orig_data], [x_i, y_i, interp_data])


def main(opts):
    """
    Main definition just delegates work to the correct areas.  Also does some 
    real basic error checking.  
    """
    if opts.list_vars:
        list_vars()
    if not os.path.exists(opts.original):
        print("ERROR: Could not find " + opts.original)
        exit()
    if not os.path.exists(opts.interp):
        print("ERROR: Could not find " + opts.interp)
        exit()
    data_functs = { 'a' : read_exp_a,
                    'c' : read_exp_c,
                    'f' : read_exp_f}
    orig_exp = os.path.basename(opts.original)[4]
    interp_exp = os.path.basename(opts.interp)[4]
    if not orig_exp == interp_exp and orig_exp in ['a','c','f']:
        print("ERROR: Data files not the same experiment!")
        exit()
    data_functs[orig_exp](opts.original, opts.interp, opts.var)


if __name__ == '__main__':
    opts = parse_args(sys.argv[1:])
    main(opts)


