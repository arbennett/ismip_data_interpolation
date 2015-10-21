# ismip_data_interpolation
Provides scripts for interpolating ISMIP-HOM datasets onto new coordinates.

## Usage
All of the base data is found in the `ismip_all` directory.  To interpolate these datasets onto their respective CISM grids simply run `python interpolate_ismip.py`.  This will output a directory called `ismip_interp`.  If you wish to further distil these datasets into single files you can run `python assimilate_ismip.py`.  This will produce the `ismip_$EXP_$RES_$MODEL_combined.txt` files.  These are csv files where the data in each column is specified by the header.
