from street import street
import os

# set the current dir as workplace
os.getcwd()

# the coarse (slave) and fine (master) data
slave_nc_file = './example/slave_tropomi_LA.nc'
master_nc_file = './example/master_tropomi_LA.nc'

# the name of vars pretaining to z,x,y
slave_nc_vars = ["values","lon","lat"]

# the master file has the same var
master_nc_vars = slave_nc_vars

# calling the object
street_obj = street(slave_nc_file,slave_nc_vars,master_nc_file,master_nc_vars,maxlag=5,n_bins=200)
# cal semivariograms and save their plots
street_obj.cal_semivar(do_plot = True, random_selection_n=4000)
# est spatial representation error for a large range of length scales
street_obj.error_estimator(do_plot=True)
# street_obj.length_scale and spatial_rep_err are the goal and are accessible from the object
# est spatial representation error for a given length scale (50 km)
street_obj.error_estimator(length_scale=50)
print(street_obj.spatial_rep_err_spc)

