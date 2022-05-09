from STREET import street

slave_nc_file = '/Users/asouri/Documents/Box_Aura_Ozone_study/STREET/STREET/example/slave_tropomi_LA.nc'
master_nc_file = '/Users/asouri/Documents/Box_Aura_Ozone_study/STREET/STREET/example/master_tropomi_LA.nc'

slave_nc_vars = ["values","lon","lat"]

master_nc_vars = slave_nc_vars

street_obj = street(slave_nc_file,slave_nc_vars,master_nc_file,master_nc_vars,maxlag=5,n_bins=200)
street_obj.cal_semivar(random_selection=4000)
street_obj.error_estimator(do_plot=True)

