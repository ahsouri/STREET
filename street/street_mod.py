# UTF-
# Estimate the Spatial Representation Error based on the principle of geostatistics
# Reference: Souri et al., 2022, Characterization of Errors in Satellite-based HCHO/NO2, 
#                                Tropospheric Column Ratios with Respect  to  Chemistry, 
#                                Column to PBL  Translation, Spatial Representation, and
#                                Retrieval Uncertainties 
# Amir Souri (ahsouri@cfa.harvard.edu;ahsouri@gmail.com)

import skgstat as skg
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import seaborn as sns
import os

class street(object):

    def __init__(self,slave_nc_file,slave_nc_vars,master_nc_file,master_nc_vars,semivar_model_slave = 1,
                 semivar_model_master = 1, maxlag = 5, minlag = 0.25, n_bins = 100):
        '''
           Initialize the street object with primary inputs
           ARGS:
                slave_nc_file(char): the file address with regards to the target field (slave)
                slave_nc_vars(list): the name of vars in the slave file about ,
                                  slave_nc_vars[0] = 2D Field
                                  slave_nc_vars[1] = Longitude
                                  slave_nc_vars[2] = Latitude
                                  example: slave_nc_vars[0] = "OMI_TNO2"
                                  currently the package doesn't support group-based nc files
                master_nc_file: the file address with regards to the reference field (master)
                master_nc_vars(list): the name of vars in the master file,
                                  master_nc_vars[0] = 2D Field
                                  master_nc_vars[1] = Longitude
                                  master_nc_vars[2] = Latitude
                                  example: master_nc_vars[0] = "TROPOMI_TNO2"
                                  currently the package doesn't support group-based nc files  
            semivar_model_slave or master (int): 1 -> Gaussian
                                 2 -> Spherical
                                 3 -> Exponential
                                 4 -> Stable Gauss
            maxlag (float): the maximum lag distance (default to 5 degree ~ 550 km)
            n_bins (int): number of bins in semivariogram 
        '''

        # read the slave file
        print("Reading the slave file: " + str(slave_nc_file))
        self.slave_field = self.read_netcdf(slave_nc_file,slave_nc_vars[0])
        self.slave_lon = self.read_netcdf(slave_nc_file,slave_nc_vars[1])
        self.slave_lat = self.read_netcdf(slave_nc_file,slave_nc_vars[2])

        # read the master file
        print("Reading the master file: " + str(master_nc_file))
        self.master_field = self.read_netcdf(master_nc_file,master_nc_vars[0])
        self.master_lon = self.read_netcdf(master_nc_file,master_nc_vars[1])
        self.master_lat = self.read_netcdf(master_nc_file,master_nc_vars[2])

        self.semivar_model_slave = semivar_model_slave
        self.semivar_model_master = semivar_model_master
        self.maxlag = maxlag
        self.minlag = minlag
        self.n_bins = n_bins
        print("the object is all set for other processes.")

    def read_netcdf(self,filename,var):
        ''' 
        Read nc format from a file without a group
        ARGS:
            filename (char): the name of file
            var (char): the target variable
        OUT:
            var (float)
        '''

        nc_f = filename
        nc_fid = Dataset(nc_f, 'r')
        var = nc_fid.variables[var][:]
        nc_fid.close()
        return np.squeeze(var)

    def cal_semivar(self,do_plot = False, plot_pngname = None, random_selection_n = None):
        ''' 
        Calculate the semivariogram
        ARGS (optional):
            do_plot (bool): whether you want to plot and save the semivariogram
            plot_pngname (char): the namefile of the plot (slave and master will be appended to it)
            random_selection_n (int): number of random samples from the field
        '''

        def cal_sem(x,y,z,semi_model,random_selection_n):
            # a mini function to cal semivariograms
            # mask bad data
            mask = np.isnan(z)
            z = z[~mask]
            x = x[~mask]
            y = y[~mask]

            # perform the random sampling if it's needed.
            if (random_selection_n is not None):
                ind = np.random.randint(np.size(z),size=random_selection_n)
                x = x[ind]
                y = y[ind]
                z = z[ind]

            # preparing x,y,z
            coords = np.empty((np.size(z),2))
            coords[:,0] = x.flatten()
            coords[:,1] = y.flatten()
            values = z.flatten()

            # choosing the semivariogram model
            if semi_model == 1:
               model = 'gaussian'
            elif semi_model== 2:
               model = 'spherical'
            elif semi_model == 3:
               model = 'exponential'
            elif semi_model == 4:
               model = 'stable'

            # running the skg (this can be memory intensive if the field is large)
            vario_obj = skg.Variogram(coords,values,n_lags=self.n_bins,estimator='matheron',
                                      model=model,maxlag=self.maxlag)
            # fitted model
            fitted_model = vario_obj.fitted_model
            # returning the goods
            return vario_obj,fitted_model

        def plotting(vario_obj,fname):
            # a mini function to plot the semivariogram
            fig, ax = plt.subplots()
            vario_obj.plot(axes = ax)
            fig.savefig(fname,dpi=300)
            plt.close()

        print("Building and modeling semivariogram for the slave")
        self.vario_obj_slave,self.fitted_model_slave = cal_sem(self.slave_lon,self.slave_lat,
                                                               self.slave_field,self.semivar_model_slave,
                                                               random_selection_n)
        print("Building and modeling semivariogram for the master")
        self.vario_obj_master,self.fitted_model_master = cal_sem(self.master_lon,self.master_lat,
                                                                self.master_field,self.semivar_model_master,
                                                                random_selection_n)

        if do_plot == True:
           if (not os.path.exists('plot_output')): os.makedirs('plot_output') 
           if (plot_pngname is not None):
              plotting(self.vario_obj_slave,  "plot_output/" + str(plot_pngname) + "_slave.png")
              plotting(self.vario_obj_master,  "plot_output/" + str(plot_pngname) + "_master.png")
           else:
              plotting(self.vario_obj_slave,"plot_output/semivariogram_slave.png")
              plotting(self.vario_obj_master,"plot_output/semivariogram_master.png")

    def error_estimator(self, do_plot = False, length_scale = None, deg2km = 110.0):
        ''' 
        Estimate the Spatial Representation Error based on Souri et al. 2022.
        ARGS (optional):
            do_plot (bool): whether you want to plot the error
            length_scale (float): unit: km, instead of estimating the error for a wide range
                                  the user can get a value for a given length scale
                                  if this is set to a value (or values), thee plotter
                                  is automatically off
            deg2km (float): degree to km conversion (default = 110 km)
        '''

        # the range of x for error estimation    
        x_value = np.arange(self.minlag,self.maxlag+0.1,0.02)

        # get their gamma values (i.e., variancc)
        var_slave = self.fitted_model_slave(x_value)
        var_master = self.fitted_model_master(x_value)

        # estimating the error for a given length scale
        if (length_scale is not None):
           self.spatial_rep_err_spc = 100.0*(1 - self.fitted_model_slave(length_scale/deg2km)\
               /self.fitted_model_master(length_scale/deg2km))
           return 1

        # the error
        self.spatial_rep_err = 1 - var_slave/var_master
        
        if do_plot:
           fig = plt.figure(figsize=(14, 8))
           sns.set_style("whitegrid")
           plt.plot(x_value*deg2km,self.spatial_rep_err*100.0,linewidth=6,color='purple')
           # Customize
           plt.title('Spatial Representation Error',fontsize=35)
           plt.ylabel('Loss of Spatial Information [%]', fontsize=30)
           plt.xlabel('Length Scale [km]', fontsize=30)
           plt.yticks(size=25)
           plt.xticks(size=25)
           plt.xlim(self.minlag*deg2km, self.maxlag*deg2km)
           plt.show()
        
        # saving the length scale
        self.length_scale = x_value









    

