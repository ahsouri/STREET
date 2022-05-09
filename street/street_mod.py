# UTF-
# Estimate the Spatial Representation Error based on the principle of geostatistics
# Reference: Souri et al., 2022, Characterization of Errors in Satellite-based HCHO/NO2, 
#                                Tropospheric Column Ratios with Respect  to  Chemistry, 
#                                Column to PBL  Translation, Spatial Representation, and
#                                Retrieval Uncertainties 
# Amir Souri (ahsouri@cfa.harvard.edu;ahsouri@gmail.com)

class street(object):

    def __init__(self,slave_nc_file,slave_nc_vars,master_nc_file,master_nc_vars,semivar_model = 1,
                 maxlag = 5, n_bins = 100):
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
            semivar_model (int): 1 -> Gaussian
                                 3 -> Spherical
                                 4 -> Exponential
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

        self.semivar_model = semivar_model
        self.maxlag = maxlag
        self.n_bins = n_bins

    def read_netcdf(self,filename,var):
        ''' 
        Read nc format from a file without a group
        ARGS:
            filename (char): the name of file
            var (char): the target variable
        OUT:
            var (float)
        '''
        from netCDF4 import Dataset
        import numpy as np
        nc_f = filename
        nc_fid = Dataset(nc_f, 'r')
        var = nc_fid.variables[var][:]
        nc_fid.close()
        return np.squeeze(var)

    def cal_semivar(self,do_plot = False, random_selection = None):

        import skgstat as skg
        import numpy as np
        import matplotlib.pyplot as plt

        def cal_sem(x,y,z,random_selection):
            mask = np.isnan(z)
            z = z[~mask]
            x = x[~mask]
            y = y[~mask]

            if (random_selection is not None):
                ind = np.random.randint(np.size(z),size=random_selection)
                x = x[ind]
                y = y[ind]
                z = z[ind]

            coords = np.empty((np.size(z),2))
            coords[:,0] = x.flatten()
            coords[:,1] = y.flatten()
            values = z.flatten()

            if self.semivar_model == 1:
               model = 'gaussian'
            elif self.semivar_model == 2:
               model = 'spherical'
            elif self.semivar_model == 3:
               model = 'exponential'       

            vario_obj = skg.Variogram(coords,values,n_lags=self.n_bins,estimator='matheron',
                                      model=model,maxlag=self.maxlag)
            fitted_model = vario_obj.fitted_model
            return vario_obj,fitted_model

        def plotting(vario_obj,fname):
            fig, ax = plt.subplots()
            vario_obj.plot(axes = ax)
            fig.savefig(fname,dpi=300)
            plt.close()

        print("Building and modeling semivariogram for the slave")
        self.vario_obj_slave,self.fitted_model_slave = cal_sem(self.slave_lon,self.slave_lat,
                                                               self.slave_field,random_selection)
        print("Building and modeling semivariogram for the master")
        self.vario_obj_master,self.fitted_model_master = cal_sem(self.master_lon,self.master_lat,
                                                                self.master_field,random_selection)

        if do_plot == True:
           plotting(self.vario_obj_slave,"/Users/asouri/Documents/STREET/STREET/example/semivariogram_slave.png")
           plotting(self.vario_obj_master,"/Users/asouri/Documents/STREET/STREET/example/semivariogram_master.png")

        
    def error_estimator(self,do_plot = False,length_scale = None):

        import numpy as np
        import matplotlib.pyplot as plt  
        import seaborn as sns

        x_value = np.arange(0,self.maxlag+0.1,0.1)

        var_slave = self.fitted_model_slave(x_value)
        var_master = self.fitted_model_master(x_value)

        self.spatial_rep_err = 1 - var_slave/var_master
        
        if do_plot:
           fig = plt.figure(figsize=(14, 8))
           sns.set_style("whitegrid")
           plt.plot(x_value*110,self.spatial_rep_err*100.0,linewidth=6,color='purple')
           # Customize
           plt.title('Spatial Representation Error',fontsize=35)
           plt.ylabel('Loss of Spatial Information [%]', fontsize=30)
           plt.xlabel('Length Scale [km]', fontsize=30)
           plt.yticks(size=25)
           plt.xticks(size=25)
           plt.xlim(0, 500)
           plt.show()

        self.length_scale = x_value
        self.spatial_rep_err_spc = 1 - var_slave(length_scale*110)/var_master(length_scale*110)







    

