# Generates basic plots using basic colorscales etc. for a biogem experiment
# ... and gathers the plot in a PDF compiled with LaTeX with one big figure per page
# Also plots the difference if only 2 exps are provided
# The idea is to obtain for each run a summary that can be compared with other runs very quickly
# Very useful to compare results of 2 exps (dev and master) before merging 2 branches on github

# To run the script, adapt the paths and filenames and just 'python basicdiag_genie_toPDF.py'
# For now, all exps must be in the same directory

# The scripts largely relies on 3 functions:
# 1. geniemap.py: plots a 2D map of a genie output
# 2. genielat.py : plots a lat-depth output
# 3. genielev.py: plots a 2D map at every depth level

# Rk: In order to accomodate cgenie experiment names with dots, we create temporary ln -s 

# update Jan 19 2021 :: add seice fraction as a map
# update May 13 2020 :: computes the diff if a single run is provided, between 2 different saved time slices
# update Mar 24 2020 :: major change in the way files are read
#                       ... allowing more flexibility in terms of input files and saved variables                       
# update Mar 23 2020 ::  script should work for grid other than 36 x 36 x 16                          
#                        possibility to change the cartopy map projection
#                        improved format of the user params
#                        improved communication between python and latex
# update Mar 6 2020 :: includes sedgem/omensed basic output
# update Dec 1 2019 :: plt.step corrected: using 'post', and appending the last value to each array plotted
# update Dec 1 2019 :: plots requested year (e.g., 9999.5) instead of the last time step
# ... also used to make sure that the run reached the expected duration
# update Nov 21 2019 :: if only 2 exps provided, also plots the difference and generates a 3rd PDF

# TODO
# avoid doing diff if grids are different and...
# allow diff for different grids except that does not plot diff of different vertical levels

import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as mticker
import netCDF4
from netCDF4 import Dataset
import numpy as np
import matplotlib.colors as colors
import os
import string
import warnings
import netCDF4
from netCDF4 import Dataset
from matplotlib.colors import BoundaryNorm
from sklearn.linear_model import LinearRegression
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

########################## USER_DEFINED OPTIONS ##########################
indir ='EXAMPLE.input'
exps = ['wardetal.2018.ECOGEM.SPIN']
times2plot = [9999.5]

# MAP PROJECTIONS [https://scitools.org.uk/cartopy/docs/latest/crs/projections.html]
projdata = ccrs.LambertCylindrical() # Andy's style
#projdata = ccrs.EqualEarth() # The trendy one
#projdata = ccrs.RotatedPole() # The funny one
#projdata = ccrs.EckertIV() # Chris' style
#projdata = ccrs.PlateCarree() # Basic and boring

# MUFFIN CONFIG
do_biogem = 'y'
do_sedgem = 'n' # are we plotting sedgem outputs?
##########################################################################

# additional options (should not be changed in general)
do_diff = 'y'
save_fig = 'y'
create_pdf_summary = 'y'
plot_stepped_outline = 'y'


# nothing to change below this line


if len(exps) > 2 or len(times2plot) > 2 or (len(exps) > 2 and len(times2plot) > 2):
    do_diff = 'n'

# using standard python frontend instead of the web browser
mpl.use('TkAgg')
 # ignore futurewarnings
warnings.filterwarnings("ignore", category=FutureWarning)

# graphical options
plot_font_size = 13
label_size = plot_font_size/2.
plot_label_size = plot_font_size+3
mpl.rc('font', size=plot_font_size+7)
# for the figure
figXsize = 6
figYsize = 4.5
# for the annotations on the maps
xlabel = 0
ylabel = 1.05
# data projection system
data_crs = ccrs.PlateCarree()
# functions
functionstoload = ['stepped_coastline_cGENIE','stepped_outline_cGENIE','light_grid','geniemap','genielat','fakealpha','custom_colormaps', 'custom_chars','genielev', 'dopdf']
for function2load in functionstoload:
    string2execute = 'source/' + function2load + '.py'
    exec(open(string2execute).read()) # python 3
diffcmap = light_centered
difflower = 'darkblue'
diffupper = 'darkred'

# Here we provide the list of variables to read
# right column = name of the variable in the ncfile
# left column = name that you want to give to this variable in the workspace
# typically, I add 'sed' to sedgem variables to avoid redundance with the oceanic component
# If a field is not present in the nc file, it will just be skipped, no worries
dict_biogem3d = {
'time':'time',
'grid_topo':'grid_topo',
'grid_mask_3d':'grid_mask_3d',
'lon':'lon',
'lat':'lat',
'lon_edges':'lon_edges',
'lat_edges':'lat_edges',
'grid_area_ocn':'grid_area',
'ocn_temp':'ocn_temp',
'ocn_sal':'ocn_sal',
'phys_ocn_rho':'phys_ocn_rho',
'ocn_O2':'ocn_O2',
'ocn_H2S':'ocn_H2S',
'ocn_PO4':'ocn_PO4',
'ocn_DIC_13C':'ocn_DIC_13C',
'misc_col_Dage':'misc_col_Dage',
}

dict_biogem2d = {
'atm_temp':'atm_temp', # basic ocean state
'phys_wspeed':'phys_wspeed',
'phys_opsi':'phys_opsi',
'phys_psi':'phys_psi',
'phys_seaice':'phys_seaice',
'phys_cost':'phys_cost',
'grid_area_atm':'grid_area',
'bio_export_POC':'bio_export_POC', # export PROD
'bio_export_POC':'bio_export_POC',
'bio_diag_k_temp':'bio_diag_k_temp',
'bio_diag_k_light':'bio_diag_k_light',
'bio_diag_k_PO4':'bio_diag_k_PO4',
'lon_psi':'lon_psi', # grids
'lon_psi_edges':'lon_psi_edges',
'lat_psi':'lat_psi',
'lat_psi_edges':'lat_psi_edges',
'lat_moc':'lat_moc',
'lat_moc_edges':'lat_moc_edges',
'zt_moc':'zt_moc',
'zt':'zt',
'zt_moc_edges':'zt_moc_edges',
'zt_edges':'zt_edges',
}

dict_sedgem2d = {
'sed_lon':'lon',
'sed_lat':'lat',
'sed_lon_edges':'lon_edges',
'sed_lat_edges':'lat_edges',
'sed_grid_topo':'sed_grid_topo',
'sed_grid_mask':'sed_grid_mask',
'sed_ocn_temp':'ocn_temp', # overlying ocean properties
'sed_ocn_sal':'ocn_sal',
'sed_ocn_DIC':'ocn_DIC',
'sed_ocn_DIC_13C':'ocn_DIC_13C',
'sed_ocn_PO4':'ocn_PO4',
'sed_ocn_O2':'ocn_O2',
'sed_ocn_ALK':'ocn_ALK',
'sed_ocn_H2S':'ocn_H2S',
'sedocn_fnet_DIC_13C':'sedocn_fnet_DIC_13C', # benthic interface exchange
'sedocn_fnet_PO4':'sedocn_fnet_PO4',
'fsed_POC':'fsed_POC', # sedimentary fluxes
'fsed_POC_frac2':'fsed_POC_frac2',
'fsed_POC_13C':'fsed_POC_13C',
'fburial_det':'fburial_det', # sediment burial flux
'fburial_POC':'fburial_POC',
'fburial_POC_13C':'fburial_POC_13C',
'OMEN_wtpct_top':'OMEN_wtpct_top', # OMENSED
'OMEN_wtpct_bot':'OMEN_wtpct_bot'
}

# function to read nc fields only if they are present in the nc file
def readncfile(inpath, indic, mydata):
    my_data=mydata
    f = Dataset(inpath)
    for x in indic:
        varname = x
        var2read = indic[x]
        if var2read in f.variables:
            my_data[varname] = f.variables[var2read][:]
    f.close()
    return(my_data)

def dovar(varname): # checking that var exists
    truefalse = varname in globals()
    return(truefalse)

globalcount = 0
expcount = 0
for exp in exps:

    timecount = 0
    for time2plot in times2plot:

        plottedT = 'yr' + str(time2plot)

        savedfiles = []
        diffsavedfiles = []
        filecount = 0
        difffilecount = 0

        if do_diff == 'y':
            if globalcount ==0:
                exp0 = exp
                plottedT0 = plottedT
            elif globalcount ==1:
                exp1 = exp
                plottedT1 = plottedT

        ########################################################################
        #                               READING DATA                           #
        ########################################################################

        # Here we loop over ncfiles, loading existing variables and skipping missing ones 
        # And gathering the data into a python ensemble called 'my_data'
        my_data={}
        if do_biogem == 'y':
            my_data = readncfile(indir + '/' + exp + '/biogem/fields_biogem_3d.nc', dict_biogem3d, my_data)
            my_data = readncfile(indir + '/' + exp + '/biogem/fields_biogem_2d.nc', dict_biogem2d, my_data)
        if do_sedgem == 'y':
            my_data = readncfile(indir + '/' + exp + '/sedgem/fields_sedgem_2d.nc', dict_sedgem2d, my_data)

        # Here we read every field in 'my_data' and send them to the workspace
        for varname in my_data:
            exec(varname + " = my_data['" + varname + "']")

        ########################################################################
        #                                  BIOGEM                              #
        ########################################################################

        if do_biogem == 'y':

            # ===================== LOADING DATA =====================
     
            # checking grids: no difference if grids are not the same (nor now...)
            #if do_diff == 'y':
            #    if globalcount ==0:
            #        zt_0 = zt
            #        lon_0 = lon
            #        lat_0 = lat
            #    elif globalcount ==1:
            #        zt_1 = zt
            #        lon_1 = lon
            #        lat_1 = lat
            #        if (zt_1 != zt_0).all or (lon_1 != lon_0).all or (lat_1 != lat_0).all:
            #            do_diff = 'n'
            #            print('Grids are different, disabling do_diff')

            # extracting the grid
            nz = np.shape(ocn_temp)[1]
            nlon = np.shape(ocn_temp)[2]
            nlat = np.shape(ocn_temp)[3]

            # extracting time slice of interest
            T = np.argwhere(time == time2plot)[0][0] # index of the time slice to plot, as a float (not an array)

            # ===================== CALCULATIONS ===================== 

            # zonal temp
            lat_ocn_temp = np.ma.mean(np.squeeze(ocn_temp[T,:,:,:]),axis=2)
            # zonal salt
            lat_ocn_sal = np.ma.mean(np.squeeze(ocn_sal[T,:,:,:]),axis=2)

            # SST
            SST = ocn_temp[T,0,:,:] # (36, 36)
            latSST = np.ma.mean(np.squeeze(SST),axis=1)
            ilat_tropics = np.where(abs(lat)<=30)
            SST_tropics = np.squeeze(SST[ilat_tropics,:]) # checked: extracts the tropics
            grid_area_ocn_tropics = np.squeeze(grid_area_ocn[ilat_tropics,:]) # checked: grid_area_ocn mashed over land
            SST_tropics_avg = np.round(np.sum(np.multiply(SST_tropics,grid_area_ocn_tropics))/np.sum(grid_area_ocn_tropics),1)

            # SSS
            SSS = ocn_sal[T,0,:,:] # (36, 36)

            if dovar('phys_ocn_rho'):
                phys_ocn_rho_surface = phys_ocn_rho[T,0,:,:]

            # building fake land-sea mask - required by stepped_coastline
            if np.shape(np.argwhere((grid_topo.mask==True)))[0] > 0: # if some land points
                landsea_mask = np.where(grid_topo.mask==True,-1,1)
                land = np.ma.masked_where(landsea_mask==1,landsea_mask)
                landflag = 1
            else: # if waterworld, needed to avoid script to crash
                landsea_mask = np.full(np.shape(grid_topo),-1)
                land = np.full(np.shape(grid_topo), np.nan)
                landflag = 0
                plot_stepped_outline = 'n'

            # ocean area
            oceanareatot = np.sum(grid_area_ocn)

            # SAT
            SAT = atm_temp[T,:,:]
            SAT_avg = np.round(np.sum(np.multiply(SAT,grid_area_atm))/np.sum(grid_area_atm),3)
            if do_diff == 'y':
                if globalcount ==0:
                    sat_0 = SAT_avg
                elif globalcount ==1:
                    sat_1 = SAT_avg
                    SAT_avg_diff = sat_1 - sat_0

            # phys_seaice
            lat_phys_seaice = np.ma.mean(phys_seaice[T,:,:],axis=1)
     
            # SI fraction
            SIfrac = phys_seaice[T,:,:]
            SIfrac50 = np.ma.masked_where(SIfrac < 50,  SIfrac)
            SATstr = 'Global SAT = ' + str(SAT_avg) + ' ' + degree_sign + 'C'

            # deep O2
            O2 = ocn_O2[T,:,:,:]*1E6 # (16, 36, 36); umol L-1)
            deepO2 = np.full(np.array([36, 36]),np.nan)
            # deepest level?
            for i in np.arange(0,36):
                for j in np.arange(0,36):
                    ix = np.squeeze(np.ma.where(grid_mask_3d[:,i,j].mask == False))
                    if ix.any(): # if land point, leave np.nan in the array
                        ideepest = ix[-1]
                        #print(ideepest)
                        deepO2[i,j] = O2[ideepest,i,j]
                    else:
                        deepO2[i,j] = O2[-1,i,j]
            deepO2 = np.ma.masked_where(grid_topo.mask==True,deepO2)
            # deepO2 on shelves
            shallow_threshold = 400 # m
            deepO2_shallowpfonly = np.ma.masked_where(grid_topo > shallow_threshold,deepO2)
            pfmask = np.full(np.array([nlon, nlat]),1)
            pfmask[grid_topo > shallow_threshold] = -1
            # area of seafloor anoxia?
            anoxia_threshold = 0
            seafloor_anoxic_area = np.ma.masked_where(deepO2 > anoxia_threshold,grid_area_ocn)
            seafloor_anoxia_abs = np.ma.sum(seafloor_anoxic_area)
            seafloor_anoxia_per = seafloor_anoxia_abs / np.ma.sum(grid_area_ocn)
            if do_diff == 'y':
                if globalcount ==0:
                    seafloor_anoxia_abs_0 = seafloor_anoxia_abs
                    seafloor_anoxia_per_0 = seafloor_anoxia_per
                elif globalcount ==1:
                    seafloor_anoxia_abs_1 = seafloor_anoxia_abs
                    seafloor_anoxia_per_1 = seafloor_anoxia_per
                    seafloor_anoxia_abs_diff = seafloor_anoxia_abs_1 - seafloor_anoxia_abs_0
                    seafloor_anoxia_per_diff = seafloor_anoxia_per_1 - seafloor_anoxia_per_0

            # ocn_O2
            lat_ocn_O2 = np.ma.mean(O2,axis=2)

            # ocn_H2S
            H2S = ocn_H2S[T,:,:,:]*1E6 # (16, 36, 36); umol L-1)
            ma_H2S = np.ma.masked_where(np.squeeze(H2S) == 0, np.squeeze(H2S))
            lat_ocn_H2S = np.ma.mean(H2S,axis=2)
            ma_lat_ocn_H2S = np.ma.masked_where(lat_ocn_H2S == 0,lat_ocn_H2S)

            if dovar('misc_col_Dage'):
                lat_misc_col_Dage = np.ma.mean(misc_col_Dage[T,:,:,:], axis=2)
           
            lat_ocn_DIC_13C = np.ma.mean(ocn_DIC_13C[T,:,:,:], axis=2) 


            # ===================== PLOTTING =====================

            if dovar('grid_topo'):
                # %%%%%%%%%%%% grid_topo %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0,5.5+1E-6,0.5)
                ticklevs = levs
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                cbartitle = 'Bathymetry (km. b.s.l.)'
                filename = exp + '_' + plottedT + '_grid_topo.png'
                # --- figure ---
                geniemap(lon_edges, lat_edges, grid_topo*1E-3, cmap, levs, ticklevs, 'max', lower, upper, projdata ,cbartitle, filename, 'n', 'none', 'none', 'none', 0, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        grid_topo_0 = grid_topo*1E-3
                    elif globalcount ==1:
                        grid_topo_1 = grid_topo*1E-3
                        diff = grid_topo_1 - grid_topo_0
                        difflevs = np.arange(-5.,5.+1E-3,0.5)
                        diffticklevs = np.arange(-5.,5.+1E-3,1)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_grid_topo.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'n', 'none', 'none', 'none', 0, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('phys_wspeed'):
                # %%%%%%%%%%%% phys_wspeed %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0,8+1E-6,0.5)
                clevs = np.arange(0,8+1E-6,1)
                ticklevs = np.arange(0,8+1E-6,2)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                cbartitle = 'Wind speed (m s$^{-1}$)'
                filename = exp + '_' + plottedT + '_phys_wspeed.png'
                # --- figure ---
                geniemap(lon_edges, lat_edges, phys_wspeed[T,:,:], cmap, levs, ticklevs, 'max', lower, upper, projdata ,cbartitle, filename, 'y', lon, lat, clevs, 0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        phys_wspeed_0 = phys_wspeed[T,:,:]
                    elif globalcount ==1:
                        phys_wspeed_1 = phys_wspeed[T,:,:]
                        diff = phys_wspeed_1 - phys_wspeed_0
                        difflevs = np.arange(-3,3+1E-3,0.25)
                        diffclevs = np.arange(-3,3+1E-3,0.5)
                        diffticklevs = np.arange(-3,3+1E-3,1)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_phys_wspeed.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'y', lon, lat, diffclevs, 0.75, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('phys_ocn_rho'):
                # %%%%%%%%%%%% SSS %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(1020,1030+1E-9,0.5) # 1E-3 is just to include the upper bound
                ticklevs = np.arange(1020,1030+1E-9,2)
                clevs = np.arange(1020,1030+1E-9,0.5)
                extend = 'both'
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap2
                cbartitle = 'rho'
                filename = exp + '_' + plottedT + '_rho.png'
                # --- figure ---
                geniemap(lon_edges, lat_edges, phys_ocn_rho_surface, cmap, levs, ticklevs, extend, lower, upper, projdata ,cbartitle, filename, 'y', lon, lat, clevs, 0.5, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        phys_ocn_rho_surface_0 = phys_ocn_rho_surface
                    elif globalcount ==1:
                        phys_ocn_rho_surface_1 = phys_ocn_rho_surface
                        diff = phys_ocn_rho_surface_1 - phys_ocn_rho_surface_0
                        difflevs = np.arange(-2.,2.+1E-9,0.25)
                        diffclevs = np.arange(-2.,2.+1E-9,0.25)
                        diffticklevs = np.arange(-2.,2.+1E-9,0.5)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_rho.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'n', lon, lat, diffclevs, 0.5, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)


            if dovar('SSS'):
                # %%%%%%%%%%%% SSS %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(31,37+1E-9,0.25) # 1E-3 is just to include the upper bound
                ticklevs = np.arange(31,37+1E-9,1)
                clevs = np.arange(31,37+1E-9,1)
                extend = 'both'
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap2
                cbartitle = 'SSS (PSU)'
                filename = exp + '_' + plottedT + '_SSS.png'
                # --- figure ---
                geniemap(lon_edges, lat_edges, SSS, cmap, levs, ticklevs, extend, lower, upper, projdata ,cbartitle, filename, 'y', lon, lat, clevs, 0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        SSS_0 = SSS
                    elif globalcount ==1:
                        SSS_1 = SSS
                        diff = SSS_1 - SSS_0
                        difflevs = np.arange(-1.,1.+1E-6,0.1)
                        diffclevs = np.arange(-1.,1.+1E-6,0.1)
                        diffticklevs = np.arange(-1.,1.+1E-6,0.5)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_SSS.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'y', lon, lat, diffclevs, 1, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('lat_ocn_sal'):
                # %%%%%%%%%%% lat-depth sal profile %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(34.6,35.2+1E-9,0.025)
                ticklevs = np.arange(34.6,35.2+1E-9,0.1)
                clevs = np.arange(34.6,35.2+1E-9,0.025)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                extend = 'both'
                cbartitle = 'Salinity (psu)'
                filename = exp + '_' + plottedT + '_lat_ocn_sal.png'
                # --- figure ---
                genielat(lat_edges,-1*zt_edges/1000.,lat_ocn_sal, cmap, levs, ticklevs, 'both', lower, upper, cbartitle, filename, 'y', lat, -1*zt/1000., levs, 0.5)
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        lat_ocn_sal_0 = lat_ocn_sal
                    elif globalcount ==1:
                        lat_ocn_sal_1 = lat_ocn_sal
                        diff = lat_ocn_sal_1 - lat_ocn_sal_0
                        difflevs = np.arange(-0.4,0.4+1E-3,0.02)
                        diffticklevs = np.arange(-0.4,0.4+1E-3,0.1)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_lat_ocn_sal.png'
                        genielat(lat_edges,-1*zt_edges/1000.,diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, cbartitle, difffilename, 'y', lat, -1*zt/1000., difflevs, 0.5)
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('SST'):
                # %%%%%%%%%%%% SST with sea ice %%%%%%%%%%%%
                 # not using standard functions because contour labels
                # --- parameters ---
                levs = np.arange(-2,36+1E-9,2) # 1E-3 is just to include the upper bound
                clevsneg = np.array([-2])
                clevspos = np.arange(2,50+1E-9,2)
                clevszero = np.array([0])
                ticklevs = levs
                extend = 'max'
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap2
                cbartitle = 'SST' + degree_sign + ' C'
                filename = exp + '_' + plottedT + '_SST_SIfrac50perc.png'
                # --- figure ---
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                fig =  plt.figure(figsize=(figXsize, figYsize))
                ax = fig.add_subplot(111,projection=projdata)
                ax = fig.gca(); ax.set(aspect=2); ax.set_aspect('auto')
                cf = plt.pcolormesh(lon_edges, lat_edges, SST, transform=data_crs, cmap=cmap, norm=norm)
                cneg = plt.contour(lon, lat, SST, clevsneg, cmap = None, colors='k',linewidths=0.5, linestyles='dashed',transform=data_crs)
                cpos = plt.contour(lon, lat, SST, clevspos, cmap = None, colors='k',linewidths=0.5,transform=data_crs)
                czero = plt.contour(lon, lat,SST, clevszero, cmap = None, colors='k',linewidths=0.75,transform=data_crs)
                plt.clabel(cneg,inline=1,inline_spacing=5, fontsize=label_size,fmt='%1.0f',colors='k')
                plt.clabel(cpos,inline=1,inline_spacing=5,fontsize=label_size,fmt='%1.0f',colors='k')
                plt.clabel(czero,inline=1,inline_spacing=5,fontsize=label_size,fmt='%1.0f',colors='k')
                cf2 = plt.pcolormesh(lon_edges, lat_edges, SIfrac50, cmap=lightgreycmap, transform=ccrs.PlateCarree())
                plt.pcolormesh(lon_edges, lat_edges, land, cmap=whitecmap, transform=ccrs.PlateCarree())
                if landflag == 1:
                    stepped_coastline_cGENIE(lon_edges,lat_edges, data_crs, grid_area_ocn,1.25)
                gl = ax.gridlines(crs=data_crs, draw_labels=False,
                              linewidth=0.5, color='k', alpha=0.2, linestyle='-', zorder=999)
                gl.xlabels_top = False; gl.ylabels_left = False; gl.xlines = True
                gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180]); gl.xformatter = LONGITUDE_FORMATTER
                gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90]); gl.yformatter = LATITUDE_FORMATTER
                gl.xlabel_style = {'size': plot_font_size, 'color': 'gray'}; gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
                cb = fig.colorbar(cf, orientation='horizontal',extend=extend,ticks=ticklevs)
                cf.cmap.set_under(lower); cf.cmap.set_over(upper)
                cb.ax.tick_params(labelsize=plot_font_size)
                cb.ax.set_title(cbartitle, weight='normal', fontsize=plot_font_size)
                plt.tight_layout()
                filename = exp + '_' + plottedT + '_SST_SIfrac50perc.png'
                if save_fig == 'y':
                    plt.savefig(filename,format='png', dpi=450)
                    filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                    os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        SST_0 = SST
                    elif globalcount ==1:
                        SST_1 = SST
                        diff = SST_1 - SST_0
                        difflevs = np.arange(-7.5,7.5+1E-3,0.5)
                        diffclevs = np.arange(-7.5,7.5+1E-3,0.5)
                        diffticklevs = np.arange(-7.5,7.5+1E-3,2.5)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_SST.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'y', lon, lat, diffclevs, 1., 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('lat_ocn_temp'):
                # %%%%%%%%%%% lat-depth temp profile %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0,30+1E-9,1)
                ticklevs = np.arange(0,30+1E-9,4)
                clevs = np.arange(0,30+1E-9,1)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                extend = 'max'
                cbartitle = 'Temperature (' + degree_sign + 'C)'
                filename = exp + '_' + plottedT + '_lat_ocn_temp.png'
                # --- figure ---
                genielat(lat_edges,-1*zt_edges/1000.,lat_ocn_temp, cmap, levs, ticklevs, 'both', lower, upper, cbartitle, filename, 'y', lat, -1*zt/1000., levs, 0.5)
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        lat_ocn_temp_0 = lat_ocn_temp
                    elif globalcount ==1:
                        lat_ocn_temp_1 = lat_ocn_temp
                        diff = lat_ocn_temp_1 - lat_ocn_temp_0
                        difflevs = np.arange(-5,5+1E-3,0.5)
                        diffticklevs = np.arange(-5,5+1E-3,2.5)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_lat_ocn_temp.png'
                        genielat(lat_edges,-1*zt_edges/1000.,diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, cbartitle, difffilename, 'y', lat, -1*zt/1000., difflevs, 0.5)
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('latSST'):
                # %%%%%%%%%%%% zonal SST profile %%%%%%%%%%%%
                # --- figure ---
                fig =  plt.figure(figsize=(figXsize, figYsize))
                ax = fig.add_subplot(111)
                ax = fig.gca()
                cf = plt.step(lat_edges, np.append(latSST,latSST[-1]), '-k', where='post')
                plt.ylim(-90, 90)
                plt.ylim(-3,38)
                major_yticks = np.arange(0, 35+1E-9, 5)
                major_xticks = np.arange(-90, 90+1E-9, 30)
                ax.set_xticks(major_xticks)
                ax.set_yticks(major_yticks)
                ax.tick_params(labelsize = plot_font_size)
                ax.grid(which='both')
                if save_fig == 'y':
                    filename = exp + '_' + plottedT + '_latSST.png' 
                    plt.savefig(filename,format='png')
                    filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                    os.system('ln -s ' + filename + ' ' + lnfile)
                    savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        latSST_0 = latSST
                    elif globalcount ==1:
                        latSST_1 = latSST
                        diff = latSST_1 - latSST_0
                        fig =  plt.figure(figsize=(figXsize, figYsize))
                        ax = fig.add_subplot(111)
                        ax = fig.gca()
                        plt.step(lat_edges, np.append(latSST_0,latSST_0[-1]), '-r', where='post')
                        plt.step(lat_edges, np.append(latSST_1,latSST_1[-1]), '-b', where='post')
                        plt.ylim(-90, 90)
                        plt.ylim(-5,40)
                        major_yticks = np.arange(-5, 40+1E-9, 5)
                        major_xticks = np.arange(-90, 90+1E-9, 30)
                        ax.set_xticks(major_xticks)
                        ax.set_yticks(major_yticks)
                        ax.tick_params(labelsize = plot_font_size)
                        ax.grid(which='both')
                        if save_fig == 'y':
                            difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_latSST.png'
                            plt.savefig(difffilename,format='png')
                            difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                            os.system('ln -s ' + difffilename + ' ' + difflnfile)
                            diffsavedfiles.append(difffilename)

            if dovar('lat_phys_seaice'):
                # %%%%%%%%%%%% phys_seaice %%%%%%%%%%%%
                # --- figure ---
                fig =  plt.figure(figsize=(figXsize, figYsize))
                ax = fig.add_subplot(111)
                ax = fig.gca(); ax.set(aspect=0.5)
                cf = plt.step(lat_edges, np.append(lat_phys_seaice,lat_phys_seaice[-1]), '-k', where='post')
                plt.ylim(-90, 90)
                plt.ylim(-0,110)
                major_yticks = np.arange(0, 100+1E-9, 20)
                major_xticks = np.arange(-90, 90+1E-9, 30)
                ax.set_xticks(major_xticks)
                ax.set_yticks(major_yticks)
                ax.tick_params(labelsize = plot_font_size)
                ax.grid(which='both')
                if save_fig == 'y':
                    filename = exp + '_' + plottedT + '_phys_seaice_fraction.png'
                    plt.savefig(filename,format='png')
                    filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                    os.system('ln -s ' + filename + ' ' + lnfile)
                    savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        lat_phys_seaice_0 = lat_phys_seaice
                    elif globalcount ==1:
                        lat_phys_seaice_1 = lat_phys_seaice
                        diff = lat_phys_seaice_1 - lat_phys_seaice_0
                        fig =  plt.figure(figsize=(figXsize, figYsize))
                        ax = fig.add_subplot(111)
                        ax = fig.gca(); ax.set(aspect=0.5)
                        plt.step(lat_edges, np.append(lat_phys_seaice_0,lat_phys_seaice_0[-1]), '-r', where='post')
                        plt.step(lat_edges, np.append(lat_phys_seaice_1,lat_phys_seaice_1[-1]), '-b', where='post')
                        plt.ylim(-90, 90)
                        plt.ylim(0,110)
                        major_yticks = np.arange(0, 100+1E-9, 20)
                        major_xticks = np.arange(-90, 90+1E-9, 30)
                        ax.set_xticks(major_xticks)
                        ax.set_yticks(major_yticks)
                        ax.tick_params(labelsize = plot_font_size)
                        ax.grid(which='both')
                        if save_fig == 'y':
                            difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_phys_seaice_fraction.png'
                            plt.savefig(difffilename,format='png', dpi=450)
                            difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                            os.system('ln -s ' + difffilename + ' ' + difflnfile)
                            diffsavedfiles.append(difffilename)

            if dovar('phys_opsi'):
                # %%%%%%%%%%%% MOC %%%%%%%%%%%%
                # --- parameters ---
                levs = np.array([-30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30])
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                extend = 'both'
                cbartitle = 'Overturning streamfunction'
                filename = exp + '_' + plottedT + '_phys_opsi.png'
                clevs = np.arange(-35,35+1E-9,5)
                # --- figure ---
                genielat(lat_moc_edges,-1*zt_moc_edges/1000.,phys_opsi[T,:,:], cmap, levs, clevs, 'both', lower, upper, cbartitle, filename, 'y', lat_moc, -1*zt_moc/1000., levs, 0.75)
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        phys_opsi_0 = phys_opsi[T,:,:]
                    elif globalcount ==1:
                        phys_opsi_1 = phys_opsi[T,:,:]
                        diff = phys_opsi_1 - phys_opsi_0
                        difflevs = np.arange(-30,30+1E-3,2.5)
                        #difflevs = np.arange(-10,10+1E-3,1)
                        diffticklevs = np.arange(-30,30+1E-3,5)
                        #diffticklevs = np.arange(-10,10+1E-3,5)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_phys_opsi.png'
                        genielat(lat_moc_edges,-1*zt_moc_edges/1000.,diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, cbartitle, difffilename, 'y', lat_moc, -1*zt_moc/1000., difflevs, 0.5)
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('phys_seaice'):
                # %%%%%%%%%%%% sea-ice cover percent %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0,100+1E-9,5.)
                ticklevs = np.arange(0,100+1E-9,10)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                extend = 'neither'
                cbartitle = 'Sea-ice cover (%)'
                filename = exp + '_' + plottedT + '_phys_seaice.png'
                # --- figure ---
                geniemap(lon_edges, lat_edges, phys_seaice[T,:,:], cmap, levs, ticklevs, 'both', lower, upper, projdata ,cbartitle, filename, 'n', 'none', 'none', 'none', 0., 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        phys_seaice_0 = phys_seaice[T,:,:]
                    elif globalcount ==1:
                        phys_seaice_1 = phys_seaice[T,:,:]
                        diff = phys_seaice_1 - phys_seaice_0
                        difflevs = np.arange(-100,100+1E-9,10)
                        diffticklevs = np.arange(-100,100+1E-9,20)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_phys_seaice.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'n', 'none', 'none', 'none', 1.25, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('phys_cost'):
                # %%%%%%%%%%%% convective adjustments %%%%%%%%%%%%
                conv_threshold = 500
                # --- parameters ---
                levs = np.arange(0,600+1E-6,50)
                if np.max(phys_cost[T,:,:]) > conv_threshold: # if new diag in nb/yr (48 biogem timesteps in 1 year)
                    timestep_factor = 1.
                else: # still old diag in nb/timestep
                    timestep_factor = 48.
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                extend = 'max'
                cbartitle = 'Ocean convection'
                filename = exp + '_' + plottedT + '_phys_cost.png'
                # --- figure ---
                geniemap(lon_edges, lat_edges, phys_cost[T,:,:]*timestep_factor, cmap, levs, levs, 'both', lower, upper, projdata ,cbartitle, filename, 'n', 'none', 'none', 'none', 0., 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        phys_cost_0 = phys_cost[T,:,:]*timestep_factor
                    elif globalcount ==1:
                        phys_cost_1 = phys_cost[T,:,:]*timestep_factor
                        diff = phys_cost_1 - phys_cost_0
                        difflevs = np.arange(-400,400+1E-3,20)
                        diffticklevs = np.arange(-400,400+1E-3,100)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_phys_cost.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'n', 'none', 'none', 'none', 1.25, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('phys_psi'):
                # %%%%%%%%%%%% phys_psi %%%%%%%%%%%%
                # not using standard functions because possibility to highlight anormal values
                # --- parameters ---
                lookfor_issues = 'n' # 'y' if you want to check runs, 'n' if you want a publication-ready figure
                thres = 50 # threshold in Sv for anormal values
                levs = np.arange(-60,60+1E-6,2.5)
                clevs = np.arange(-1000,100+1E-6,5)
                clevs2 = np.arange(-100,100+1E-6,10)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = light_centered
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                extend = 'both'
                cbartitle = 'Barotropic streamfunction (Sv)'
                # highlight anormaly high values
                if lookfor_issues == 'y':
                    ix = np.array([])
                    sqval = np.squeeze(phys_psi[T,:,:])
                    ix = np.argwhere(abs(sqval) > thres)
                    nbanomalies = np.shape(ix)[0]
                    if nbanomalies > 0:
                        print('        Value > ' + str(thres) + ' Sv detected in ' + str(nbanomalies) + ' location(s) -- DANGER')
                        anomalflag = 'y'
                    else:
                        anomalflag = 'f'
                # --- figure ---
                fig =  plt.figure(figsize=(figXsize, figYsize))
                ax = fig.add_subplot(111,projection=projdata)
                ax = fig.gca(); ax.set(aspect=2); ax.set_aspect('auto')
                cf = plt.pcolormesh(lon_psi_edges, lat_psi_edges, phys_psi[T,:,:], transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
                ct = plt.contour(lon_psi, lat_psi, np.squeeze(phys_psi[T,:,:]), clevs, transform=ccrs.PlateCarree(), colors='k', linewidths = 0.45)
                ct2 = plt.contour(lon_psi, lat_psi, np.squeeze(phys_psi[T,:,:]), clevs2, transform=ccrs.PlateCarree(), colors='k', linewidths = 0.75)
                plt.pcolormesh(lon_edges, lat_edges, land, cmap=whitecmap, transform=ccrs.PlateCarree())
                if landflag == 1:
                    stepped_coastline_cGENIE(lon_edges,lat_edges, data_crs, grid_area_ocn,1.25)
                if ((lookfor_issues == 'y') and (anomalflag == 'y')):
                    for anom in np.arange(nbanomalies):
                        anom_lon = lon_psi[ix[anom][1]]
                        anom_lat = lat_psi[ix[anom][0]]
                        r = 35
                        plt.plot(anom_lon,anom_lat, marker='o', markersize=r, markeredgecolor = 'orange', markeredgewidth = 1.5, markerfacecolor = 'none', transform=ccrs.PlateCarree(), zorder=1000)
                        print('        lon = ' + str(anom_lon) + ', lat = ' + str(anom_lat) + ', psi = ' + str(sqval[ix[anom][0], ix[anom][1]]) + ' Sv')
                gl = ax.gridlines(crs=data_crs, draw_labels=False,
                              linewidth=0.5, color='k', alpha=0.2, linestyle='-', zorder=999)
                gl.xlabels_top = False; gl.ylabels_left = False; gl.xlines = True
                gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180]); gl.xformatter = LONGITUDE_FORMATTER
                gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90]); gl.yformatter = LATITUDE_FORMATTER
                gl.xlabel_style = {'size': plot_font_size, 'color': 'gray'}; gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
                cb = fig.colorbar(cf, orientation='horizontal',extend=extend)
                cf.cmap.set_under(lower); cf.cmap.set_over(upper)
                cb.ax.tick_params(labelsize=plot_font_size)
                cb.ax.set_title(cbartitle, weight='normal', fontsize=plot_font_size)
                plt.tight_layout()
                filename = exp + '_' + plottedT + '_phys_psi.png'
                if save_fig == 'y':
                    plt.savefig(filename,format='png')
                    filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                    os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        phys_psi_0 = phys_psi[T,:,:]
                    elif globalcount ==1:
                        phys_psi_1 = phys_psi[T,:,:]
                        diff = abs(phys_psi_1) - abs(phys_psi_0)
                        difflevs = np.arange(-15,15+1E-3,1)
                        diffclevs = np.arange(-15,15+1E-3,5)
                        diffclevs2 = np.arange(-15,15+1E-3,1)
                        diffticklevs = np.arange(-15,15+1E-3,5)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_phys_psi.png'
                        diffnorm = BoundaryNorm(difflevs, ncolors=diffcmap.N, clip=False)
                        fig =  plt.figure(figsize=(figXsize, figYsize))
                        ax = fig.add_subplot(111,projection=projdata)
                        ax = fig.gca(); ax.set(aspect=2); ax.set_aspect('auto')
                        cf = plt.pcolormesh(lon_psi_edges, lat_psi_edges, diff, transform=ccrs.PlateCarree(),cmap=diffcmap, norm=diffnorm)
                        ct = plt.contour(lon_psi, lat_psi, np.squeeze(diff), diffclevs, transform=ccrs.PlateCarree(), colors='k', linewidths = 0.45)
                        ct2 = plt.contour(lon_psi, lat_psi, np.squeeze(diff), diffclevs2, transform=ccrs.PlateCarree(), colors='k', linewidths = 0.75)
                        plt.pcolormesh(lon_edges, lat_edges, land, cmap=whitecmap, transform=ccrs.PlateCarree())
                        if landflag == 1:
                            stepped_coastline_cGENIE(lon_edges,lat_edges, data_crs, grid_area_ocn,1.25)
                        gl = ax.gridlines(crs=data_crs, draw_labels=False,
                                      linewidth=0.5, color='k', alpha=0.2, linestyle='-', zorder=999)
                        gl.xlabels_top = False; gl.ylabels_left = False; gl.xlines = True
                        gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180]); gl.xformatter = LONGITUDE_FORMATTER
                        gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90]); gl.yformatter = LATITUDE_FORMATTER
                        gl.xlabel_style = {'size': plot_font_size, 'color': 'gray'}; gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
                        cb = fig.colorbar(cf, orientation='horizontal',extend=extend)
                        cf.cmap.set_under(difflower); cf.cmap.set_over(diffupper)
                        cb.ax.tick_params(labelsize=plot_font_size)
                        cb.ax.set_title(cbartitle, weight='normal', fontsize=plot_font_size)
                        plt.tight_layout()
                        if save_fig == 'y':
                            plt.savefig(difffilename,format='png', dpi=450)
                            difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                            os.system('ln -s ' + difffilename + ' ' + difflnfile)
                            diffsavedfiles.append(difffilename)

            if dovar('bio_export_POC'):
                # %%%%%%%%%%%% bio export POC %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0,7.5+1E-9,0.5)
                ticklevs = np.arange(0,7+1E-9,1)
                clevs = np.arange(0,7.5+1E-9,2)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                extend = 'max'
                cbartitle = 'Biological export -- POC (mol m$^{-2}$ yr$^{-1}$)'
                filename = exp + '_' + plottedT + '_bio_export_POC.png'
                # --- figure ---
                geniemap(lon_edges, lat_edges, bio_export_POC[T,:,:], cmap, levs, ticklevs, 'both', lower, upper, projdata ,cbartitle, filename, 'y', lon, lat, clevs, 0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        bio_export_POC_0 = bio_export_POC[T,:,:]
                    elif globalcount ==1:
                        bio_export_POC_1 = bio_export_POC[T,:,:]
                        diff = bio_export_POC_1 - bio_export_POC_0
                        difflevs = np.arange(-2.,2.+1E-3,0.2)
                        diffclevs = np.arange(-2.,2.+1E-3,0.5)
                        diffticklevs = np.arange(-2.,2.+1E-3,0.5)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_bio_export_POC.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'y', lon, lat, difflevs, 0.75, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('bio_diag_k_PO4'):
                # %%%%%%%%%%%% biological productivity control - k_PO4 %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0,1+1E-9,0.1)
                ticklevs = np.arange(0,1+1E-9,0.5)
                clevs = np.arange(0,1+1E-9,0.5)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                extend = 'max'
                cbartitle = 'Biological productivity control - k_PO4'
                filename = exp + '_' + plottedT + '_bio_diag_k_PO4.png'
                # --- figure ---
                geniemap(lon_edges, lat_edges, bio_diag_k_PO4[T,:,:], cmap, levs, ticklevs, 'both', lower, upper, projdata ,cbartitle, filename, 'y', lon, lat, clevs, 0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        bio_diag_k_PO4_0 = bio_diag_k_PO4[T,:,:]
                    elif globalcount ==1:
                        bio_diag_k_PO4_1 = bio_diag_k_PO4[T,:,:]
                        diff = bio_diag_k_PO4_1 - bio_diag_k_PO4_0
                        difflevs = np.arange(-1,1+1E-3,0.1)
                        diffclevs = np.arange(-1,1+1E-3,0.25)
                        diffticklevs = np.arange(-1,1+1E-3,0.5)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_bio_diag_kPO4.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'y', lon, lat, diffclevs, 0.75, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('bio_diag_k_temp'):
                # %%%%%%%%%%%% biological productivity control - k_temp %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0,5+1E-9,0.25)
                ticklevs = np.arange(0,5+1E-9,1)
                clevs = np.arange(0,5+1E-9,0.5)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                extend = 'max'
                cbartitle = 'Biological productivity control - k_temp'
                filename = exp + '_' + plottedT + '_bio_diag_k_temp.png'
                # --- figure ---
                geniemap(lon_edges, lat_edges, bio_diag_k_temp[T,:,:], cmap, levs, ticklevs, 'both', lower, upper, projdata ,cbartitle, filename, 'y', lon, lat, clevs, 0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        bio_diag_k_temp_0 = bio_diag_k_temp[T,:,:]
                    elif globalcount ==1:
                        bio_diag_k_temp_1 = bio_diag_k_temp[T,:,:]
                        diff = bio_diag_k_temp_1 - bio_diag_k_temp_0
                        difflevs = np.arange(-1,1+1E-3,0.1)
                        diffclevs = np.arange(-1,1+1E-3,0.25)
                        diffticklevs = np.arange(-1,1+1E-3,0.5)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_bio_diag_k_temp.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'y', lon, lat, diffclevs, 0.75, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('bio_diag_k_light'):
                # %%%%%%%%%%%% biological productivity control - k_light %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0,1+1E-9,0.1)
                difflevs = np.arange(-1,1+1E-3,0.1)
                diffclevs = np.arange(-1,1+1E-3,0.25)
                diffticklevs = np.arange(-1,1+1E-3,0.5)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                extend = 'max'
                cbartitle = 'Biological productivity control - k_light'
                filename = exp + '_' + plottedT + '_bio_diag_k_light.png'
                # --- figure ---
                geniemap(lon_edges, lat_edges, bio_diag_k_light[T,:,:], cmap, levs, ticklevs, 'both', lower, upper, projdata ,cbartitle, filename, 'y', lon, lat, clevs, 0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        bio_diag_k_light_0 = bio_diag_k_light[T,:,:]
                    elif globalcount ==1:
                        bio_diag_k_light_1 = bio_diag_k_light[T,:,:]
                        diff = bio_diag_k_light_1 - bio_diag_k_light_0
                        difflevs = np.arange(-1,1+1E-3,0.05)
                        diffclevs = np.arange(-1,1+1E-3,0.2)
                        diffticklevs = np.arange(-1,1+1E-3,0.5)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_bio_diag_k_light.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'y', lon, lat, diffclevs, 0.75, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('ocn_PO4'):
                # %%%%%%%%%%%% PO4 %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0,0.3+1E-9,0.025)
                ticklevs = np.arange(0,0.3+1E-9,0.05)
                clevs = np.arange(0,0.3+1E-9,0.05)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                extend = 'max'
                cbartitle = 'PO$_4$ ($\mu$mol kg$^{-1}$)'
                filename = exp + '_' + plottedT + '_ocn_PO4.png'
                # --- figure ---
                geniemap(lon_edges, lat_edges, ocn_PO4[T,0,:,:]*1E6, cmap, levs, ticklevs, 'max', lower, upper, projdata ,cbartitle, filename, 'y', lon, lat, clevs, 0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        ocn_PO4_0 = ocn_PO4[T,0,:,:]*1E6
                    elif globalcount ==1:
                        ocn_PO4_1 = ocn_PO4[T,0,:,:]*1E6
                        diff = ocn_PO4_1 - ocn_PO4_0
                        difflevs = np.arange(-0.16,0.16+1E-3,0.02)
                        diffclevs = np.arange(-0.16,0.16+1E-3,0.04)
                        diffticklevs = np.arange(-0.16,0.16+1E-3,0.08)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_ocn_PO4.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'y', lon, lat, difflevs, 0.75, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('deepO2'):
                # %%%%%%%%%%%% deep O2 %%%%%%%%%%%%
                # --- parameters ---
                levs = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250])
                clevs = np.array([0.])
                ticklevs = levs
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cbartitle='O2 $\mu$mol L$^{-1}$'
                cmap = fzcmap_alpha065
                filename = exp + '_' + plottedT + '_O2_deep.png' 
                # --- figure ---
                geniemap(lon_edges, lat_edges, deepO2, cmap, levs, ticklevs, 'both', lower, upper, projdata ,cbartitle, filename, 'y', lon, lat, clevs, 1.2, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                if do_diff == 'y':
                    if globalcount ==0:
                        deepO2_0 = deepO2
                    elif globalcount ==1:
                        deepO2_1 = deepO2
                        diff = deepO2_1 - deepO2_0
                        difflevs = np.arange(-100,100+1E-3,10)
                        diffclevs = np.arange(-100,100+1E-3,20)
                        diffticklevs = np.arange(-100,100+1E-3,50)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_deepO2.png'
                        geniemap(lon_edges, lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'y', lon, lat, diffclevs, 1., 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('deepO2_shallowpfonly'):
                # %%%%%%%%%%%% deep O2 shallow-platforms only %%%%%%%%%%%%
                # not using standard function because outlines the land mask, not the variable mask
                # --- parameters ---
                levs = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250])
                clevs = np.array([0.])
                ticklevs = levs
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cbartitle='O2 $\mu$mol L$^{-1}$'
                cmap = fzcmap_alpha065
                filename = exp + '_' + plottedT + '_O2_deep_shallowpfonly' + str(shallow_threshold)+ 'm.png'
                # --- figure ---
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                fig =  plt.figure(figsize=(figXsize, figYsize))
                ax = fig.add_subplot(111,projection=projdata)
                ax = fig.gca(); ax.set(aspect=2); ax.set_aspect('auto')
                cf = plt.pcolormesh(lon_edges, lat_edges, deepO2_shallowpfonly, transform=ccrs.PlateCarree(),cmap=cmap, norm=norm)
                ct = plt.contour(lon, lat, np.squeeze(deepO2_shallowpfonly), clevs, transform=ccrs.PlateCarree(), colors='k', linewidths = 1.2)
                plt.pcolormesh(lon_edges, lat_edges, land, cmap=whitecmap, transform=ccrs.PlateCarree())
                if landflag == 1:
                    stepped_coastline_cGENIE(lon_edges,lat_edges, data_crs, grid_area_ocn,1.25)
                    stepped_coastline_cGENIE(lon_edges,lat_edges, data_crs, pfmask,0.7)
                gl = ax.gridlines(crs=data_crs, draw_labels=False,
                              linewidth=0.5, color='k', alpha=0.2, linestyle='-', zorder=999)
                gl.xlabels_top = False; gl.ylabels_left = False; gl.xlines = True
                gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180]); gl.xformatter = LONGITUDE_FORMATTER
                gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90]); gl.yformatter = LATITUDE_FORMATTER
                gl.xlabel_style = {'size': plot_font_size, 'color': 'gray'}; gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
                cb = fig.colorbar(cf, orientation='horizontal',extend='both',ticks=ticklevs)
                cf.cmap.set_under(lower); cf.cmap.set_over(upper)
                cb.ax.tick_params(labelsize=plot_font_size)
                cb.ax.set_title(cbartitle, weight='normal', fontsize=plot_font_size)
                plt.tight_layout()
                if save_fig == 'y':
                    plt.savefig(filename,format='png')
                    filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                    os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        deepO2_shallowpfonly_0 = deepO2_shallowpfonly
                    elif globalcount ==1:
                        deepO2_shallowpfonly_1 = deepO2_shallowpfonly
                        diff = deepO2_shallowpfonly_1 - deepO2_shallowpfonly_0
                        difflevs = np.arange(-100,100+1E-3,10)
                        diffclevs = np.arange(-100,100+1E-3,20)
                        diffticklevs = np.arange(-100,100+1E-3,50)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_deepO2_shallowpfonly' + str(shallow_threshold)+ 'm.png'
                        diffnorm = BoundaryNorm(difflevs, ncolors=diffcmap.N, clip=False)
                        fig =  plt.figure(figsize=(figXsize, figYsize))
                        ax = fig.add_subplot(111,projection=projdata)
                        ax = fig.gca(); ax.set(aspect=2); ax.set_aspect('auto')
                        cf = plt.pcolormesh(lon_edges, lat_edges, diff, transform=ccrs.PlateCarree(),cmap=diffcmap, norm=diffnorm)
                        ct = plt.contour(lon, lat, np.squeeze(diff), diffclevs, transform=ccrs.PlateCarree(), colors='k', linewidths = 1.)
                        plt.pcolormesh(lon_edges, lat_edges, land, cmap=whitecmap, transform=ccrs.PlateCarree())
                        if landflag == 1:
                            stepped_coastline_cGENIE(lon_edges,lat_edges, data_crs, grid_area_ocn,1.25)
                            stepped_coastline_cGENIE(lon_edges,lat_edges, data_crs, pfmask,0.7)
                        gl = ax.gridlines(crs=data_crs, draw_labels=False,
                                      linewidth=0.5, color='k', alpha=0.2, linestyle='-', zorder=999)
                        gl.xlabels_top = False; gl.ylabels_left = False; gl.xlines = True
                        gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180]); gl.xformatter = LONGITUDE_FORMATTER
                        gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90]); gl.yformatter = LATITUDE_FORMATTER
                        gl.xlabel_style = {'size': plot_font_size, 'color': 'gray'}; gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
                        cb = fig.colorbar(cf, orientation='horizontal',extend=diffextend,ticks=diffticklevs)
                        cf.cmap.set_under(difflower); cf.cmap.set_over(diffupper)
                        cb.ax.tick_params(labelsize=plot_font_size)
                        cb.ax.set_title(cbartitle, weight='normal', fontsize=plot_font_size)
                        plt.tight_layout()
                        if save_fig == 'y':
                            plt.savefig(difffilename,format='png', dpi=450)
                            difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                            os.system('ln -s ' + difffilename + ' ' + difflnfile)
                            diffsavedfiles.append(difffilename)

            if dovar('lat_ocn_O2'):
                # %%%%%%%%%%%% zonal O2 %%%%%%%%%%%%        
                # --- parameters ---
                levs = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250])
                clevs = np.array([0.])
                ticklevs = levs
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap_alpha065
                cbartitle='O2 $\mu$mol L$^{-1}$'
                filename = exp + '_' + plottedT + '_zonal_O2.png'
                # --- figure ---
                genielat(lat_edges,-1*zt_edges/1000.,lat_ocn_O2, cmap, levs, ticklevs, 'both', lower, upper, cbartitle, filename, 'y', lat, -1*zt/1000., clevs, 1.)        
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        lat_ocn_O2_0 = lat_ocn_O2
                    elif globalcount ==1:
                        lat_ocn_O2_1 = lat_ocn_O2
                        diff = lat_ocn_O2_1 - lat_ocn_O2_0
                        difflevs = np.arange(-100,100+1E-3,10)
                        diffticklevs = np.arange(-100,100+1E-3,50)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_lat_ocn_O2.png'
                        genielat(lat_edges,-1*zt_edges/1000.,diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, cbartitle, difffilename, 'y', lat, -1*zt/1000., difflevs, 0.75)
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('O2'):
                # %%%%%%%%%%%% O2 per level %%%%%%%%%%%%
                # --- parameters ---
                levs = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250]) # default
                #levs = np.arange(0,100+1E-9,10) # For 0.4 O2 0.4 PO4
                clevs = np.array([0.])
                ticklevs =  np.array([0, 25, 50, 100, 150, 200, 250]) # default
                #ticklevs =  np.arange(0,100+1E-9,20) # For 0.4 O2 0.4 PO4
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap3 #fzcmap_alpha065
                filename = exp + '_' + plottedT + '_O2_lev.png'
                # --- figure ---
                genielev(lon_edges, lat_edges, np.squeeze(O2), zt, cmap, levs, ticklevs, 'both', lower, upper, projdata ,'O$_2$ concentration ($\mu$mol kg$^{-1}$)', filename, 'y', lon, lat, clevs)
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        O2_0 = np.squeeze(O2)
                    elif globalcount ==1:
                        O2_1 = np.squeeze(O2)
                        diff = O2_1 - O2_0
                        difflevs = np.arange(-80,80+1E-3,10)
                        diffclevs = np.arange(-80,80+1E-3,20)
                        diffticklevs = np.arange(-80,80+1E-3,40)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_O2_lev.png'
                        genielev(lon_edges, lat_edges, diff, zt, diffcmap, difflevs, diffticklevs, 'both', difflower, diffupper, projdata ,'O$_2$ concentration ($\mu$mol kg$^{-1}$)', difffilename, 'y', lon, lat, diffclevs)
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('ma_lat_ocn_H2S'):
                # %%%%%%%%%%%% zonal H2S %%%%%%%%%%%%i
                # --- parameters ---
                levs = np.arange(0, 100+1E-6, 10)
                clevs = np.array([1.])
                ticklevs = levs
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap_alpha065
                cbartitle = 'H2S $\mu$mol L$^{-1}$'
                filename = exp + '_' + plottedT + '_zonal_H2S.png'
                # --- figure ---
                genielat(lat_edges,-1*zt_edges/1000.,ma_lat_ocn_H2S, cmap, levs, ticklevs, 'both', lower, upper, cbartitle, filename, 'y', lat, -1*zt/1000., clevs, 1.)
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        lat_ocn_H2S_0 = lat_ocn_H2S
                    elif globalcount ==1:
                        lat_ocn_H2S_1 = lat_ocn_H2S
                        diff = lat_ocn_H2S_1 - lat_ocn_H2S_0
                        difflevs = np.arange(-50,50+1E-3,5)
                        diffclevs = np.arange(-50,50+1E-3,10)
                        diffticklevs = np.arange(-50,50+1E-3,25)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_lat_ocn_H2S.png'
                        genielat(lat_edges,-1*zt_edges/1000.,diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, cbartitle, difffilename, 'y', lat, -1*zt/1000., diffclevs, 0.5)
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('ma_H2S'):
                # %%%%%%%%%%%% H2S per level %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0, 100+1E-6, 10)
                clevs = np.array([1.])
                ticklevs = levs
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap_alpha065
                filename = exp + '_' + plottedT + '_H2S_lev.png'
                # --- figure ---
                genielev(lon_edges, lat_edges, ma_H2S, zt, cmap, levs, ticklevs, 'both', lower, upper, projdata ,'H2S $\mu$mol L$^{-1}$', filename, 'y', lon, lat, clevs)
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        maH2S_0 = ma_H2S
                    elif globalcount ==1:
                        maH2S_1 = ma_H2S
                        diff = maH2S_1 - maH2S_0
                        difflevs = np.arange(-50,50+1E-3,5)
                        diffclevs = np.arange(-50,50+1E-3,10)
                        diffticklevs = np.arange(-50,50+1E-3,25)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_H2S_lev.png'
                        genielev(lon_edges, lat_edges, diff, zt, diffcmap, difflevs, diffticklevs, 'both', difflower, diffupper, projdata ,'H2S $\mu$mol L$^{-1}$', difffilename, 'y', lon, lat, diffclevs)
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('misc_col_Dage'):
                # %%%%%%%%%%%% ventilation age per level %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0, 2000+1E-6, 200)
                clevs = np.arange(0, 2000+1E-6, 400)
                ticklevs = np.arange(0, 2000+1E-6, 400)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap_alpha065
                filename = exp + '_' + plottedT + '_ventilation_age.png'
                # --- figure ---
                genielev(lon_edges, lat_edges,np.squeeze(misc_col_Dage[T,:,:,:]), zt, cmap, levs, ticklevs, 'max', lower, upper, projdata ,'Ventilation age (years)', filename, 'n', lon, lat, clevs)
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        misc_col_Dage_0 = np.squeeze(misc_col_Dage[T,:,:,:])
                    elif globalcount ==1:
                        misc_col_Dage_1 = np.squeeze(misc_col_Dage[T,:,:,:])
                        diff = misc_col_Dage_1 - misc_col_Dage_0
                        difflevs = np.arange(-500,500+1E-3,50)
                        diffclevs = np.arange(-500,500+1E-3,100)
                        diffticklevs = np.arange(-500,500+1E-3,250)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_ventilation_age.png'
                        genielev(lon_edges, lat_edges, diff, zt, diffcmap, difflevs, diffticklevs, 'both', difflower, diffupper, projdata ,'Ventilation age (years)', difffilename, 'y', lon, lat, diffclevs)
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('lat_misc_col_Dage'):
                # %%%%%%%%%%%% zonal ventilation age %%%%%%%%%%%%i
                # --- parameters ---
                levs = np.arange(0, 2000+1E-6, 100)
                clevs = np.arange(0, 2000+1E-6, 200)
                ticklevs = np.arange(0, 2000+1E-6, 400)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap_alpha065
                cbartitle = 'Ventilation age (yrs)'
                filename = exp + '_' + plottedT + '_zonal_ventilation_age.png'
                # --- figure ---
                genielat(lat_edges,-1*zt_edges/1000.,lat_misc_col_Dage, cmap, levs, ticklevs, 'both', lower, upper, cbartitle, filename, 'y', lat, -1*zt/1000., clevs, 0.75)
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        lat_misc_col_Dage_0 = lat_misc_col_Dage
                    elif globalcount ==1:
                        lat_misc_col_Dage_1 = lat_misc_col_Dage
                        diff = lat_misc_col_Dage_1 - lat_misc_col_Dage_0
                        difflevs = np.arange(-500,500+1E-3,50)
                        diffclevs = np.arange(-500,500+1E-3,100)
                        diffticklevs = np.arange(-500,500+1E-3,250)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_zonal_ventilation_age.png'
                        genielat(lat_edges,-1*zt_edges/1000.,diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, cbartitle, difffilename, 'y', lat, -1*zt/1000., diffclevs, 0.5)
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('ocn_DIC_13C'):
                # %%%%%%%%%%%% d13C of DIC surface %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0.5, 3+1E-6, 0.1)
                clevs = np.arange(0.5, 3+1E-6, 0.2)
                ticklevs = np.arange(0.5, 3+1E-6, 0.5)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap_alpha065
                extend = 'both'
                cbartitle = 'd13C of DIC (permil)'
                filename = exp + '_' + plottedT + '_d13C_DIC_surface.png'
                # --- figure ---
                geniemap(lon_edges, lat_edges, ocn_DIC_13C[T,0,:,:], cmap, levs, ticklevs, extend, lower, upper, projdata ,cbartitle, filename, 'y', lon, lat, clevs,0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        ocn_DIC_13C_surf_0 = ocn_DIC_13C[T,0,:,:]
                    elif globalcount ==1:
                        ocn_DIC_13C_surf_1 = ocn_DIC_13C[T,0,:,:]
                        diff = ocn_DIC_13C_surf_1 - ocn_DIC_13C_surf_0
                        difflevs = np.arange(-1.5, 1.5+1E-6, 0.1)
                        diffclevs = np.arange(-1.5, 1.5+1E-6, 0.1)
                        diffticklevs = np.arange(-1.5, 1.5+1E-6, 1)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_ocn_DIC_13C_surf.png'
                        geniemap(lon_edges, lat_edges,  diff, diffcmap, difflevs, diffticklevs, diffextend, lower, upper, projdata ,cbartitle, difffilename, 'y', lon, lat, diffclevs, 0.75, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('lat_ocn_DIC_13C'):
                # %%%%%%%%%%%% zonal d13C DIC %%%%%%%%%%%%i
                # --- parameters ---
                levs = np.arange(-2, 2+1E-6, 0.25)
                clevs = np.arange(-2, 2+1E-6, 0.5)
                ticklevs = np.arange(-2, 2+1E-6, 0.5)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap_alpha065
                extend = 'both'
                cbartitle = 'd13C of DIC (permil)'
                filename = exp + '_' + plottedT + '_zonal_d13C_DIC.png'
                # --- figure ---
                genielat(lat_edges,-1*zt_edges/1000.,lat_ocn_DIC_13C, cmap, levs, ticklevs, 'both', lower, upper, cbartitle, filename, 'y', lat, -1*zt/1000., clevs, 0.75)
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        lat_ocn_DIC_13C_0 = lat_ocn_DIC_13C
                    elif globalcount ==1:
                        lat_ocn_DIC_13C_1 = lat_ocn_DIC_13C
                        diff = lat_ocn_DIC_13C_1 - lat_ocn_DIC_13C_0
                        difflevs = np.arange(-2, 2+1E-6, 0.1)
                        diffclevs = np.arange(-2, 2+1E-6, 0.2)
                        diffticklevs = np.arange(-2, 2+1E-6, 1)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_zonal_ocn_DIC_13C.png'
                        genielat(lat_edges,-1*zt_edges/1000.,diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, cbartitle, difffilename, 'y', lat, -1*zt/1000., diffclevs, 0.5)
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('ocn_DIC_13C'):
                # %%%%%%%%%%%% d13C of DIC per level %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(-2, 2+1E-6, 0.1)
                clevs = np.arange(-2, 2+1E-6, 0.2)
                ticklevs = np.arange(-2, 2+1E-6, 1.)
                extend = 'both'
                filename = exp + '_' + plottedT + '_ocn_DIC_13C_lev.png'
                # --- figure ---
                genielev(lon_edges, lat_edges, np.squeeze(ocn_DIC_13C[T,:,:,:]), zt, diffcmap, levs, ticklevs, extend, difflower, diffupper, projdata ,'d13C of DIC (permil)', filename, 'y', lon, lat, clevs)
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        ocn_DIC_13C_0 = ocn_DIC_13C[T,:,:,:]
                    elif globalcount ==1:
                        ocn_DIC_13C_1 = ocn_DIC_13C[T,:,:,:]
                        diff = ocn_DIC_13C_1 - ocn_DIC_13C_0
                        difflevs = np.arange(-2, 2+1E-6, 0.1)
                        diffclevs = np.arange(-2, 2+1E-6, 0.2)
                        diffticklevs = np.arange(-2, 2+1E-6, 1)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_ocn_DIC_13C_lev.png'
                        genielev(lon_edges, lat_edges, np.squeeze(diff), zt, diffcmap, difflevs, diffticklevs, extend, difflower, diffupper, projdata ,'d13C of DIC (permil)', difffilename, 'y', lon, lat, diffclevs)
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

        ########################################################################
        #                                  SEDGEM                              #
        ########################################################################

        if do_sedgem == 'y':

            # ===================== CALCULATIONS =====================

            Tsed = 0 # sedgem only saves the final output

            # building fake land-sea mask - required by stepped_coastline
            sed_landsea_mask = np.where(sed_grid_topo.mask==True,-1,1)
            sed_land = np.ma.masked_where(sed_landsea_mask==1,sed_landsea_mask)

            # ===================== PLOTTING =====================

            if dovar('sed_grid_topo'):
                # %%%%%%%%%%%% sed_grid_topo %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0., 5.5+1E-6,0.5)
                ticklevs = levs
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                cbartitle = 'Bathymetry (km. b.s.l.)'
                filename = exp + '_' + plottedT + '_sed_grid_topo.png'
                # --- figure ---
                geniemap(sed_lon_edges, sed_lat_edges, (-1)*sed_grid_topo*1E-3, cmap, levs, ticklevs, 'max', lower, upper, projdata ,cbartitle, filename, 'n', 'none', 'none', 'none', 0, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        sed_grid_topo_0 = (-1)*sed_grid_topo*1E-3
                    elif globalcount ==1:
                        sed_grid_topo_1 = (-1)*sed_grid_topo*1E-3
                        diff = sed_grid_topo_1 - sed_grid_topo_0
                        difflevs = np.arange(-5.,5.+1E-3,0.5)
                        diffticklevs = np.arange(-5.,5.+1E-3,1)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_sed_grid_topo.png'
                        geniemap(sed_lon_edges, sed_lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'n', 'none', 'none', 'none', 0, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('sed_ocn_O2'):
                # %%%%%%%%%%%% sed_ocn_O2 %%%%%%%%%%%%
                # --- parameters ---
                levs = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250])
                clevs = np.array([0.])
                ticklevs = levs
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cbartitle='O2 $\mu$mol L$^{-1}$'
                cmap = fzcmap_alpha065
                filename = exp + '_' + plottedT + '_sed_ocn_O2.png'
                # --- figure ---
                geniemap(sed_lon_edges, sed_lat_edges, sed_ocn_O2[Tsed,:,:]*1E6, cmap, levs, ticklevs, 'both', lower, upper, projdata ,cbartitle, filename, 'y', sed_lon, sed_lat, clevs, 1.2, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                if do_diff == 'y':
                    if globalcount ==0:
                        sed_ocn_O2_0 = sed_ocn_O2[Tsed,:,:]
                    elif globalcount ==1:
                        sed_ocn_O2_1 = sed_ocn_O2[Tsed,:,:]
                        diff = sed_ocn_O2_1 - sed_ocn_O2_0
                        difflevs = np.arange(-100,100+1E-3,10)
                        diffclevs = np.arange(-100,100+1E-3,20)
                        diffticklevs = np.arange(-100,100+1E-3,50)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_sed_ocn_O2.png'
                        geniemap(sed_lon_edges, sed_lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'y', sed_lon, sed_lat, diffclevs, 1., 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('sed_ocn_DIC_13C'):
                # %%%%%%%%%%%% ocn_DIC_13C %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(-10., 1.+1E-6, 0.5)
                clevs = np.arange(-10., 1.+1E-6, 2.)
                ticklevs = np.arange(-10., 1.+1E-6, 1)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap_alpha065
                extend = 'both'
                cbartitle = 'sed_ocn_DIC_13C (permil)'
                filename = exp + '_' + plottedT + '_sed_ocn_DIC_13C.png'
                # --- figure ---
                geniemap(sed_lon_edges, sed_lat_edges, sed_ocn_DIC_13C[Tsed,:,:], cmap, levs, ticklevs, extend, lower, upper, projdata ,cbartitle, filename, 'y', sed_lon, sed_lat, clevs,0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        sed_ocn_DIC_13C_0 = sed_ocn_DIC_13C[Tsed,:,:]
                    elif globalcount ==1:
                        sed_ocn_DIC_13C_1 = sed_ocn_DIC_13C[Tsed,:,:]
                        diff = sed_ocn_DIC_13C_1 - sed_ocn_DIC_13C_0
                        difflevs = np.arange(-1.5, 1.5+1E-6, 0.1)
                        diffclevs = np.arange(-1.5, 1.5+1E-6, 0.1)
                        diffticklevs = np.arange(-1.5, 1.5+1E-6, 1)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_sed_ocn_DIC_13C_surf.png'
                        geniemap(sed_lon_edges, sed_lat_edges,  diff, diffcmap, difflevs, diffticklevs, diffextend, lower, upper, projdata ,cbartitle, difffilename, 'y', sed_lon, sed_lat, diffclevs, 0.75, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('OMEN_wtpct_top'):
                # %%%%%%%%%%%% OMEN_wtpct_top %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0,25+1E-6,1)
                ticklevs = np.arange(0,25+1E-6,5)
                clevs = np.arange(0,25+1E-6,5)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                cbartitle = 'OMEN_wtpct_top'
                filename = exp + '_' + plottedT + '_OMEN_wtpct_top.png'
                # --- figure ---
                geniemap(sed_lon_edges, sed_lat_edges, OMEN_wtpct_top[Tsed,:,:], cmap, levs, ticklevs, 'max', lower, upper, projdata ,cbartitle, filename, 'y', sed_lon, sed_lat, clevs, 0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        OMEN_wtpct_top_0 = OMEN_wtpct_top[Tsed,:,:]
                    elif globalcount ==1:
                        OMEN_wtpct_top_1 = OMEN_wtpct_top[Tsed,:,:]
                        diff = OMEN_wtpct_top_1 - OMEN_wtpct_top_0
                        difflevs = np.array([-100, -90, -80, -70, -60, -50, -40, -30, -20, -10, -5, -2.5, -1, 0, 1, 2.5, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
                        diffticklevs = np.array([-100, -50, -20, -10, -5, -2.5, -1, 1, 2.5, 5, 10, 20, 50, 100])
                        diffextend = 'both'
                        diffclevs = np.array([-100, -50, -20, -10, -5, -2.5, 0, 2.5, 5, 10, 20, 50, 100])
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_OMEN_wtpct_top.png'
                        geniemap(sed_lon_edges, sed_lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'n', 'none', 'none', 'none', 0, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('OMEN_wtpct_bot'):
                # %%%%%%%%%%%% OMEN_wtpct_bot %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0,25+1E-6,1)
                ticklevs = np.arange(0,25+1E-6,5)
                clevs = np.arange(0,25+1E-6,5)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.90)
                cmap = fzcmap_alpha065
                norm = BoundaryNorm(levs, ncolors=cmap.N, clip=False)
                cbartitle = 'OMEN_wtpct_bot'
                filename = exp + '_' + plottedT + '_OMEN_wtpct_bot.png'
                # --- figure ---
                geniemap(sed_lon_edges, sed_lat_edges, OMEN_wtpct_bot[Tsed,:,:], cmap, levs, ticklevs, 'max', lower, upper, projdata ,cbartitle, filename, 'y', sed_lon, sed_lat, clevs, 0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        OMEN_wtpct_bot_0 = OMEN_wtpct_bot[Tsed,:,:]
                    elif globalcount ==1:
                        OMEN_wtpct_bot_1 = OMEN_wtpct_bot[Tsed,:,:]
                        diff = OMEN_wtpct_bot_1 - OMEN_wtpct_bot_0
                        difflevs = np.array([-100, -90, -80, -70, -60, -50, -40, -30, -20, -10, -5, -2.5, -1, 0, 1, 2.5, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
                        diffticklevs = np.array([-100, -50, -20, -10, -5, -2.5, -1, 1, 2.5, 5, 10, 20, 50, 100])
                        diffextend = 'both'
                        diffclevs = np.array([-100, -50, -20, -10, -5, -2.5, 0, 2.5, 5, 10, 20, 50, 100])
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_OMEN_wtpct_bot.png'
                        geniemap(sed_lon_edges, sed_lat_edges, diff, diffcmap, difflevs, diffticklevs, diffextend, difflower, diffupper, projdata ,cbartitle, difffilename, 'n', 'none', 'none', 'none', 0, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('sedocn_fnet_DIC_13C'):
                # %%%%%%%%%%%% sedocn_fnet_DIC_13C %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(-33.5, -31.5+1E-6, 0.25)
                clevs = np.arange(-33.5, -31.5+1E-6, 0.25)
                ticklevs = np.arange(-33.5, -31.5+1E-6, 0.5)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap_alpha065
                extend = 'both'
                cbartitle = 'DIC_13C benthic interface exchange flux (permil)'
                filename = exp + '_' + plottedT + '_sedocn_fnet_DIC_13C.png'
                # --- figure ---
                geniemap(sed_lon_edges, sed_lat_edges, sedocn_fnet_DIC_13C[Tsed,:,:], cmap, levs, ticklevs, extend, lower, upper, projdata ,cbartitle, filename, 'y', sed_lon, sed_lat, clevs,0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        sedocn_fnet_DIC_13C_0 = sedocn_fnet_DIC_13C[Tsed,:,:]
                    elif globalcount ==1:
                        sedocn_fnet_DIC_13C_1 = sedocn_fnet_DIC_13C[Tsed,:,:]
                        diff = sedocn_fnet_DIC_13C_1 - sedocn_fnet_DIC_13C_0
                        difflevs = np.arange(-1.5, 1.5+1E-6, 0.1)
                        diffclevs = np.arange(-1.5, 1.5+1E-6, 0.1)
                        diffticklevs = np.arange(-1.5, 1.5+1E-6, 1)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_sedocn_fnet_DIC_13C.png'
                        geniemap(sed_lon_edges, sed_lat_edges,  diff, diffcmap, difflevs, diffticklevs, diffextend, lower, upper, projdata ,cbartitle, difffilename, 'y', sed_lon, sed_lat, diffclevs, 0.75, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('sedocn_fnet_PO4'):
                # %%%%%%%%%%%% sedocn_fnet_PO4 %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0., 1E-5, 0.1E-5)
                clevs = np.arange(0., 1E-5, 0.1E-5)
                ticklevs = np.arange(0., 1E-5, 2.E-5)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap_alpha065
                extend = 'both'
                cbartitle = 'benthic interface exchange flux (permil)'
                filename = exp + '_' + plottedT + '_sedocn_fnet_PO4.png'
                # --- figure ---
                geniemap(sed_lon_edges, sed_lat_edges, sedocn_fnet_PO4[Tsed,:,:], cmap, levs, ticklevs, extend, lower, upper, projdata ,cbartitle, filename, 'y', sed_lon, sed_lat, clevs,0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        sedocn_fnet_PO4_0 = sedocn_fnet_PO4[Tsed,:,:]
                    elif globalcount ==1:
                        sedocn_fnet_PO4_1 = sedocn_fnet_PO4[Tsed,:,:]
                        diff = sedocn_fnet_PO4_1 - sedocn_fnet_PO4_0
                        difflevs = np.arange(-1.5, 1.5+1E-6, 0.1)
                        diffclevs = np.arange(-1.5, 1.5+1E-6, 0.1)
                        diffticklevs = np.arange(-1.5, 1.5+1E-6, 1)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_sedocn_fnet_PO4.png'
                        geniemap(sed_lon_edges, sed_lat_edges,  diff, diffcmap, difflevs, diffticklevs, diffextend, lower, upper, projdata ,cbartitle, difffilename, 'y', sed_lon, sed_lat, diffclevs, 0.75, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('fburial_POC_13C'):
                # %%%%%%%%%%%% fburial_POC_13C %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(-33.5, -31.5+1E-6, 0.25)
                clevs = np.arange(-33.5, -31.5+1E-6, 0.25)
                ticklevs = np.arange(-33.5, -31.5+1E-6, 0.5)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap_alpha065
                extend = 'both'
                cbartitle = 'sediment burial flux - POC_13C'
                filename = exp + '_' + plottedT + '_fburial_POC_13C.png'
                # --- figure ---
                geniemap(sed_lon_edges, sed_lat_edges, fburial_POC_13C[Tsed,:,:], cmap, levs, ticklevs, extend, lower, upper, projdata ,cbartitle, filename, 'y', sed_lon, sed_lat, clevs,0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        fburial_POC_13C_0 = fburial_POC_13C[Tsed,:,:]
                    elif globalcount ==1:
                        fburial_POC_13C_1 = fburial_POC_13C[Tsed,:,:]
                        diff = fburial_POC_13C_1 - fburial_POC_13C_0
                        difflevs = np.arange(-1.5, 1.5+1E-6, 0.1)
                        diffclevs = np.arange(-1.5, 1.5+1E-6, 0.1)
                        diffticklevs = np.arange(-1.5, 1.5+1E-6, 1)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_fburial_POC_13C.png'
                        geniemap(sed_lon_edges, sed_lat_edges,  diff, diffcmap, difflevs, diffticklevs, diffextend, lower, upper, projdata ,cbartitle, difffilename, 'y', sed_lon, sed_lat, diffclevs, 0.75, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                        diffsavedfiles.append(difffilename)

            if dovar('fburial_det'):
                # %%%%%%%%%%%% fburial_det %%%%%%%%%%%%
                # --- parameters ---
                levs = np.arange(0, 10E-6, 1E-6)
                ticklevs = np.arange(0, 10E-6, 5E-6)
                lower = fakealpha(mpl.colors.to_rgba('darkblue')[0:3],0.65)
                upper = fakealpha(mpl.colors.to_rgba('darkred')[0:3],0.75)
                cmap = fzcmap_alpha065
                extend = 'both'
                cbartitle = 'sediment burial flux - det'
                filename = exp + '_' + plottedT + '_fburial_det.png'
                # --- figure ---
                geniemap(sed_lon_edges, sed_lat_edges, fburial_det[Tsed,:,:], cmap, levs, ticklevs, extend, lower, upper, projdata ,cbartitle, filename, 'n', 'none', 'none', 'none',0.75, 'n', 'none', 'none','png')
                filecount += 1; lnfile = 'file' + str(filecount) + '.png'
                os.system('ln -s ' + filename + ' ' + lnfile)
                savedfiles.append(filename)
                # --- diff ---
                if do_diff == 'y':
                    if globalcount ==0:
                        fburial_det_0 = fburial_det[Tsed,:,:]
                    elif globalcount ==1:
                        fburial_det_1 = fburial_det[Tsed,:,:]
                        diff = fburial_det_1 - fburial_det_0
                        difflevs = np.arange(-1.5, 1.5+1E-6, 0.1)
                        diffclevs = np.arange(-1.5, 1.5+1E-6, 0.1)
                        diffticklevs = np.arange(-1.5, 1.5+1E-6, 1)
                        diffextend = 'both'
                        difffilename = exp1 + '_' + plottedT1 + '_minus_' + exp0 + '_' + plottedT0 + '_fburial_det.png'
                        geniemap(sed_lon_edges, sed_lat_edges,  diff, diffcmap, difflevs, diffticklevs, diffextend, lower, upper, projdata ,cbartitle, difffilename, 'y', sed_lon, sed_lat, diffclevs, 0.75, 'n', 'none', 'none','png')
                        difffilecount += 1; difflnfile = 'difffile' + str(difffilecount) + '.png'
                        os.system('ln -s ' + difffilename + ' ' + difflnfile)
                    diffsavedfiles.append(difffilename)

        ########################################################################
        #                                  LATEX                               #
        ########################################################################

        if create_pdf_summary == 'y':
            sentence_to_use= ('Year plotted: ' + str(time2plot) + ' \linebreak '
            'Global SAT: ' + str(SAT_avg) + '$^\circ$C \linebreak '
            'Seafloor anoxic area: ' + str(np.round(seafloor_anoxia_abs/1E6/1E6,2)) + ' Mkm$^2$ \linebreak '
            'Seafloor anoxic area percentage: ' + str(np.round(seafloor_anoxia_per*100,2)) + ' \%')
            dopdf(savedfiles, 'file', filecount, exp + '_' + plottedT, sentence_to_use)

        if (do_diff == 'y'and globalcount == 1):
            diff_sentence_to_use= ('Year plotted: ' + str(time2plot) + ' \linebreak '
            'Global SAT: ' + str(SAT_avg_diff) + '$^\circ$C \linebreak '
            'Seafloor anoxic area: ' + str(np.round(seafloor_anoxia_abs_diff/1E6/1E6,2)) + ' Mkm$^2$ \linebreak '
            'Seafloor anoxic area percentage: ' + str(np.round(seafloor_anoxia_per_diff*100,2)) + ' \%')
            dopdf(diffsavedfiles, 'difffile', difffilecount, exp1 + '_' + plottedT1+ '_minus_' + exp0 + '_' + plottedT0, diff_sentence_to_use ) 
        timecount += 1
        globalcount += 1
    expcount += 1
    
# The End


