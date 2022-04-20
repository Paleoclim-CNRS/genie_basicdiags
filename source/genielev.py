def genielev(londata_edge, latdata_edge, zdata, ztdata, cmapdata, levdata, tickdata, extendata, lowercolor, uppercolor, projdata, titledata, fnamedata, contourplotflag, londata, latdata, clevdata):

    # zdata should have 3 dimsn: levs, lon, lat

    # update Dec 1 2019 :: accounts for different grids, i.e., empty levels in depth

    data_crs = ccrs.PlateCarree()
    norm = BoundaryNorm(levdata, ncolors=cmapdata.N, clip=False)
    nlevels = np.shape(zdata)[0]

    nonemptylevels = 0 # how many levels filled with data ?
    for k in np.arange(nlevels):
        ixx = np.argwhere(zdata[k,:,:].mask==False)
        if ixx.size != 0:
            nonemptylevels += 1
    
    nlevels = nonemptylevels

    # number of plots and figure size
    nbofplots = nlevels
    plotsperline = 4 # must be an integer
    nboflines = int(nbofplots/plotsperline)
    if nboflines == 0:
        nboflines += 1
    if nboflines*plotsperline < nbofplots:
        nboflines += 1
    if nbofplots >= plotsperline:
        nbofcol = plotsperline
    else:
        nbofcol = nbofplots
    nbofunusedsubplots = nboflines*nbofcol - nbofplots

    plot_font_size = 36

    # initializing figure
    fig, axes =  plt.subplots(figsize=(nbofcol*3.2, nboflines*2.075),nrows=nboflines, ncols=nbofcol)

    for k in np.arange(nlevels):

        if plot_stepped_outline == 'y':
            # building fake land-sea mask - required by steppe_coastline
            landsea_mask = np.where(zdata[k,:,:].mask==False,-1,1)
            land = np.ma.masked_where(landsea_mask==(-1),landsea_mask)

        # subplot
        ax = plt.subplot(nboflines,nbofcol,k+1,projection=projdata)
        ax = fig.gca()
        ax.set_aspect('auto')
        cf = plt.pcolormesh(londata_edge, latdata_edge, zdata[k,:,:], transform=data_crs,cmap=cmapdata, norm=norm)
        ct = plt.contour(londata, latdata, zdata[k,:,:], clevdata, cmap = None, colors='k',linewidths=0.75,transform=data_crs)
        if plot_stepped_outline == 'y':
            stepped_coastline_cGENIE(londata_edge,latdata_edge, ccrs.PlateCarree(), landsea_mask,1.5)
        if landflag == 'y':
            plt.pcolormesh(londata_edge, latdata_edge, land, cmap=whitecmap, transform=data_crs)
        plt.annotate('(' + str(k+1) + ') ' + str(np.round(ztdata[k],0)) + ' m', xy=(xlabel, ylabel), xycoords='axes fraction',fontsize=plot_font_size*0.5, color='k')
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                      linewidth=0.5, color='k', alpha=1, linestyle='-', zorder=999)
        gl.xlabels_top = False
        gl.ylabels_left = False
        gl.xlines = True
        gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': plot_font_size, 'color': 'gray'}
        gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
        #fig.subplots_adjust(left=0., bottom=-0.01, right=1., top=0.99, wspace=0.02, hspace=0.005)        # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.subplots_adjust.html
    fig.subplots_adjust(left=0., bottom=0.11, right=1., top=0.96, wspace=0.02, hspace=0.22)        # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.subplots_adjust.html
    cbar_ax = fig.add_axes([0.35, 0.04, 0.3, 0.02]) # list [x0, y0, width, height]
    cb = fig.colorbar(cf, cax=cbar_ax,orientation='horizontal',extend=extendata,ticks=tickdata)
    cf.cmap.set_under(lowercolor)
    cf.cmap.set_over(uppercolor)
    cb.ax.tick_params(labelsize=plot_font_size*0.5)
    cb.ax.set_title(titledata, weight='normal', fontsize=plot_font_size*0.5)

    # deleting unused subplots
    for i in range(-nbofunusedsubplots,0):
        axes[-1,i].set_axis_off()

    if save_fig == 'y':
        fname = fnamedata
        plt.savefig(fname,format='png', dpi=450)
    return()
