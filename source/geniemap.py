def geniemap(londata_edge, latdata_edge, zdata, cmapdata, levdata, tickdata, extendata, lowercolor, uppercolor, projdata, titledata, fnamedata, contourplotflag, londata, latdata, clevdata, lwdata, overlay_flag, overlay_zdata, overlay_cmapdata, outfmtdata):


    nlon = np.shape(zdata)[0]
    nlat = np.shape(zdata)[1]
    data_crs = ccrs.PlateCarree()
    landseamask = np.full((nlon, nlat), 1)
    landseamask = np.ma.masked_where(zdata.mask==False, landseamask)
    norm = BoundaryNorm(levdata, ncolors=cmapdata.N, clip=False)
    fig =  plt.figure(figsize=(figXsize, figYsize))
    ax = fig.add_subplot(111,projection=projdata)
    #ax = plt.axes([0.3, 0.1, 0.4, 0.8], projection=projdata)
    #ax.set_aspect('auto')
    ax = fig.gca()
    ax.set(aspect=2)
    ax.set_aspect('auto')
    cf = plt.pcolormesh(londata_edge, latdata_edge, zdata, transform=data_crs, cmap=cmapdata, norm=norm) # https://stackoverflow.com/questions/20678817/pyplot-pcolormesh-confused-when-alpha-not-1
    if contourplotflag == 'y':
        ct = plt.contour(londata, latdata, zdata, clevdata, cmap = None, colors='k',linewidths= lwdata,transform=ccrs.PlateCarree())
    if overlay_flag == 'y': # overlay another shading (used for SI)
        cf2 = plt.pcolormesh(londata_edge, latdata_edge, overlay_zdata, transform=data_crs, cmap=overlay_cmapdata)
    # upper-left annotation is the global sum over the whole land area
    if (plot_stepped_outline == 'y'): # continental outline
        stepped_coastline_cGENIE(londata_edge,latdata_edge, data_crs, landseamask,1.25)
        plt.pcolormesh(londata_edge, latdata_edge, landseamask, cmap=whitecmap, transform=data_crs)
    gl = ax.gridlines(crs=data_crs, draw_labels=False,
                  linewidth=0.5, color='k', alpha=0.2, linestyle='-', zorder=999)
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.xlines = True
    gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': plot_font_size, 'color': 'gray'}
    gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
    cb = fig.colorbar(cf, orientation='horizontal',extend=extendata,ticks=tickdata)
    cf.cmap.set_under(lowercolor)
    cf.cmap.set_over(uppercolor)
    cb.ax.tick_params(labelsize=plot_font_size)
    cb.ax.set_title(titledata, weight='normal', fontsize=plot_font_size)
    plt.tight_layout()

    if save_fig == 'y':
        fname = fnamedata
        plt.savefig(fname,format=outfmtdata, dpi=450)
    return()
