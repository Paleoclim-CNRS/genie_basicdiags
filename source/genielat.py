def genielat(latdata_edge, zdata_edge, data, cmapdata, levdata, tickdata, extendata, lowercolor, uppercolor, titledata, fnamedata, contourplotflag, latdata, zdata, clevdata, lwdata):
    nx = np.shape(data)[0]
    ny = np.shape(data)[1]
    mask = np.full((nx, ny), 1)
    mask = np.ma.masked_where(data.mask==False, mask)
    norm = BoundaryNorm(levdata, ncolors=cmapdata.N, clip=False)
    fig =  plt.figure(figsize=(figXsize, figYsize))
    ax = fig.add_subplot(111)
    ax = fig.gca()
    ax.set(aspect=2)
    ax.set_aspect('auto')
    cf = plt.pcolormesh(latdata_edge, zdata_edge, data, cmap=cmapdata, norm=norm) # https://stackoverflow.com/questions/20678817/pyplot-pcolormesh-confused-when-alpha-not-1
    if contourplotflag == 'y':
        ct = plt.contour(latdata, zdata, np.squeeze(data), clevdata, cmap = None, colors='k',linewidths= lwdata)
    # upper-left annotation is the global sum over the whole land area
    stepped_outline_cGENIE(latdata_edge,zdata_edge, mask,1.25)
    plt.pcolormesh(latdata_edge, zdata_edge, mask, cmap=whitecmap)
    plt.ylim(-90, 90)
    plt.ylim(-5.500, 0)
    major_yticks = np.arange(-5, 0+1E-3, 1)
    major_xticks = np.arange(-90, 90+1E-3, 30)
    ax.set_xticks(major_xticks)
    ax.set_yticks(major_yticks)
    ax.tick_params(labelsize = plot_font_size)
    ax.grid(which='both')
    cb = fig.colorbar(cf, orientation='horizontal',extend=extendata,ticks=tickdata)
    cf.cmap.set_under(lowercolor)
    cf.cmap.set_over(uppercolor)
    cb.ax.tick_params(labelsize=plot_font_size)
    cb.ax.set_title(titledata, weight='normal', fontsize=plot_font_size)

    if save_fig == 'y':
        fname = fnamedata
        plt.savefig(fname,format='png', dpi=450)
    return()
