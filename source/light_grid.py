#!/usr/bin/env python

def light_grid():
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

    gl = ax.gridlines(crs=projdata, draw_labels=False,
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

