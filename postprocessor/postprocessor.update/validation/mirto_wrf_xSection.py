import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from netCDF4 import Dataset
import wrf
from wrf import to_np, getvar, CoordPair, vertcross
from matplotlib.ticker import FormatStrFormatter
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import glob
import os
import sys
from argparse import ArgumentParser
from datetime import datetime, timedelta
from matplotlib.cm import ScalarMappable

__author__ = "Paolo Antonelli, Paolo Scaccia and Alessandra Valletti "
__copyright__ = "Copyright 2017, AdaptiveMeteo S.r.l."
__credits__ = ["Paolo Antonelli", "Paolo Scaccia", "Alessndra Valletti"]
__license__ = "--"
__version__ = "0.0.1"
__maintainer__ = "Paolo Antonelli"
__email__ = "paolo.antonelli@adaptivemeteo.com"
__status__ = "Development"


def get_wspd_colormap():

    import numpy as np
    from scipy import interpolate

    Source = np.asarray([[255,255,255], [246,246,248], [229,230,245], [201,216,236], [101,156,186], [33,157,74], [122,222,0], [245,225,1], [243,142,7], [244,80,37], [248,66,114], [212,61,212], [139,39,139], [71,20,71], [19,9,20], [43,43,43], [86,86,86], [126,126,126], [170,170,170], [213,213,213]], dtype=np.int)

    x = np.arange(0, Source.shape[0])

    fit = interpolate.interp1d(x, Source, axis=0)
    
    Target = fit(np.linspace(0, Source.shape[0]-1, 256))

    #print('Target {}'.format(Target))

    return Target/255

def get_pv_colormap():

    import numpy as np
    from scipy import interpolate

    Source = np.asarray([[71,17,99], [73,45,124], [73,45,124], [72,48,125], [74,61,132], [76,77,142], [82,98,152], [97,125,164], [122,154,182], [156,187,202], [196,217,224], [231,241,242], [252,253,253], [229,244,241], [194,231,221], [156,218,194], [127,210,166], [114,207,139], [115,208,113], [128,211,91], [149,216,69], [174,220,50], [201,224,30], [227,227,25]], dtype=np.int)

    x = np.arange(0, Source.shape[0])

    fit = interpolate.interp1d(x, Source, axis=0)

    Target = fit(np.linspace(0, Source.shape[0]-1, 256))

    #print('Target {}'.format(Target))

    return Target/255

def get_w_data(ncfile):

    # Extract the model height and wind speed
    z = wrf.g_geoht.get_height(ncfile, timeidx=0, method='cat', squeeze=True, cache=None, meta=True, _key=None, msl=True, units='m')
    #w = getvar(ncfile, "wa")
    ter = getvar(ncfile, "ter", timeidx=-1)
    w = wrf.g_wind.get_w_destag(ncfile, timeidx=0, method='cat', squeeze=True, cache=None, meta=True, _key=None, units='m s-1')
    print("W extremes:{} {}".format(np.max(to_np(w)),np.min(to_np(w))))

    return w, z, ter

def get_wspd_data(ncfile):

    # Extract the model height and wind speed
    #z = getvar(ncfile, "z")
    z = wrf.g_geoht.get_height(ncfile, timeidx=0, method='cat', squeeze=True, cache=None, meta=True, _key=None, msl=True, units='m')
    print("Z size: {}".format(z.shape))
    wspd_dir = wrf.g_wind.get_destag_wspd_wdir(ncfile, timeidx=0, method='cat', squeeze=True, cache=None, meta=True, _key=None, units='m s-1') 
    wspd = wspd_dir[0,:,:,:]
    wdir = wspd_dir[1,:,:,:]
    print("WSPD size: {}".format(wspd.shape))

    return wspd, wdir, z

def get_pvo_data(ncfile):

    # Extract the model height and Potential Vorticity
    #z = getvar(ncfile, "z")
    z = wrf.g_geoht.get_height(ncfile, timeidx=0, method='cat', squeeze=True, cache=None, meta=True, _key=None, msl=True, units='m')
    print("Z size: {}".format(z.shape))
    pvo = wrf.g_vorticity.get_pvo(ncfile, timeidx=0, method='cat', squeeze=True, cache=None, meta=True, _key=None)
    print("PV size: {}".format(pvo.shape))

    return pvo, z

def get_rh_data(ncfile):

    # Extract the model height and Potential Vorticity
    #z = getvar(ncfile, "z")
    z = wrf.g_geoht.get_height(ncfile, timeidx=0, method='cat', squeeze=True, cache=None, meta=True, _key=None, msl=True, units='m')
    print("Z size: {}".format(z.shape))
    ter = getvar(ncfile, "ter", timeidx=-1)
    rh = wrf.g_rh.get_rh(ncfile, timeidx=0, method='cat', squeeze=True, cache=None, meta=True, _key=None)

    return rh, z, ter

def make_wspd_xsection(ncfile,wspd,z,lat1,lon1,lat2,lon2,runTypeStr,tltStr):

    # Create the start point and end point for the cross section
    start_point = CoordPair(lat=lat1, lon=lon1)
    end_point = CoordPair(lat=lat2, lon=lon2)

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section.
    #print("WSPD: {}".format(to_np(wspd)))
    wspd_cross = vertcross(wspd, z, wrfin=ncfile, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)

    # Create the figure
    fig = plt.figure(figsize=(12,4))
    ax = plt.axes()
    #ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))

    # Make the contour plot
    vmin = 0.0
    vmax = 15.0
    #print("WSPD CROSS: {}".format(to_np(wspd_cross)))

    #wspd_contours = ax.contourf(to_np(wspd_cross), 40, levels = np.arange(vmin,vmax,(vmax-vmin)/50), cmap=ListedColormap(get_wspd_colormap()))
    wspd_contours = ax.contourf(to_np(wspd_cross), 40, levels = np.arange(vmin,vmax,(vmax-vmin)/50), cmap='jet',extend='both')

    # Set the x-ticks to use latitude and longitude labels.
    coord_pairs = to_np(wspd_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in to_np(coord_pairs)]
    #print("HTick: {}".format(x_ticks))
    #print("XLABELS: {}".format(x_labels))
    ax.set_xticks(x_ticks[::4])
    ax.set_xticklabels(x_labels[::4], rotation=45, fontsize=8)

    # Set the y-ticks to be height.
    vert_vals = to_np(wspd_cross.coords["vertical"])
    v_ticks = np.arange(vert_vals.shape[0])
    #print("VTick: {}".format(v_ticks))
    #print("VERT VALS: {}".format(vert_vals))

    ax.set_yticks(v_ticks[0:31:3])
    ax.set_yticklabels([str(round(float(label)/1000, 0)) for label in vert_vals[0:31:3]], fontsize=8)

    ax.set_ylim(ymin=1, ymax=31)

    # Set the x-axis and  y-axis labels
    ax.set_xlabel("Latitude, Longitude", fontsize=12)
    ax.set_ylabel("Height (km)", fontsize=12)

    tltStr_row1 = "Vertical Cross Section of Wind Speed [m/s]\n"
    tltStr_row2 = "{} {}".format(runTypeStr,tltStr)

    plt.title(tltStr_row1 + tltStr_row2)
    # Add the color bar
    aspect = 20
    pad_fraction = 0.5

    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=0.5)
    clabel = "Wind Speed [m/s]"
    #cbar = plt.colorbar(ScalarMappable(norm=wspd_contours.norm,cmap=wspd_contours.cmap),cax=cax,label=clabel,format='%2d',orientation="horizontal")
    cbar = plt.colorbar(ScalarMappable(norm=wspd_contours.norm,cmap=wspd_contours.cmap),cax=cax,label=clabel,format='%2d',orientation="vertical")
    cbar.set_ticks(np.arange(0,38,2))

    return plt

def make_w_xsection(ncfile,w,pvo,z,ter,lat1,lon1,lat2,lon2,runTypeStr,tltStr):

    # Create the start point and end point for the cross section
    start_point = CoordPair(lat=lat1, lon=lon1)
    end_point = CoordPair(lat=lat2, lon=lon2)
    
    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section.
    #print("W: {}".format(to_np(w)))
    w_cross = vertcross(w, z, wrfin=ncfile, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)
    pvo_cross = vertcross(pvo, z, wrfin=ncfile, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)


    xs = np.arange(0, w_cross.shape[-1], 1)
    print('XS SHAPE {}'.format(xs.shape))

    # Create the figure
    fig = plt.figure(figsize=(12,4))
    ax = plt.axes()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    # Make the contour plot
    vmin = -0.1
    vmax = 0.1
    #print("W CROSS: {}".format(to_np(w_cross)))

    w_contours = ax.contourf(to_np(w_cross), 10, levels = np.arange(vmin,vmax,(vmax-vmin)/10), cmap='jet',extend='both')
    #w_contours = ax.contourf(to_np(w_cross), 40, cmap='jet',extend='both')
    #w_contours = ax.contourf(to_np(w_cross), 40, cmap='jet',extend='both')
   
    # Add the PVO contours
    levels = np.arange(-5.0, 16.0, 2.)
    contours = ax.contour(pvo_cross, levels=levels, colors="black", linewidths=0.5)
    ax.clabel(contours, inline=True, fontsize=5, fmt = '%.1f')
 
    # Set the x-ticks to use latitude and longitude labels.
    coord_pairs = to_np(w_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in to_np(coord_pairs)]
    #print("HTick: {}".format(x_ticks))
    #print("XLABELS: {}".format(x_labels))
    ax.set_xticks(x_ticks[::4])
    ax.set_xticklabels(x_labels[::4], rotation=45, fontsize=8)

    # Set the y-ticks to be height.
    vert_vals = to_np(w_cross.coords["vertical"])
    v_ticks = np.arange(vert_vals.shape[0])
    #print("VTick: {}".format(v_ticks))
    #print("VERT VALS: {}".format(vert_vals))

    ax.set_yticks(v_ticks[0:39:3])
    ax.set_yticklabels([str(round(float(label)/1000, 0)) for label in vert_vals[0:39:3]], fontsize=8)

    ax.set_ylim(ymin=1, ymax=39)


    # Set the x-axis and  y-axis labels
    ax.set_xlabel("Latitude, Longitude", fontsize=12)
    ax.set_ylabel("Height (km)", fontsize=12)
 
    tltStr_row1 = "Vertical Cross Section of Vertical Wind [m/s] + PVO lines\n"
    tltStr_row2 = "{} {}".format(runTypeStr,tltStr)

    #Add Terrain
    # Get the terrain heights along the cross section line
    ter_line = wrf.interpline(ter, wrfin=ncfile, start_point=start_point, end_point=end_point)

    # Fill in the mountain area
    ht_fill = ax.fill_between(xs, 1, to_np(ter_line)/1000.0, facecolor="saddlebrown")

    plt.title(tltStr_row1 + tltStr_row2)
    # Add the color bar
    aspect = 20
    pad_fraction = 0.5

    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=0.5)
    clabel = "Vertical Wind [m/s]"
    #cbar = plt.colorbar(ScalarMappable(norm=wspd_contours.norm,cmap=wspd_contours.cmap),cax=cax,label=clabel,format='%2d',orientation="horizontal")
    cbar = plt.colorbar(ScalarMappable(norm=w_contours.norm,cmap=w_contours.cmap),cax=cax,label=clabel,format='%.2f',orientation="vertical")
    cbar.set_ticks(np.arange(-.1,.1,.02))

    return plt

def make_wdir_xsection(ncfile,wdir,z,lat1,lon1,lat2,lon2,runTypeStr,tltStr):

    # Create the start point and end point for the cross section
    start_point = CoordPair(lat=lat1, lon=lon1)
    end_point = CoordPair(lat=lat2, lon=lon2)

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section.
    #print("WSPD: {}".format(to_np(wdir)))
    wdir_cross = vertcross(wdir, z, wrfin=ncfile, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)

    # Create the figure
    fig = plt.figure(figsize=(12,4))
    ax = plt.axes()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    # Make the contour plot
    vmin = 0.0
    vmax = 360.0
    #print("WSPD CROSS: {}".format(to_np(wdir_cross)))

    #wdir_contours = ax.contourf(to_np(wdir_cross), 40, levels = np.arange(vmin,vmax,(vmax-vmin)/50), cmap=ListedColormap(get_wdir_colormap()))
    wdir_contours = ax.contourf(to_np(wdir_cross), 40, levels = np.arange(vmin,vmax,(vmax-vmin)/50), cmap='hsv',extend='both')

    # Set the x-ticks to use latitude and longitude labels.
    coord_pairs = to_np(wdir_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in to_np(coord_pairs)]
    #print("HTick: {}".format(x_ticks))
    #print("XLABELS: {}".format(x_labels))
    ax.set_xticks(x_ticks[::4])
    ax.set_xticklabels(x_labels[::4], rotation=45, fontsize=8)

    vert_vals = to_np(wdir_cross.coords["vertical"])
    v_ticks = np.arange(vert_vals.shape[0])
    #print("VTick: {}".format(v_ticks))
    #print("VERT VALS: {}".format(vert_vals))

    ax.set_yticks(v_ticks[0:31:3])
    ax.set_yticklabels([str(round(float(label)/1000, 0)) for label in vert_vals[0:31:3]], fontsize=8)

    ax.set_ylim(ymin=1, ymax=31)

    # Set the x-axis and  y-axis labels
    ax.set_xlabel("Latitude, Longitude", fontsize=12)
    ax.set_ylabel("Height (km)", fontsize=12)

    tltStr_row1 = "Vertical Cross Section of Wind Dir [Degrees]\n"
    tltStr_row2 = "{} {}".format(runTypeStr,tltStr)

    plt.title(tltStr_row1 + tltStr_row2)
    # Add the color bar
    aspect = 20
    pad_fraction = 0.5

    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=0.5)
    clabel = "Wind Dir [Degrees]"
    #cbar = plt.colorbar(ScalarMappable(norm=wdir_contours.norm,cmap=wdir_contours.cmap),cax=cax,label=clabel,format='%2d',orientation="horizontal")
    cbar = plt.colorbar(ScalarMappable(norm=wdir_contours.norm,cmap=wdir_contours.cmap),cax=cax,label=clabel,format='%2d',orientation="vertical")
    cbar.set_ticks(np.arange(0,360,30))

    return plt

def make_pvo_xsection(ncfile,pvo,z,lat1,lon1,lat2,lon2,runTypeStr,tltStr):

    # Create the start point and end point for the cross section
    #start_point = CoordPair(lat=70.0, lon=15.0)
    #end_point = CoordPair(lat=85.47, lon=15.47)
    start_point = CoordPair(lat=lat1, lon=lon1)
    end_point = CoordPair(lat=lat2, lon=lon2)

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section.
    pvo_cross = vertcross(pvo, z, wrfin=ncfile, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)

    # Create the figure
    fig = plt.figure(figsize=(12,6))
    ax = plt.axes()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    # Make the contour plot
    vmin = -5.0
    vmax = 16.0
    #pvo_contours = ax.contourf(to_np(pvo_cross), 40, levels = np.arange(vmin,vmax,(vmax-vmin)/50), cmap=ListedColormap(get_pv_colormap()))
    pvo_contours = ax.contourf(to_np(pvo_cross), 40, levels = np.arange(vmin,vmax,(vmax-vmin)/50), cmap='jet')

    # Add the color bar
    aspect = 20
    pad_fraction = 0.5

    # Set the x-ticks to use latitude and longitude labels.
    coord_pairs = to_np(pvo_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in to_np(coord_pairs)]
    ax.set_xticks(x_ticks[::10])
    ax.set_xticklabels(x_labels[::10], rotation=45, fontsize=8)

    # Set the y-ticks to be height.
    vert_vals = to_np(pvo_cross.coords["vertical"])
    v_ticks = np.arange(vert_vals.shape[0])
    #print("VTick: {}".format(v_ticks))
    #print("VERT VALS: {}".format(vert_vals))

    ax.set_yticks(v_ticks[0:39:3])
    ax.set_yticklabels([str(round(float(label)/1000, 0)) for label in vert_vals[0:39:3]], fontsize=8)

    ax.set_ylim(ymin=1, ymax=39)


    # Set the x-axis and  y-axis labels
    ax.set_xlabel("Latitude, Longitude", fontsize=12)
    ax.set_ylabel("Height (km)", fontsize=12)

    tltStr_row1 = "Vertical Cross Section of Potential Vorticity  $[m^{2}s^{-1}KKg^{-1}]$ \n"
    tltStr_row2 = "{} {}".format(runTypeStr,tltStr)


    #plt.title("CONV+TR_CLR+PHYS_CLD Vertical Cross Section of Wind Speed (kt)")
    plt.title(tltStr_row1 + tltStr_row2)

    # Add the color bar
    aspect = 20
    pad_fraction = 0.5

    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    #cax = divider.append_axes("bottom", size=width, pad=1.0)
    cax = divider.append_axes("bottom", size=width, pad=1.0)
    clabel = r"Potential Vorticity $[m^{2}s^{-1}KKg^{-1}]$"
    cbar = plt.colorbar(ScalarMappable(norm=pvo_contours.norm,cmap=pvo_contours.cmap),cax=cax,label=clabel,format='%2d',orientation="horizontal")
    cbar.set_ticks(np.arange(-24,24,2))

    return plt

def make_rh_xsection(ncfile,rh,pvo,z,ter,lat1,lon1,lat2,lon2,runTypeStr,tltStr):

    # Create the start point and end point for the cross section
    #start_point = CoordPair(lat=70.0, lon=15.0)
    #end_point = CoordPair(lat=85.47, lon=15.47)
    start_point = CoordPair(lat=lat1, lon=lon1)
    end_point = CoordPair(lat=lat2, lon=lon2)

    # Compute the vertical cross-section interpolation.  Also, include the
    # lat/lon points along the cross-section.
    rh_cross = vertcross(rh, z, wrfin=ncfile, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)
    pvo_cross = vertcross(pvo, z, wrfin=ncfile, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)


    # Create the figure
    fig = plt.figure(figsize=(12,6))
    ax = plt.axes()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


    # Make the contour plot
    vmin = 0.0
    vmax = 105.0
    #rh_contours = ax.contourf(to_np(rh_cross), 40, levels = np.arange(vmin,vmax,(vmax-vmin)/50), cmap=ListedColormap(get_pv_colormap()))
    rh_contours = ax.contourf(to_np(rh_cross), 40, levels = np.arange(vmin,vmax,(vmax-vmin)/50), cmap='jet')

    xs = np.arange(0, rh_cross.shape[-1], 1)
    print('XS SHAPE {}'.format(xs.shape))

    # Add the PVO contours
    #levels = np.arange(-5.0, 16.0, 2.)
    #contours = ax.contour(pvo_cross, levels=levels, colors="black", linewidths=0.5)
    #ax.clabel(contours, inline=True, fontsize=5, fmt = '%.1f')

    # Add the color bar
    aspect = 20
    pad_fraction = 0.5

    # Set the x-ticks to use latitude and longitude labels.
    coord_pairs = to_np(rh_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in to_np(coord_pairs)]
    ax.set_xticks(x_ticks[::10])
    ax.set_xticklabels(x_labels[::10], rotation=45, fontsize=8)

    # Set the y-ticks to be height.
    vert_vals = to_np(rh_cross.coords["vertical"])
    v_ticks = np.arange(vert_vals.shape[0])
    #print("VTick: {}".format(v_ticks))
    #print("VERT VALS: {}".format(vert_vals))

    ax.set_yticks(v_ticks[0:39:3])
    ax.set_yticklabels([str(round(float(label)/1000, 0)) for label in vert_vals[0:39:3]], fontsize=8)

    ax.set_ylim(ymin=1, ymax=39)


    # Set the x-axis and  y-axis labels
    ax.set_xlabel("Latitude, Longitude", fontsize=12)
    ax.set_ylabel("Height (km)", fontsize=12)

    tltStr_row1 = "Vertical Cross Section of Relative Humidity [%] \n"
    tltStr_row2 = "{} {}".format(runTypeStr,tltStr)

    #Add Terrain
    # Get the terrain heights along the cross section line
    ter_line = wrf.interpline(ter, wrfin=ncfile, start_point=start_point, end_point=end_point)

    # Fill in the mountain area
    ht_fill = ax.fill_between(xs, 0, to_np(ter_line)/1000.0, facecolor="saddlebrown")

    #plt.title("CONV+TR_CLR+PHYS_CLD Vertical Cross Section of Wind Speed (kt)")
    plt.title(tltStr_row1 + tltStr_row2)


    # Add the color bar
    aspect = 20
    pad_fraction = 0.5

    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    #cax = divider.append_axes("bottom", size=width, pad=1.0)
    cax = divider.append_axes("bottom", size=width, pad=1.0)
    clabel = r"Relative Humidity [%]"
    cbar = plt.colorbar(ScalarMappable(norm=rh_contours.norm,cmap=rh_contours.cmap),cax=cax,label=clabel,format='%2d',orientation="horizontal")
    cbar.set_ticks(np.arange(0,100,10))

    return plt

#####################################   MAIN  ##################
def main(argv):

    parser = ArgumentParser()
    parser.add_argument('--outdir', '-o',default='./',
                    help='Output plot directory')
    parser.add_argument('--wrf_dir',required=True,default='/galaxy/wrf/arctic/2022/con_tr_cloudy/',
                    help='Output plot directory')
    parser.add_argument('--assimilation_h','-ah', default=None, type=str,
                    help='Starting time of Forecating Cycle ex. 2022080812')
    parser.add_argument('--forecast_date_time','-fdt', default=None, type=str,
                    help='Validating forecast date and time in given forecasting cycle ex. 2022-08-07_00:00:00')
    parser.add_argument('--start_point','-stp', default=None, nargs='+', type=float,
                    help='Array with starting point lat lon')
    parser.add_argument('--end_point','-enp', default=None, nargs='+', type=float, 
                    help='Array with ending point lat lon')
  
    argv = parser.parse_args()

    lat1, lon1 = argv.start_point
    lat2, lon2 = argv.end_point


    #build selected wrf filename
    wrf_f_fmt = '%Y-%m-%d_%H:%M:%S'
    wrf_o_fmt = 'wrfout_d01_%Y-%m-%d_%H:%M:%S'
    wrf_tit_fmt = 'Date %Y-%m-%d time %H:%M UTC'

    if 'conv_tr_cloudy' in argv.wrf_dir:
        runTypeStr = 'CONV+TR_CLR+PHY_CLD'
    elif  'conv_tr'  in argv.wrf_dir:
        runTypeStr = 'CONV+TR_CLR'
    elif  'conv_cloudy'  in argv.wrf_dir:
        runTypeStr = 'CONV+PHY_CLD'
    elif  'ecmwf'  in argv.wrf_dir:
        runTypeStr = 'ECMWF'
    else:
        runTypeStr = 'CONV'

    if argv.assimilation_h is not None and argv.forecast_date_time is not None:
        wrffilenames = glob.glob(argv.wrf_dir + '/' + argv.assimilation_h + '/' + datetime.strptime(argv.forecast_date_time, wrf_f_fmt).strftime(wrf_o_fmt))
        print("WRF files: {}".format(wrffilenames))
    elif argv.assimilation_h is not None:
        wrffilenames = glob.glob(argv.wrf_dir + '/' + argv.assimilation_h + '/wrfout*')
    else: 
        wrfdirs = glob.glob(argv.wrf_dir + '/202*')
        #print("WRF directories: {}".format(wrfdirs))
        wrffilenames = []
        for dire in wrfdirs:
            #print('Derectory: {}',format(dire))
            wrffilenames.append(glob.glob(dire + '/wrfout*'))
            #print("WRF files: {}".format(wrffilenames))


    if len(wrffilenames) > 1:
        flat_wrffilenames = [item for sublist in wrffilenames for item in sublist]
    else: 
        flat_wrffilenames = wrffilenames

    print("WRF flat files: {}".format(flat_wrffilenames))
    for filename in flat_wrffilenames:
        forecast = os.path.basename(filename)
        print("Forecst Time and Date: {}".format(forecast))
        plotTitleStr = datetime.strptime(forecast, wrf_o_fmt).strftime(wrf_tit_fmt)


        #filename = "/galaxy/wrf/arctic/2022/con_tr_cloudy/2022080812/wrfout_d01_2022-08-08_12:00:00"
        #filename = "/galaxy/wrf/arctic/2022/conv/2022080812/wrfout_d01_2022-08-08_12:00:00"

        ncfile = Dataset(filename)

        wspd, wdir, z = get_wspd_data(ncfile)
        plt = make_wspd_xsection(ncfile,wspd,z,lat1,lon1,lat2,lon2,runTypeStr,plotTitleStr)

        outfilename = forecast.replace('.nc','.png')
        outfile = argv.outdir + '/WSPDxSection_' + outfilename
        plt.savefig(outfile,dpi=300,bbox_inches='tight')
        plt.close('all')


        pvo, z = get_pvo_data(ncfile)

        plt = make_pvo_xsection(ncfile,pvo,z,lat1,lon1,lat2,lon2,runTypeStr,plotTitleStr)
        outfilename = forecast.replace('.nc','.png')
        outfile = argv.outdir + '/PVxSection_' + outfilename
        plt.savefig(outfile,dpi=300,bbox_inches='tight')
        plt.close('all')

        rh, z, ter = get_rh_data(ncfile)

        plt = make_rh_xsection(ncfile,rh,pvo,z,ter,lat1,lon1,lat2,lon2,runTypeStr,plotTitleStr)
        outfilename = forecast.replace('.nc','.png')
        outfile = argv.outdir + '/RHxSection_' + outfilename
        plt.savefig(outfile,dpi=300,bbox_inches='tight')
        plt.close('all')

        w, z, ter = get_w_data(ncfile)
        plt = make_w_xsection(ncfile,w,pvo,z,ter,lat1,lon1,lat2,lon2,runTypeStr,plotTitleStr)

        outfilename = forecast.replace('.nc','.png')
        outfile = argv.outdir + '/WxSection_' + outfilename
        plt.savefig(outfile,dpi=300,bbox_inches='tight')
        plt.close('all')

        plt = make_wdir_xsection(ncfile,wdir,z,lat1,lon1,lat2,lon2,runTypeStr,plotTitleStr)

        outfilename = forecast.replace('.nc','.png')
        outfile = argv.outdir + '/WDIRxSection_' + outfilename
        plt.savefig(outfile,dpi=300,bbox_inches='tight')
        plt.close('all')



    return


#########################
if __name__=="__main__":
    main(sys.argv[1:])

