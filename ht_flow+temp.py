# -*- coding: utf-8 -*-
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd
import math
import subprocess
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from netCDF4 import Dataset
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm
import numpy.ma as ma
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator

from itertools import product

class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))
"""
Imaging HT results. 
"""
Workdir = '../ht-oiso_output/'
Workdir += 'perm_expexp_-16_8-11_8_zh300_spx1_200x320_irr-20200830'
Workdir += '/'


nx = 200
ny = 320

xmax = 30.
ymax=12.

plot_log = False
plot_log = True

dx = xmax/nx
dy = ymax/ny

x_edges =np.loadtxt(Workdir+"x_edges.txt")/1e3
y_edges =np.loadtxt(Workdir+"y_edges.txt")/1e3
x =np.loadtxt(Workdir+"x.txt")/1e3
y =np.loadtxt(Workdir+"y.txt")/1e3
dataU=np.loadtxt(Workdir+"vmx.txt")
dataV=np.loadtxt(Workdir+"vmy.txt")
dataP=np.loadtxt(Workdir+"pres.txt")
data=np.loadtxt(Workdir+"cd2.txt")
dataQx=np.loadtxt(Workdir+"qx.txt")
dataQy=np.loadtxt(Workdir+"qy.txt")
rho=np.loadtxt(Workdir+"rho.txt")
temp = np.loadtxt(Workdir+"temp.txt")

# this is the case without vabs file
yr2sec=60.*60.*24.*365.25
poro=0.05
vmx = dataU/rho*yr2sec/poro
vmy = dataV/rho*yr2sec/poro

vabs = np.sqrt(vmx*vmx+vmy*vmy)
##vabs=np.loadtxt(Workdir+"vabs.txt")

#####

vabs_max = np.log10(vabs)[data!=0].max()
vabs_min = np.log10(vabs)[data!=0].min()
vabs[data==0]=1e-10

obstacle_mask = False

no_background = False
##no_background = True

chk_calc = True
chk_calc = False

datam = ma.ones((ny,nx),dtype=np.float)

datam[:,:]=dataQx[:,:]
if chk_calc:
    datam[:,:]=dataQx[1:241,:]


##dataQy[data==0]=0
##dataQx[data==0]=0
datam[:,:]=dataQy[:,:]
##datam[:,:]=dataP[:,:]
obstacle_mask = True

if chk_calc:
    datam[:,:]=dataQx[1:241,:]

if obstacle_mask:
    datam[data==0]=-1e100
    datam=ma.masked_values(datam,-1e100)

plt.rcParams['font.family'] = 'Arial' 
plt.rcParams['font.size'] = 20

linewidth = 1.5

plt.rcParams['axes.linewidth'] = linewidth

plt.rcParams['xtick.major.width'] = linewidth
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.minor.width'] = linewidth
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.major.pad'] = 8

plt.rcParams['ytick.major.width'] = linewidth
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.minor.width'] = linewidth
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.major.pad'] = 8

plt.rcParams['axes.labelpad'] = 8

plt.rcParams['xtick.direction']='out'
plt.rcParams['ytick.direction']='out'

plt.tick_params(top=True)
plt.tick_params(right=True)

##cmap = cm.jet
##cmap = cm.inferno
##cmap = cm.plasma

figsize = (10,8)

##fig = plt.figure(figsize=figsize)

fig = plt.figure(figsize=figsize)

fs = 14

nx = 2
ny = 2

cmap = cm.gnuplot2
##cmap = cm.viridis
##cmap = cm.Pastel1
cmap = cm.terrain
cmap = cm.tab20c
cmap = cm.nipy_spectral_r
cmap = cm.binary_r
cmap=cm.jet
cmap.set_under('0.9')
##cmap = cm.Pastel1_r
##cmap = cm.cool 

levels = np.array([-1.8,-1.2,-0.6,0.,1.2,1.8,2.4,3.0,3.6])
levels = np.linspace(datam.min()*1e3,datam.max()*1e3,101)
levels = np.linspace(-4,-2,15)

levels1 = np.array([2,10, 30, 60, 100, 300, 600, 1200])
levels2 = [3, 10, 30, 100, 300, 1000]
##levels = levels*0.1
norm = MidpointNormalize(midpoint=0)

##x = np.linspace(0,xmax,nx,endpoint=True)
##y=np.linspace(0,ymax,ny,endpoint=True)
lons, lats = np.meshgrid(x_edges,y_edges)
x2, y2 = (lons, lats)
##x2 = np.linspace(0-dx/2.,xmax+dx/2.,nx+1,endpoint=True)
##y2=np.linspace(0-dy/2.,ymax+dy/2.,ny+1,endpoint=True)
##plt.axes().set_aspect('equal')

speed = np.sqrt(dataU*dataU + dataV*dataV)

##Q=plt.pcolormesh(x, y, datam)
ax1 = plt.subplot2grid((ny,nx), (0,0))

Qr=ax1.pcolor(x2,y2,np.log10(vabs),cmap=cmap,vmax=vabs_max
               ,vmin=vabs_min)
##cbar=plt.colorbar(Q,format='%.0f'
##                  )
ax1.contour(x,y,np.log10(datam), 21,corner_mask=True
            ,colors='w'
            ,linewidths=0.6
            ,linestyles='solid'
            ,levels=levels
            )

ax1.set_xlim(0,xmax)
ax1.set_ylim(0,ymax)
ax1.xaxis.set_minor_locator(MultipleLocator(5))
ax1.yaxis.set_minor_locator(MultipleLocator(1))
ax1.invert_yaxis()


ax2 = plt.subplot2grid((ny,nx), (0,1))

Q=ax2.pcolor(x2,y2,temp,cmap=cmap
                  ,norm=LogNorm()
             )
##cbar=plt.colorbar(Q,format='%.0f'
##                  ,norm=LogNorm()
##                  ,ticks=levels2
##                  )
cs = ax2.contour(x,y,temp, colors='k',linestyles='solid', levels=levels1)
# plt.clabel(cs, fontsize=12,inline=False,fmt='%.0f',colors='k'
           # )
ax2.set_ylim(0,ymax)
ax2.xaxis.set_minor_locator(MultipleLocator(5))
ax2.yaxis.set_minor_locator(MultipleLocator(1))
ax2.invert_yaxis()
ax2.set_xlim(0,xmax)

ax3 = plt.subplot2grid((ny,nx), (1,0))

Qr=ax3.pcolor(x2,y2,np.log10(vabs),cmap=cmap,vmax=vabs_max
               ,vmin=vabs_min)
##cbar=plt.colorbar(Q,format='%.0f'
##                  )
ax3.contour(x,y,np.log10(datam), 21,corner_mask=True
            ,colors='w'
            ,linewidths=0.6
            ,linestyles='solid'
            ,levels=levels
            )

ax3.set_xlim(xmax*1e-4,xmax)
ax3.set_ylim(1e-3,ymax)
ax3.set_yscale('log')
ax3.set_xscale('log')
locmaj = ticker.LogLocator(base=10,numticks=4) 
ax3.xaxis.set_major_locator(locmaj)
locmin = ticker.LogLocator(base=10.0
                           ,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
                           ,numticks=4)
ax3.xaxis.set_minor_locator(locmin)
ax3.xaxis.set_minor_formatter(ticker.NullFormatter())
ax3.invert_yaxis()


ax4 = plt.subplot2grid((ny,nx), (1,1))

Q=ax4.pcolor(x2,y2,temp,cmap=cmap
                  ,norm=LogNorm()
             )
##cbar=plt.colorbar(Q,format='%.0f'
##                  ,norm=LogNorm()
##                  ,ticks=levels2
##                  )
cs = ax4.contour(x,y,temp, colors='k',linestyles='solid', levels=levels1)
# plt.clabel(cs, fontsize=12,inline=False,fmt='%.0f',colors='k'
           # )
##xminors = np.log10(([k*m for k,m in product([2,3,4,5,6,7,8,9],[1e-3,100])]))
ax4.set_xlim(xmax*1e-4,xmax)
ax4.set_ylim(1e-3,ymax)
ax4.set_yscale('log')
ax4.set_xscale('log')
locmaj = ticker.LogLocator(base=10,numticks=4) 
ax4.xaxis.set_major_locator(locmaj)
locmin = ticker.LogLocator(base=10.0
                           ,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
                           ,numticks=4)
ax4.xaxis.set_minor_locator(locmin)
ax4.xaxis.set_minor_formatter(ticker.NullFormatter())
##ax.xaxis.set_minor_locator(ticker.FixedLocator(xminors))
ax4.invert_yaxis()

cbaxes = fig.add_axes([0.4, 0.39, 0.01, 0.3]) 
fig.colorbar(Qr, ax=[ax1,ax3],cax = cbaxes
             ,format='%.0f'
             )
cbaxes = fig.add_axes([0.85, 0.39, 0.01, 0.3]) 
fig.colorbar(Q, ax=[ax2,ax4],cax = cbaxes
             ,format='%.0f'
                  ,norm=LogNorm()
                  ,ticks=levels2
             )
fig.subplots_adjust(left=0.10,bottom=0.20,wspace=0.7,hspace=0.2,right=0.82) 
outfilename = Workdir+"flow_temp_v2.svg"
plt.savefig(outfilename, transparent=True)
subprocess.call('"C:\Program Files\Inkscape\inkscape.exe" -z -f ' \
                + outfilename + ' --export-emf '+outfilename+\
                '.emf',shell=True)
plt.show()
plt.clf()
plt.close()


