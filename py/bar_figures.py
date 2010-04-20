###############################################################################
#   bar_figures.py: make nice figures for the bar
###############################################################################
import os, os.path
import sys
import cPickle as pickle
import math as m
import scipy as sc
import bovy_plot as plot
from matplotlib import pyplot
from calc_veldist_2d import calc_veldist_2d
_degtorad= m.pi/180.
_XWIDTH= 1.8*8/10/1.8
_YWIDTH= 1.15*8/10/1.8
def veldist_2d_Rphi(plotfilename,nx=10,ny=6,dx=_XWIDTH/20.,dy=_YWIDTH/20.,
                    nsx=2,nsy=2,ngrid=51,rrange=[0.7,1.3],
                    phirange=[-m.pi/2.,m.pi/2.],
                    saveDir='../bar/2d/'):
    """
    NAME:
       veldist_2d_Rphi
    PURPOSE:
       plot how the velocity distribution changes as a function of R and phi
    INPUT:
       nx - number of plots in the x-direction
       ny - number of plots in the y direction
       dx - x-spacing
       dy - y-spacing
       nsx, nsy - number of subplots in the middle plot
       ngrid - number of gridpoints to evaluate the density on
       rrange - range of Galactocentric radii to consider
       phirange - range of Galactic azimuths to consider
       saveDir - directory to save the pickles in
    OUTPUT:
       plot!
    HISTORY:
       2010-04-19 - Written - Bovy (NYU)
    """
    levels= sc.array([2,6,12,21,33,50,68,80,90,95,99,99.9])/100.
    cntrcolors= ['w' for ii in range(len(levels)) if levels[ii] <= .5]
    cntrcolors+= ['k' for ii in range(len(levels)) if levels[ii] > .5]

    ulinspace= (-0.9,0.9,ngrid)
    vlinspace= (-.7,.45,ngrid)

    picklebasename= '2d_%i_%i_%i_%i_%i_%.1f_%.1f_%.1f_%.1f' % (nx,ny,nsx,nsy,ngrid,rrange[0],rrange[1],phirange[0],phirange[1])
    if not os.path.exists(saveDir):
        os.mkdir(saveDir)
    left, bottom = 0.1, 0.1
    width= nx*_XWIDTH+(nx-1)*dx+2*left
    height= ny*_YWIDTH+(ny-1)*dy+2*bottom
    plot.bovy_print(fig_width=width,fig_height=height,
                    xtick_major_size=2.,ytick_major_size=2.,
                    xtick_minor_size=0.,ytick_minor_size=0.)
    fig= pyplot.figure()
    for ii in range(nx):
        for jj in range(ny):
            #Middle plot
            if ii == nx/2-1 and jj == ny/2-1:
                thisax= fig.add_axes([(left+ii*(_XWIDTH+dx))/width,
                                      (bottom+jj*(_YWIDTH+dy))/height,
                                      (nsx*_XWIDTH+(nsx-1)*dx)/width,
                                      (nsy*_YWIDTH+(nsy-1)*dy)/height])
                thisR= (rrange[0]+(rrange[1]-rrange[0])/
                        (ny*_YWIDTH+(ny-1)*dy)*(jj*(_YWIDTH+dy)+_YWIDTH+dy/2.))
                thisphi= (phirange[0]+(phirange[1]-phirange[0])/
                          (nx*_XWIDTH+(nx-1)*dx)*(ii*(_XWIDTH+dx)+_XWIDTH+dx/2.))
            elif (ii >= nx/2-1 and ii < nx/2+nsx-1 and 
                  jj >= ny/2-1 and jj < ny/2 +nsy-1):
                continue
            else:
                thisax= fig.add_axes([(left+ii*(_XWIDTH+dx))/width,
                                      (bottom+jj*(_YWIDTH+dy))/height,
                                      _XWIDTH/width,_YWIDTH/height])
                thisR= (rrange[0]+(rrange[1]-rrange[0])/
                        (ny*_YWIDTH+(ny-1)*dy)*(jj*(_YWIDTH+dy)+_YWIDTH/2.))
                thisphi= (phirange[0]+(phirange[1]-phirange[0])/
                          (nx*_XWIDTH+(nx-1)*dx)*(ii*(_XWIDTH+dx)+_XWIDTH/2.))
            thissavefilename= os.path.join(saveDir,picklebasename+'_%i_%i.sav' %(ii,jj))
            if os.path.exists(thissavefilename):
                print "Restoring velocity distribution at %.1f, %.1f ..." %(thisR,thisphi)
                savefile= open(thissavefilename,'r')
                thisveldist= pickle.load(savefile)
                savefile.close()
            else:
                print "Calculating velocity distribution at %.1f, %.1f ..." %(thisR,thisphi)
                thisveldist= calc_veldist_2d(ulinspace,
                                             vlinspace,
                                             R=thisR,
                                             potparams=(0.9,0.01,
                                                        25.*_degtorad-thisphi,
                                                        .8,None))
                savefile= open(thissavefilename,'w')
                pickle.dump(thisveldist,savefile)
                savefile.close()
            fig.sca(thisax)
            plot.bovy_dens2d(thisveldist.T,origin='lower',
                             cmap='gist_yarg',
                             xrange=[ulinspace[0],ulinspace[1]],
                             yrange=[vlinspace[0],vlinspace[1]],
                             contours=True,cntrmass=True,
                             levels=levels,cntrcolors=cntrcolors,
                             overplot=True)
            thisax.xaxis.set_ticklabels('')
            thisax.yaxis.set_ticklabels('')
    plot.bovy_end_print(plotfilename)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Must provide a filename for the figure"
        sys.exit(-1)
    veldist_2d_Rphi(sys.argv[1])
