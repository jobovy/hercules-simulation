import os, os.path
import sys
import cPickle as pickle
import math as m
import scipy as sc
from scipy import interpolate, stats
from optparse import OptionParser
import bovy_plot as plot
import probDistance
from matplotlib.pyplot import contour, clabel
_RADTODEG= 180./m.pi
_XWIDTH= 1.8*80/100/1.8
_YWIDTH= 1.15*80/20/1.8
def bar_detectability(parser,
                      dx=_XWIDTH/20.,dy=_YWIDTH/20.,
                      nx=100,ny=20,
                      ngrid=201,rrange=[0.7,1.3],
                      phirange=[-m.pi/2.,m.pi/2.],
                      saveDir='../bar/1dLarge/'):
    """
    NAME:
       bar_detectability
    PURPOSE:
       analyze the detectability of the Hercules moving group in the 
       los-distribution around the Galaxy
    INPUT:
       nx - number of plots in the x-direction
       ny - number of plots in the y direction
       dx - x-spacing
       dy - y-spacing
       ngrid - number of gridpoints to evaluate the density on
       rrange - range of Galactocentric radii to consider
       phirange - range of Galactic azimuths to consider
       saveDir - directory to save the pickles in
    OUTPUT:
       plot in plotfilename
    HISTORY:
       2010-05-09 - Written - Bovy (NYU)
    """
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return 
    
    if not options.convolve == None:
        bar_detectability_convolve(parser,dx=dx,dy=dy,nx=nx,ny=ny,ngrid=ngrid,
                                   rrange=rrange,phirange=phirange,
                                   saveDir=saveDir)
        return

    vloslinspace= (-.9,.9,ngrid)
    vloss= sc.linspace(*vloslinspace)

    picklebasename= '1d_%i_%i_%i_%.1f_%.1f_%.1f_%.1f' % (nx,ny,ngrid,rrange[0],rrange[1],phirange[0],phirange[1])

    detect= sc.zeros((nx,ny))
    losd= sc.zeros((nx,ny))
    gall= sc.zeros((nx,ny))
    for ii in range(nx):
        for jj in range(ny):
            thisR= (rrange[0]+(rrange[1]-rrange[0])/
                    (ny*_YWIDTH+(ny-1)*dy)*(jj*(_YWIDTH+dy)+_YWIDTH/2.))
            thisphi= (phirange[0]+(phirange[1]-phirange[0])/
                      (nx*_XWIDTH+(nx-1)*dx)*(ii*(_XWIDTH+dx)+_XWIDTH/2.))
            thissavefilename= os.path.join(saveDir,picklebasename+'_%i_%i.sav' %(ii,jj))
            if os.path.exists(thissavefilename):
                print "Restoring los-velocity distribution at %.2f, %.2f ..." %(thisR,thisphi)
                savefile= open(thissavefilename,'r')
                vlosd= pickle.load(savefile)
                axivlosd= pickle.load(savefile)
                savefile.close()
            else:
                print "Did not find the los-velocity distribution at at %.2f, %.2f ..." %(thisR,thisphi)
                print "returning ..."
                return
            ddx= 1./sc.sum(axivlosd)
            #skipCenter
            if not options.skipCenter == 0.:
                skipIndx= (sc.fabs(vloss) < options.skipCenter)
                indx= (sc.fabs(vloss) >= options.skipCenter)
                vlosd= vlosd/sc.sum(vlosd[indx])/ddx
                axivlosd= axivlosd/sc.sum(axivlosd[indx])/ddx
                vlosd[skipIndx]= 1.
                axivlosd[skipIndx]= 1.
            vlosd_zeroindx= (vlosd == 0.)
            axivlosd_zeroindx= (axivlosd == 0.)
            vlosd[vlosd_zeroindx]= 1.
            axivlosd[vlosd_zeroindx]= 1.
            vlosd[axivlosd_zeroindx]= 1.
            axivlosd[axivlosd_zeroindx]= 1.
            detect[ii,jj]= probDistance.kullbackLeibler(vlosd,axivlosd,ddx,nan=True)
            #los distance and Galactic longitude
            d= m.sqrt(thisR**2.+1.-2.*thisR*m.cos(thisphi))
            losd[ii,jj]= d
            if 1./m.cos(thisphi) < thisR and m.cos(thisphi) > 0.:
                l= m.pi-m.asin(thisR/d*m.sin(thisphi))
            else:
                l= m.asin(thisR/d*m.sin(thisphi))
            gall[ii,jj]= l

    #Find maximum, further than 3 kpc away
    detectformax= detect.flatten()
    detectformax[losd.flatten() < 3./8.2]= 0.
    x= sc.argmax(detectformax)
    indx = sc.unravel_index(x,detect.shape)
    maxR= (rrange[0]+(rrange[1]-rrange[0])/
           (ny*_YWIDTH+(ny-1)*dy)*(indx[1]*(_YWIDTH+dy)+_YWIDTH/2.))
    maxphi= (phirange[0]+(phirange[1]-phirange[0])/
                      (nx*_XWIDTH+(nx-1)*dx)*(indx[0]*(_XWIDTH+dx)+_XWIDTH/2.))
    print maxR, maxphi, losd[indx[0],indx[1]], detect[indx[0],indx[1]], gall[indx[0],indx[1]]*180./sc.pi

    #Now plot
    plot.bovy_print()
    plot.bovy_dens2d(detect.T,origin='lower',#interpolation='nearest',
                     xlabel=r'$\mathrm{Galactocentric\ azimuth}\ [\mathrm{deg}]$',
                     ylabel=r'$\mathrm{Galactocentric\ radius}\ /R_0$',
                     cmap='gist_yarg',xrange=sc.array(phirange)*_RADTODEG,
                     yrange=rrange,
                     aspect=(phirange[1]-phirange[0])*_RADTODEG/(rrange[1]-rrange[0]))
    #contour the los distance and gall
    #plot.bovy_text(-22.,1.1,r'$\mathrm{apogee}$',color='w',
    #                rotation=105.)
    plot.bovy_text(-18.,1.1,r'$\mathrm{APOGEE}$',color='w',
                    rotation=285.)
    levels= [2/8.2*(ii+1/2.) for ii in range(10)]
    contour(losd.T,levels,colors='0.25',origin='lower',linestyles='--',
            aspect=(phirange[1]-phirange[0])*_RADTODEG/(rrange[1]-rrange[0]),
            extent=(phirange[0]*_RADTODEG,phirange[1]*_RADTODEG,
                    rrange[0],rrange[1]))
    gall[gall < 0.]+= sc.pi*2.
    levels= [0.,sc.pi/2.,sc.pi,3.*sc.pi/2.]
    contour(gall.T,levels,colors='w',origin='lower',linestyles='--',
            aspect=(phirange[1]-phirange[0])*_RADTODEG/(rrange[1]-rrange[0]),
            extent=(phirange[0]*_RADTODEG,phirange[1]*_RADTODEG,
                    rrange[0],rrange[1]))
    levels= [-5/180.*sc.pi,250/180.*sc.pi]
    contour(gall.T,levels,colors='w',origin='lower',linestyles='-.',
            aspect=(phirange[1]-phirange[0])*_RADTODEG/(rrange[1]-rrange[0]),
            extent=(phirange[0]*_RADTODEG,phirange[1]*_RADTODEG,
                    rrange[0],rrange[1]))
    if options.skipCenter == 0.:
        plot.bovy_text(r'$\mathrm{KL\ divergence\ / \ all}\ v_{\mathrm{los}}$',
                       title=True)
    else:
        plot.bovy_text(r'$\mathrm{KL\ divergence\ / }\ |v_{\mathrm{los}}| \geq %.2f \ v_0$' % options.skipCenter,
                       title=True)
    plot.bovy_end_print(args[0])

def bar_detectability_convolve(parser,nconvsamples=1000,
                               dx=_XWIDTH/20.,dy=_YWIDTH/20.,
                               nx=100,ny=20,
                               ngrid=201,rrange=[0.7,1.3],
                               phirange=[-m.pi/2.,m.pi/2.],
                               saveDir='../bar/1dLarge/',
                               saveDirConv='../bar/1dLargeConv/'):
    """
    NAME:
       bar_detectability_convolve
    PURPOSE:
       analyze the detectability of the Hercules moving group in the 
       los-distribution around the Galaxy, convolving with distance 
       uncertainties
    INPUT:
       parser - from optparse
       nconvsamples - number of samples to take to perform the convolution
       nx - number of plots in the x-direction
       ny - number of plots in the y direction
       dx - x-spacing
       dy - y-spacing
       ngrid - number of gridpoints to evaluate the density on
       rrange - range of Galactocentric radii to consider
       phirange - range of Galactic azimuths to consider
       saveDir - directory to save the pickles in
    OUTPUT:
       plot in plotfilename
    HISTORY:
       2010-05-09 - Written - Bovy (NYU)
    """
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        return 

    vloslinspace= (-.9,.9,ngrid)
    vloss= sc.linspace(*vloslinspace)

    picklebasename= '1d_%i_%i_%i_%.1f_%.1f_%.1f_%.1f' % (nx,ny,ngrid,rrange[0],rrange[1],phirange[0],phirange[1])

    #First load all of the precalculated velocity distributions
    
    phis= sc.zeros(nx)
    rs= sc.zeros(ny)
    bisplphis= sc.zeros(nx*ny)
    bisplrs= sc.zeros(nx*ny)
    vlosds= sc.zeros((nx*ny,ngrid))
    axivlosds= sc.zeros((nx*ny,ngrid))
    for ii in range(nx):
        for jj in range(ny):
            thisR= (rrange[0]+(rrange[1]-rrange[0])/
                    (ny*_YWIDTH+(ny-1)*dy)*(jj*(_YWIDTH+dy)+_YWIDTH/2.))
            thisphi= (phirange[0]+(phirange[1]-phirange[0])/
                      (nx*_XWIDTH+(nx-1)*dx)*(ii*(_XWIDTH+dx)+_XWIDTH/2.))
            phis[ii]= thisphi
            rs[jj]= thisR
            bisplphis[ii*ny+jj]= thisphi
            bisplrs[ii*ny+jj]= thisR
            thissavefilename= os.path.join(saveDir,picklebasename+'_%i_%i.sav' %(ii,jj))
            if os.path.exists(thissavefilename):
                print "Restoring los-velocity distribution at %.2f, %.2f ..." %(thisR,thisphi)
                savefile= open(thissavefilename,'r')
                vlosd= pickle.load(savefile)
                axivlosd= pickle.load(savefile)
                savefile.close()
                vlosds[ii*ny+jj,:]= vlosd
                axivlosds[ii*ny+jj,:]= axivlosd
            else:
                print "Did not find the los-velocity distribution at at %.2f, %.2f ..." %(thisR,thisphi)
                print "returning ..."
                return
    
    #Now convolve and calculate Kullback-Leibler divergence
    picklebasename= '1d_%i_%i_%i_%.1f_%.1f_%.1f_%.1f_%.1f' % (nx,ny,ngrid,rrange[0],rrange[1],phirange[0],phirange[1],options.convolve)
    detect= sc.zeros((nx,ny))
    losd= sc.zeros((nx,ny))
    for ii in range(nx):
        for jj in range(ny):
            if ii == 45 and jj == 13:
                continue#BOVY: FIX FOR NOW
            thisR= rs[jj]
            thisphi= phis[ii]
            thissavefilename= os.path.join(saveDirConv,picklebasename+'_%i_%i.sav' %(ii,jj))
            if os.path.exists(thissavefilename):
                print "Restoring convolved los-velocity distribution at %.2f, %.2f ..." %(thisR,thisphi)
                savefile= open(thissavefilename,'r')
                convvlosd= pickle.load(savefile)
                convaxivlosd= pickle.load(savefile)
                savefile.close()
            else:
                print "Calculating convolved los-velocity distribution at %.2f, %.2f ..." %(thisR,thisphi)
                #los distance
                losd[ii,jj]= m.sqrt(thisR**2.+1.-2.*thisR*m.cos(thisphi))
                thislosd= losd[ii,jj]
                #Galactic longitude
                if 1./m.cos(thisphi) < thisR and m.cos(thisphi) > 0.:
                    thisl= m.pi-m.asin(thisR/thislosd*m.sin(thisphi))
                else:
                    thisl= m.asin(thisR/thislosd*m.sin(thisphi))
                convvlosd= sc.zeros(ngrid)
                convaxivlosd= sc.zeros(ngrid)
                broke= False
                for kk in range(ngrid):
                    weights= invdist2(bisplrs[ii*ny+jj],bisplphis[ii*ny+jj],
                                      bisplrs,bisplphis)
                    weights[sc.isinf(weights)]= sc.amax(weights[sc.isfinite(weights)])
                    splindx= (weights > 1./(thislosd*options.convolve*5.)**2.)
                    try:
                        tck= interpolate.bisplrep(bisplphis[splindx],
                                                  bisplrs[splindx],
                                                  vlosds[splindx,kk],
                                                  w=weights[splindx])
                    except TypeError:
                        #Trye interp2d?
                        broke= True
                        break
                    except ValueError:
                        splindx= (weights > 1./(2.*thislosd*options.convolve*5.)**2.)
                        try:
                            tck= interpolate.bisplrep(bisplphis[splindx],
                                                      bisplrs[splindx],
                                                      vlosds[splindx,kk],
                                                      w=weights[splindx])
                        except TypeError:
                            broke= True
                            break                      
                        except ValueError:
                            broke= True
                            break
                    #Now convolve
                    thisnsamples= 0
                    for ll in range(nconvsamples):
                        samplelosd= thislosd+stats.norm.rvs()*thislosd*options.convolve
                        sampleR, samplephi= dlToRphi(samplelosd,thisl)
                        addvlos= interpolate.bisplev(samplephi,sampleR,tck)
                        convvlosd[kk]+= addvlos
                        if not addvlos == 0.:
                            thisnsamples+= 1
                    convvlosd[kk]/= thisnsamples
                    try:
                        tck= interpolate.bisplrep(bisplphis[splindx],
                                                  bisplrs[splindx],
                                                  axivlosds[splindx,kk],
                                                  w=weights[splindx])
                    except TypeError:
                        broke= True
                        break
                    except ValueError:
                        splindx= (weights > 1./(2.*thislosd*options.convolve*5.)**2.)
                        try:
                            tck= interpolate.bisplrep(bisplphis[splindx],
                                                      bisplrs[splindx],
                                                      axivlosds[splindx,kk],
                                                      w=weights[splindx])
                        except TypeError:
                            broke= True
                            break                      
                        except ValueError:
                            broke= True
                            break
                    #Now convolve
                    thisnsamples= 0
                    for ll in range(nconvsamples):
                        samplelosd= thislosd+stats.norm.rvs()*thislosd*options.convolve
                        sampleR, samplephi= dlToRphi(samplelosd,thisl)
                        addvlos= interpolate.bisplev(samplephi,sampleR,tck)
                        convaxivlosd[kk]+= addvlos
                        if not addvlos == 0.:
                            thisnsamples+= 1
                    convaxivlosd[kk]/= thisnsamples
                savefile= open(thissavefilename,'w')
                pickle.dump(convvlosd,savefile)
                pickle.dump(convaxivlosd,savefile)
                savefile.close()
                if broke:
                    continue#BOVY: FIX FOR NOW
            ddx= 1./sc.sum(axivlosd)
            #skipCenter
            if not options.skipCenter == 0.:
                skipIndx= (sc.fabs(vloss) < options.skipCenter)
                indx= (sc.fabs(vloss) >= options.skipCenter)
                convvlosd= convvlosd/sc.sum(convvlosd[indx])/ddx
                convaxivlosd= convaxivlosd/sc.sum(convaxivlosd[indx])/ddx
                convvlosd[skipIndx]= 1.
                convaxivlosd[skipIndx]= 1.
            convvlosd_zeroindx= (convvlosd == 0.)
            convaxivlosd_zeroindx= (convaxivlosd == 0.)
            convvlosd[convvlosd_zeroindx]= 1.
            convaxivlosd[convvlosd_zeroindx]= 1.
            convvlosd[convaxivlosd_zeroindx]= 1.
            convaxivlosd[convaxivlosd_zeroindx]= 1.
            detect[ii,jj]= probDistance.kullbackLeibler(convvlosd,convaxivlosd,
                                                        ddx,nan=True)
    detect[(detect < 0)]= 0.
    detect[(detect > 0.07)]= 0.
    #Now plot
    plot.bovy_print()
    plot.bovy_dens2d(detect.T,origin='lower',#interpolation='nearest',
                     xlabel=r'$\mathrm{Galactocentric\ azimuth}\ [\mathrm{deg}]$',
                     ylabel=r'$\mathrm{Galactocentric\ radius}\ /R_0$',
                     cmap='gist_yarg',xrange=sc.array(phirange)*_RADTODEG,
                     yrange=rrange,
                     aspect=(phirange[1]-phirange[0])*_RADTODEG/(rrange[1]-rrange[0]))
    #contour the los distance
    levels= [2/8.2*(ii+1/2.) for ii in range(10)]
    contour(losd.T,levels,colors='0.25',origin='lower',linestyles='--',
            aspect=(phirange[1]-phirange[0])*_RADTODEG/(rrange[1]-rrange[0]),
            extent=(phirange[0]*_RADTODEG,phirange[1]*_RADTODEG,
                    rrange[0],rrange[1]))
    plot.bovy_end_print(args[0])

def dlToRphi(d,l):
    """Convert los distance and Galactic longitude into Galactocentric 
    radius and azimuth"""
    R= m.sqrt(1.+d**2.-2.*d*m.cos(l))
    asinarg= d/R*m.sin(l)
    if asinarg < -1.:
        asinarg= -1.
    elif asinarg > 1.:
        asinarg= 1.
    if 1./m.cos(l) < d and m.cos(l) > 0.:
        theta= m.pi-m.asin(asinarg)
    else:
        theta= m.asin(asinarg)
    return (R,theta)

def invdist2(r1,phi1,r2,phi2):
    """return the inverse square distance between two points in 
    Galactocentric coordinates"""
    dist2= r1**2.+r2**2.-2.*r1*r2*sc.cos(phi1-phi2)
    return 1./dist2

def get_options():
    usage = "usage: %prog [options] <plotfilename>\n\nplotfilename= name of the file that the figure will be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("--skipCenter", dest="skipCenter",type='float',
                      default=0.,
                      help="skip the central part of the velocity distribution, expressed in terms of vo")
    parser.add_option("--convolve", dest="convolve",type='float',
                      default=None,
                      help="Convolve with relative distance uncertainties of this magnitude")
    
    return parser

if __name__ == '__main__':
    bar_detectability(get_options())
