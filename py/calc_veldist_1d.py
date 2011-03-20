import os, os.path
import cPickle as pickle
import math as m
import numpy as nu
import scipy as sc
import scipy.integrate as integrate
from integrate_orbits import uvToELz
from interpret_as_df import dehnenDF, shuDF
import galpy.util.bovy_plot as plot
_DEBUG=False
_degtorad= sc.pi/180.
_NCORRECT=20
_MAXITER= 20
_EPSREL= 1.e-02
_EPSABS= 1.e-05
_NSIGMA=5.
def plotVlos(vloslinspace,l=0.,d=1.,t=-4.,distCoord='GC',
             pot='bar',beta=0.,potparams=(0.9,0.01,25.*_degtorad,.8,None),
             dfparams=(1./3.,1.,0.2),dftype='dehnen',correct=True,
             plotfilename='../bar/vlosdist.ps',
             savefilename='../bar/vlosdist.sav',
             normalize=True,overplotAxiDF=True,labelPos=True):
    """
    NAME:
       plotVlos
    PURPOSE:
       predict and plot the line-of-sight velocity distribution in a given 
       direction and at a given distance
    INPUT:
       plotfilename - filename for plot
       savefilename - filename to save the vlos distribution
       normalize - if True, normalize to integrate to 1
       overplotAxiDF - if True, overplot the distribution of the axisymmetric 
                        DF
       labelPos - label the position of this direction in the Galaxy
       vloslinspace - los velocity grid to get the marginalized probability at
       l - Galactic longitude (rad)
       d - distance (from Sun or from GC)
       distCoord - 'GC' or 'Sun' (origin of d)
       t - time to integrate backwards for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-dependent potential
       beta - power-law index of rotation curve
       potparams - parameters for this potential
       dfparams - parameters of the DF (xD,xS,Sro)
       dftype - type of DF ('dehnen' or 'shu')
       correct - If True, correct the DF
    OUTPUT:
       [nvlos] array with the predicted vlos-distribution
    HISTORY:
       2010-04-14 - Written - Bovy (NYU)
    """
    if os.path.exists(savefilename):
        savefile= open(savefilename,'r')
        vlosd= pickle.load(savefile)
        savefile.close()
    else:
        vlosd= predictVlos(vloslinspace,l=l,d=d,t=t,distCoord=distCoord,
                           pot=pot,beta=beta,potparams=potparams,
                           dfparams=dfparams,dftype=dftype,correct=correct)
        savefile= open(savefilename,'w')
        pickle.dump(vlosd,savefile)
        savefile.close()
    vloss= sc.linspace(*vloslinspace)
    if normalize:
        vlosd= vlosd/(nu.sum(vlosd)*(vloss[1]-vloss[0]))
    if overplotAxiDF:
        Rolr,alpha,phi,chi,t1= potparams
        alpha= 0.
        potparams= (Rolr,alpha,phi,chi,t1)
        vlosdAxi= predictVlos(vloslinspace,l=l,d=d,t=0.0000001,
                              distCoord=distCoord,
                              pot=pot,beta=beta,potparams=potparams,
                              dfparams=dfparams,dftype=dftype,correct=correct)
        if normalize:
            vlosdAxi= vlosdAxi/(nu.sum(vlosdAxi)*(vloss[1]-vloss[0]))
    plot.bovy_print()
    plot.bovy_plot(vloss,vlosd,xlabel=r'$v_r/v_0$',
                   ylabel=r'$p(v_r)$',color='k',zorder=2)
    plot.bovy_plot(vloss,vlosdAxi,color='k',ls='--',zorder=1,overplot=True)
    if labelPos:
        if distCoord.lower() == 'gc':
            R= d
            if R < 1.:
                theta= m.asin(m.sin(l)/R)-l
            else:
                theta= m.pi-m.asin(m.sin(l)/R)-l
            d= m.sqrt(1.+R**2.-2.*R*m.cos(theta))
        llabel= nu.round(l/m.pi*180.)
        plot.bovy_text(r'$l=%i^\circ$' % llabel+'\n'+r'$d/R_0 = %.1f$' % d,
                       top_left=True)
    plot.bovy_end_print(plotfilename)

def predictVlosConvolve(vloslinspace,l=0.,d=1.,t=-4.,distCoord='GC',
                pot='bar',beta=0.,potparams=(0.9,0.01,25.*_degtorad,.8,None),
                dfparams=(1./3.,1.,0.2),dftype='dehnen',correct=True,
                convolve=0.1,nconvolve=101,sigmaconvolve=3):
    """
    NAME:
       predictVlosConvolve
    PURPOSE:
       predict the line-of-sight velocity distribution in a given direction
       and at a given distance, convolving with distance uncertainties
    INPUT:
       vloslinspace - los velocity grid to get the marginalized probability at
       l - Galactic longitude (rad)
       d - distance (from Sun or from GC)
       distCoord - 'GC', 'GCGC' or 'Sun' (origin of d, if GCGC l is Galactic 
                    too)
       t - time to integrate backwards for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-dependent potential
       beta - power-law index of rotation curve
       potparams - parameters for this potential
       dfparams - parameters of the DF (xD,xS,Sro)
       dftype - type of DF ('dehnen' or 'shu')
       correct - If True, correct the DF
       convolve - convolve with distance errors of this relative size 
                  (fraction)
       nconvolve - number of convolution gridpoints to use
       sigmaconvolve - number of sigmas to go out in the convolution
    OUTPUT:
       [nvlos] array with the predicted vlos-distribution
    HISTORY:
       2010-05-13 - Written - Bovy (NYU)
    """
    if distCoord.lower() == 'sun':
        R= m.sqrt(1.+d**2.-2.*d*m.cos(l))
        if 1./m.cos(l) < d and m.cos(l) > 0.:
            theta= m.pi-m.asin(d/R*m.sin(l))
        else:
            theta= m.asin(d/R*m.sin(l))
    elif distCoord.lower() == 'gcgc':
        R= d
        theta= l
        d= m.sqrt(R**2.+1.-2.*R*m.cos(theta))
        if 1./m.cos(theta) < R and m.cos(theta) > 0.:
            l= m.pi-m.asin(R/d*m.sin(theta))
        else:
            l= m.asin(R/d*m.sin(theta))
    else:
        R= d
        if R < 1.:
            if m.cos(l) < 0.:
                raise ValueError("Cannot probe R < 1. for 90. < l < 270.")
            theta= m.asin(m.sin(l)/R)-l
        else:
            theta= m.pi-m.asin(m.sin(l)/R)-l
        d= m.sqrt(R**2.+1.-2.*R*m.cos(theta))
    #d is now definitely the distance from the Sun, l is Galactic longitude
    dgrid= sc.linspace(-sigmaconvolve,sigmaconvolve,nconvolve)
    dgrid*= d*convolve
    dgrid+= d
    nvlos= vloslinspace[2]
    vlosd= sc.zeros(nvlos)
    norm= 0.
    for ii in range(nconvolve):
        vlosd+= (predictVlos(vloslinspace,l=l,d=dgrid[ii],t=t,distCoord='sun',
                            pot=pot,beta=beta,potparams=potparams,
                            dfparams=dfparams,dftype=dftype,correct=correct)
                 *m.exp(-(d-dgrid[ii])**2./2./d**2./convolve**2.))
        norm+= m.exp(-(d-dgrid[ii])**2./2./d**2./convolve**2.)
    return vlosd/norm

def predictVlos(vloslinspace,l=0.,d=1.,t=-4.,distCoord='GC',
                pot='bar',beta=0.,potparams=(0.9,0.01,25.*_degtorad,.8,None),
                dfparams=(1./3.,1.,0.2),dftype='dehnen',correct=True,
                convolve=None,nconvolve=101,sigmaconvolve=3):
    """
    NAME:
       predictVlos
    PURPOSE:
       predict the line-of-sight velocity distribution in a given direction
       and at a given distance
    INPUT:
       vloslinspace - los velocity grid to get the marginalized probability at
       l - Galactic longitude (rad)
       d - distance (from Sun or from GC)
       distCoord - 'GC', 'GCGC' or 'Sun' (origin of d, if GCGC l is Galactic 
                    too)
       t - time to integrate backwards for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-dependent potential
       beta - power-law index of rotation curve
       potparams - parameters for this potential
       dfparams - parameters of the DF (xD,xS,Sro)
       dftype - type of DF ('dehnen' or 'shu')
       correct - If True, correct the DF
    OUTPUT:
       [nvlos] array with the predicted vlos-distribution
    HISTORY:
       2010-04-14 - Written - Bovy (NYU)
    """
    if distCoord.lower() == 'sun':
        R= m.sqrt(1.+d**2.-2.*d*m.cos(l))
        if 1./m.cos(l) < d and m.cos(l) > 0.:
            theta= m.pi-m.asin(d/R*m.sin(l))
        else:
            theta= m.asin(d/R*m.sin(l))
    elif distCoord.lower() == 'gcgc':
        R= d
        theta= l
        d= m.sqrt(R**2.+1.-2.*R*m.cos(theta))
        if 1./m.cos(theta) < R and m.cos(theta) > 0.:
            l= m.pi-m.asin(R/d*m.sin(theta))
        else:
            l= m.asin(R/d*m.sin(theta))
    else:
        R= d
        if R < 1.:
            if m.cos(l) < 0.:
                raise ValueError("Cannot probe R < 1. for 90. < l < 270.")
            theta= m.asin(m.sin(l)/R)-l
        else:
            theta= m.pi-m.asin(m.sin(l)/R)-l
        d= m.sqrt(R**2.+1.-2.*R*m.cos(theta))
    #d is now definitely the distance from the Sun
    alphalos= theta+l
    Rolr,alpha,phi,chi,t1= potparams
    phi-=theta
    if _DEBUG:
        print "Theta: ", theta
        print "alpha: ", alphalos
        print "R: ", R
        print "d: ", d
    potparams= (Rolr,alpha,phi,chi,t1)
    return marginalizeAngleGrid(vloslinspace,alphalos,R=R,t=t,pot=pot,
                                beta=beta,potparams=potparams,
                                dfparams=dfparams,dftype=dftype,
                                correct=correct)

def marginalizeAngleGrid(vloslinspace,alpha,R=1.,t=-4.,pot='bar',beta=0.,
                         potparams=(0.9,0.01,25.*_degtorad,.8,None),
                         dfparams=(1./3.,1.,0.2),dftype='dehnen',
                         correct=True):
    """
    NAME:
       marginalizeAngleGrid
    PURPOSE:
       marginalize the 2d velocity distribution onto the direction making an
       angle alpha with v=0 of alpha, on a grid in vlos
    INPUT:
       vloslinspace - los velocity grid to get the marginalized probability at
       alpha - angle between los and v=0 (in radians)
       R - Galactocentric Radius
       t - time to integrate backwards for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-dependent potential
       beta - power-law index of rotation curve
       potparams - parameters for this potential
       dfparams - parameters of the DF (xD,xS,Sro)
       dftype - type of DF ('dehnen' or 'shu')
       correct - If True, correct the DF
    OUTPUT:
       value of marginalized DF
    HISTORY:
       2010-04-04 - Written - Bovy (NYU)
    """
    vloss= sc.linspace(*vloslinspace)
    out= sc.zeros(len(vloss))
    for ii in range(len(vloss)):
        out[ii]= marginalizeAngle(vloss[ii],alpha,R=R,t=t,pot=pot,beta=beta,
                                  potparams=potparams,
                                  dfparams=dfparams,dftype=dftype,
                                  correct=correct)
        if _DEBUG:
            print ii, out[ii]
    return out

def marginalizeAngle(vlos,alpha,R=1.,t=-4.,pot='bar',beta=0.,
                    potparams=(0.9,0.01,25.*_degtorad,.8,None),
                    dfparams=(1./3.,1.,0.2),dftype='dehnen',
                    correct=True):
    """
    NAME:
       marginalizeAngle
    PURPOSE:
       marginalize the 2d velocity distribution onto the direction making an
       angle alpha with v=0 of alpha
    INPUT:
       vlos - los velocity to get the marginalized probability at
       alpha - angle between los and v=0 (in radians)
       R - Galactocentric Radius
       t - time to integrate backwards for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-dependent potential
       beta - power-law index of rotation curve
       potparams - parameters for this potential
       dfparams - parameters of the DF (xD,xS,Sro)
       dftype - type of DF ('dehnen' or 'shu')
       correct - If True, correct the DF
    OUTPUT:
       value of marginalized DF
    HISTORY:
       2010-04-04 - Written - Bovy (NYU)
    """
    if dftype.lower() == 'dehnen':
        thisDF= dehnenDF
    elif dftype.lower() == 'shu':
        thisDF= shuDF
    df= thisDF(profileParams=dfparams,beta=beta,correct=correct,
               niter=_NCORRECT)
    intLimit= _NSIGMA*m.sqrt(df.targetSigma2(R))
    if m.fabs(m.sin(alpha)) < m.sqrt(1./2.):
        cosalpha= m.cos(alpha)
        tanalpha= m.tan(alpha)
        #return integrate.quadrature(integrandSinAlphaSmall,-intLimit,intLimit,
        #                            args=(cosalpha,tanalpha,vlos,R,t,potparams,
        #                                  beta,df,pot),
        #                            vec_func=False,maxiter=_MAXITER)[0].real/m.fabs(cosalpha)
        return integrate.quad(integrandSinAlphaSmall,-intLimit,intLimit,
                              args=(cosalpha,tanalpha,vlos,R,t,potparams,
                                    beta,df,pot),epsrel=_EPSREL,epsabs=_EPSABS,
                              limit=_MAXITER)[0]/m.fabs(cosalpha)
    else:
        sinalpha= m.sin(alpha)
        cotalpha= 1./m.tan(alpha)
        #return integrate.quadrature(integrandSinAlphaLarge,-intLimit,intLimit,
        #                            args=(sinalpha,cotalpha,vlos,R,t,potparams,
        #                                  beta,df,pot),
        #                            vec_func=False,maxiter=_MAXITER)[0].real/m.fabs(sinalpha)
        return integrate.quad(integrandSinAlphaLarge,-intLimit,intLimit,
                              args=(sinalpha,cotalpha,vlos,R,t,potparams,
                                    beta,df,pot),epsrel=_EPSREL,epsabs=_EPSABS,
                              limit=_MAXITER)[0]/m.fabs(sinalpha)

def integrandSinAlphaLarge(u,sinalpha,cotalpha,vlos,R,t,potparams,beta,df,pot):
    return df.eval(*uvToELz(UV=(-u,-cotalpha*u+vlos/sinalpha),R=R,t=t,pot=pot,potparams=potparams,beta=beta))

def integrandSinAlphaSmall(u,cosalpha,tanalpha,vlos,R,t,potparams,beta,df,
                            pot):
    return df.eval(*uvToELz(UV=(tanalpha*u-vlos/cosalpha,u),R=R,t=t,pot=pot,potparams=potparams,beta=beta))

