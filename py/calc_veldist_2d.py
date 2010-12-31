###############################################################################
#   calc_veldist_2d: 2d velocity distribution functions
#
#   Includes:
#      plot_veldist_2d
#      calc_veldist_2d
#      calc_meanR_2d: calculate the mean radii as a function of u,v
#      plot_meanR_2d
###############################################################################
import os, os.path
import cPickle as pickle
import scipy as sc
from integrate_orbits import uvToELz
from interpret_as_df import dehnenDF
import galpy.util.bovy_plot as plot
_degtorad= sc.pi/180.
_NCORRECT=20
def plot_veldist_2d(ulinspace=(-0.9,0.9,201),
                    vlinspace=(-.7,.45,201),
                    R=1.,t=-4.,pot='bar',beta=0.,
                    potparams=(0.9,0.01,20.*_degtorad,.8,None),
                    dfparams=(1./3.,1.,0.2),dftype='dehnen',
                    correct=True,plotfilename='../bar/veldist.ps',
                    savefilename='../bar/veldist.sav',
                    label=None):
    """
    NAME:
       plot_veldist_2d
    PURPOSE:
       plot the velocity distribution at a given R on a grid in (u,v)
    INPUT:
       ulinspace, vlinspace - build the grid using scipy's linspace with
                              these arguments
       R - Galactocentric Radius
       t - time to integrate backwards for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-dependent potential
       beta - power-law index of rotation curve
       potparams - parameters for this potential
       dfparams - parameters of the DF (xD,xS,Sro)
       dftype - type of DF ('dehnen' or 'shu')
       correct - If True, correct the DF
       savefilename - filename for savefile
       plotfilename - filename for plot
       label - a label to put in the upper-left corner of the plot
    OUTPUT:
       f(u,v) on the grid
    HISTORY:
       2010-03-29 - Written - Bovy (NYU)
    """
    if os.path.exists(savefilename):
        savefile= open(savefilename,'r')
        vd= pickle.load(savefile)
        savefile.close()
    else:
        vd= calc_veldist_2d(ulinspace,vlinspace,R=R,t=t,pot=pot,beta=beta,
                            potparams=potparams,dfparams=dfparams,
                            dftype=dftype,correct=correct)
        savefile= open(savefilename,'w')
        pickle.dump(vd,savefile)
        savefile.close()
    levels= sc.array([2,6,12,21,33,50,68,80,90,95,99,99.9])/100.
    cntrcolors= ['w' for ii in range(len(levels)) if levels[ii] <= .5]
    cntrcolors+= ['k' for ii in range(len(levels)) if levels[ii] > .5]
    plot.bovy_print()
    plot.bovy_dens2d(vd.T,origin='lower',cmap='gist_yarg',
                     xrange=[ulinspace[0],ulinspace[1]],
                     yrange=[vlinspace[0],vlinspace[1]],
                     contours=True,cntrmass=True,
                     levels=levels,cntrcolors=cntrcolors,
                     xlabel=r'$u/v_0$',ylabel='$v/v_0$')
    if not label == None:
        plot.bovy_text(label,top_left=True)
    plot.bovy_end_print(plotfilename)

def calc_veldist_2d(ulinspace,vlinspace,R=1.,t=-4.,pot='bar',beta=0.,
                    potparams=(0.9,0.01,25.*_degtorad,.8,None),
                    dfparams=(1./3.,1.,0.2),dftype='dehnen',
                    correct=True):
    """
    NAME:
       calc_veldist_2d
    PURPOSE:
       calculate the velocity distribution at a given R on a grid in (u,v)
    INPUT:
       ulinspace, vlinspace - build the grid using scipy's linspace with
                              these arguments
       R - Galactocentric Radius
       t - time to integrate backwards for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-dependent potential
       beta - power-law index of rotation curve
       potparams - parameters for this potential
       dfparams - parameters of the DF (xD,xS,Sro)
       dftype - type of DF ('dehnen' or 'shu')
    OUTPUT:
       f(u,v) on the grid
    HISTORY:
       2010-03-08 - Written - Bovy (NYU)
    """
    us= sc.linspace(*ulinspace)
    vs= sc.linspace(*vlinspace)
    nus= len(us)
    nvs= len(vs)
    out= sc.zeros((nus,nvs))
    if dftype == 'dehnen':
        thisDF= dehnenDF
    df= thisDF(profileParams=dfparams,beta=beta,correct=correct,
               niter=_NCORRECT)
    for ii in range(nus):
        for jj in range(nvs):
            E,L= uvToELz(UV=(-us[ii],vs[jj]),R=R,t=t,pot=pot,potparams=potparams,beta=beta)
            out[ii,jj]= df.eval(E,L)
    return out

def plot_meanR_2d(ulinspace,vlinspace,R=1.,t=-4.,pot='bar',beta=0.,
                  potparams=(0.9,0.01,25.*_degtorad,.8,None),
                  plotfilename='../talk-figures/meanR.ps',
                  savefilename='../talk-figures/meanR.sav',
                  label=None,
                  correct=True,xE=True):
    """
    NAME:
       plot_meanR_2d
    PURPOSE:
       plot the mean-radii distribution at a given R on a grid in (u,v)
    INPUT:
       ulinspace, vlinspace - build the grid using scipy's linspace with
                              these arguments
       R - Galactocentric Radius
       t - time to integrate backwards for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-dependent potential
       beta - power-law index of rotation curve
       potparams - parameters for this potential
       xE - if True, calculate mean radius based on energy, 
            else calculate mean radius based on angular momentum 
            (default: True)
       savefilename - filename for savefile
       plotfilename - filename for plot
       label - a label to put in the upper-left corner of the plot
    OUTPUT:
       Rmean(u,v) on the grid
    HISTORY:
       2010-05-03 - Written - Bovy (NYU)
    """
    if os.path.exists(savefilename):
        savefile= open(savefilename,'r')
        vd= pickle.load(savefile)
        savefile.close()
    else:
        vd= calc_meanR_2d(ulinspace,vlinspace,R=R,t=t,pot=pot,beta=beta,
                          potparams=potparams,correct=correct,xE=xE)
        savefile= open(savefilename,'w')
        pickle.dump(vd,savefile)
        savefile.close()
    vd[vd > 1.3]= 1.3
    levels= [0.6+ii*0.1 for ii in range(7)]
    cntrlabelcolors= ['w' for ii in range(4)]
    cntrlabelcolors+= ['k' for ii in range(3)]
    cntrcolors= cntrlabelcolors
    plot.bovy_print()
    plot.bovy_dens2d(vd.T,origin='lower',cmap='gist_gray',
                     xrange=[ulinspace[0],ulinspace[1]],
                     yrange=[vlinspace[0],vlinspace[1]],
                     contours=True,cntrmass=False,cntrlabel=True,
                     cntrinline=True,cntrcolors=cntrcolors,
                     levels=levels,cntrlabelcolors=cntrlabelcolors,
                     xlabel=r'$u/v_0$',ylabel='$v/v_0$')
    if not label == None:
        plot.bovy_text(label,top_left=True)
    plot.bovy_end_print(plotfilename)

def calc_meanR_2d(ulinspace,vlinspace,R=1.,t=-4.,pot='bar',beta=0.,
                  potparams=(0.9,0.01,25.*_degtorad,.8,None),
                  correct=True,xE=True):
    """
    NAME:
       calc_meanR_2d
    PURPOSE:
       calculate the mean-radii distribution at a given R on a grid in (u,v)
    INPUT:
       ulinspace, vlinspace - build the grid using scipy's linspace with
                              these arguments
       R - Galactocentric Radius
       t - time to integrate backwards for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-dependent potential
       beta - power-law index of rotation curve
       potparams - parameters for this potential
       xE - if True, calculate mean radius based on energy, 
            else calculate mean radius based on angular momentum 
            (default: True)
    OUTPUT:
       f(u,v) on the grid
    HISTORY:
       2010-05-03 - Written - Bovy (NYU)
    """
    us= sc.linspace(*ulinspace)
    vs= sc.linspace(*vlinspace)
    nus= len(us)
    nvs= len(vs)
    out= sc.zeros((nus,nvs))
    for ii in range(nus):
        for jj in range(nvs):
            E,L= uvToELz(UV=(-us[ii],vs[jj]),R=R,t=t,pot=pot,potparams=potparams,beta=beta)
            out[ii,jj]= calc_meanR(E,L,beta,xE=xE)
    return out

def calc_meanR(E,L,beta,xE=True):
    """
    NAME:
       calc_meanR
    PURPOSE:
       calculate the mean radius of an orbit based on its energy and angular 
       momentum
    INPUT:
       E - Enerhy
       L - angular momentum
       beta - exponent of the rotation curve
       xE - if True, calculate mean radius based on energy, 
            else calculate mean radius based on angular momentum 
            (default: True)
    OUTPUT:
       mean radius
    HISTORY:
       2010-05-03 - Written - Bovy (NYU)
    """
    if xE:
        if beta == 0.:
            xE= sc.exp(E-.5)
        else: #non-flat rotation curve
            xE= (2.*E/(1.+1./beta))**(1./2./beta)
    else:
        xE= L**(1./(beta+1.))
    return xE

if __name__ == '__main__':
    testN= 51
    #print calc_veldist_2d((-.5,.5,testN),(-.5,.5,testN))[testN/2,testN/2]
    print calc_veldist_2d((-.5,.5,testN),(-.5,.5,testN))
