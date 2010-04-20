import os, os.path
import cPickle as pickle
import scipy as sc
from integrate_orbits import uvToELz
from interpret_as_df import dehnenDF
import bovy_plot as plot
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

if __name__ == '__main__':
    testN= 51
    #print calc_veldist_2d((-.5,.5,testN),(-.5,.5,testN))[testN/2,testN/2]
    print calc_veldist_2d((-.5,.5,testN),(-.5,.5,testN))
