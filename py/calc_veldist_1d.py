import math as m
import scipy as sc
import scipy.integrate as integrate
from integrate_orbits import uvToELz
from interpret_as_df import dehnenDF
_degtorad= sc.pi/180.
_NCORRECT=20
_MAXITER= 50
_NSIGMA=5.
def marginalizeAngleGrid(vloslinspace,alpha,R=1.,t=-4.,pot='bar',beta=0.,
                         potparams=(0.9,0.01,20.*_degtorad,.8,None),
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
        print ii, out[ii]
    return out

def marginalizeAngle(vlos,alpha,R=1.,t=-4.,pot='bar',beta=0.,
                    potparams=(0.9,0.01,20.*_degtorad,.8,None),
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
    if dftype == 'dehnen':
        thisDF= dehnenDF
    df= thisDF(profileParams=dfparams,beta=beta,correct=correct,
               niter=_NCORRECT)
    intLimit= _NSIGMA*m.sqrt(df.targetSigma2(R))
    if m.fabs(m.sin(alpha)) < m.sqrt(1./2.):
        cosalpha= m.cos(alpha)
        tanalpha= m.tan(alpha)
        return integrate.quadrature(integrandSinAlphaSmall,-intLimit,intLimit,
                                    args=(cosalpha,tanalpha,vlos,R,t,potparams,
                                          beta,df,pot),
                                    vec_func=False,maxiter=_MAXITER)[0].real/m.fabs(cosalpha)
    else:
        sinalpha= m.sin(alpha)
        cotalpha= 1./m.tan(alpha)
        return integrate.quadrature(integrandSinAlphaLarge,-intLimit,intLimit,
                                    args=(sinalpha,cotalpha,vlos,R,t,potparams,
                                          beta,df,pot),
                                    vec_func=False,maxiter=_MAXITER)[0].real/m.fabs(sinalpha)

def integrandSinAlphaLarge(u,sinalpha,cotalpha,vlos,R,t,potparams,beta,df,pot):
    return df.eval(*uvToELz(UV=(-u,-cotalpha*u+vlos/sinalpha),R=R,t=t,pot=pot,potparams=potparams,beta=beta))

def  integrandSinAlphaSmall(u,cosalpha,tanalpha,vlos,R,t,potparams,beta,df,
                            pot):
    return df.eval(*uvToELz(UV=(tanalpha*u-vlos/cosalpha,u),R=R,t=t,pot=pot,potparams=potparams,beta=beta))

