import numpy as nu
from scipy import stats
from scipy import interpolate

def sample_vd(vloslinspace,vd,nsample=1):
    """
    NAME:
       sample_vd
    PURPOSE:
       sample points from a velocity distribution
    INPUT:
       vloslinspace - argument for linspace
       vd - velocity distribution over this range
       nsample - number of samples
    OUTPUT:
       list of samples
    HISTORY:
       2010-09-12 - Written - Bovy (NYU)
    """
    vloss= nu.linspace(*vloslinspace)
    vdmax= nu.amax(vd)
    vd= interpolate.InterpolatedUnivariateSpline(vloss,vd)
    out= []
    while len(out) < nsample:
        prop= stats.uniform.rvs()*(vloss[-1]-vloss[0])+vloss[0]
        if stats.uniform.rvs()*vdmax < vd(prop):
            out.append(prop)
    return out
