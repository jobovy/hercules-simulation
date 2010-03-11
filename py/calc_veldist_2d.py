import scipy as sc
from integrate_orbits import uvToELz
from interpret_as_df import interpret_as_df, distF
_degtorad= sc.pi/180.
def calc_veldist_2d(ulinspace,vlinspace,R=1.,t=-4.,pot='bar',beta=0.,
                    potparams=(0.9,0.01,20.*_degtorad,.8,None),
                    dfparams=(0.33,1.,0.2),dftype='dehnen'):
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
    df= distF(dftype=dftype,dfparams=dfparams,beta=beta)
    for ii in range(nus):
        for jj in range(nvs):
            E,L= uvToELz(UV=(-us[ii],vs[jj]),R=R,t=t,pot=pot,potparams=potparams,beta=beta)
            out[ii,jj]= df.eval(E,L)
    return out

if __name__ == '__main__':
    testN= 51
    #print calc_veldist_2d((-.5,.5,testN),(-.5,.5,testN))[testN/2,testN/2]
    print calc_veldist_2d((-.5,.5,testN),(-.5,.5,testN))
