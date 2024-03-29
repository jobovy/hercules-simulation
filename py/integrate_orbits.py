###############################################################################
#   integrate_orbits.py: module that integrates orbits in a given 
#                        non-axisymmetric, time-dependent potential
###############################################################################
import math as m
import scipy as sc
import scipy.integrate as integrate
_degtorad= sc.pi/180.
def uvToELz_grid(ulinspace,vlinspace,R=1.,t=-4.,pot='bar',beta=0.,
            potparams=(0.9,0.01,25.*_degtorad,.8,None)):
    """
    NAME:
       uvToELz_grid
    PURPOSE:
       calculate uvToLz on a grid in (u,v)
    INPUT:
       ulinspace, vlinspace - build the grid using scipy's linspace with
                              these arguments
       R - Galactocentric Radius
       t - time to integrate backwards for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-dependent potential
       beta - power-law index of rotation curve
       potparams - parameters for this potential
    OUTPUT:
       final (E,Lz) on grid [nus,nvs,2]
       E=E/vo^2; Lz= Lz/Ro/vo
    HISTORY:
       2010-03-01 - Written - Bovy (NYU)
    """
    us= sc.linspace(*ulinspace)
    vs= sc.linspace(*vlinspace)
    nus= len(us)
    nvs= len(vs)
    out= sc.zeros((nus,nvs,2))
    for ii in range(nus):
        for jj in range(nvs):
            tmp_out= uvToELz(UV=(us[ii],vs[jj]),R=R,t=t,pot=pot,beta=beta,potparams=potparams)
            out[ii,jj,0]= tmp_out[0]
            out[ii,jj,1]= tmp_out[1]
    return out

def uvToELz(UV=(0.,0.),R=1.,t=-4.,pot='bar',beta=0.,
            potparams=(0.9,0.01,25.*_degtorad,.8,None)):
    """
    NAME:
       uvToELz
    PURPOSE:
       calculate initial (E.Lz) for final (u,v)
    INPUT:
       (u,v) - final radial and tangential velocity, divided by vcirc
               ACTUALLY -U; U towards GC!!
       R = Galactocentric radius
       t - time to integrate backwards for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-dependent potential
       beta - power-law index of rotation curve
       potparams - parameters for this potential
                   (Rolr,alpha,phi,chi,t1)
    OUTPUT:
       final (E,Lz)
       E=E/vo^2; Lz= Lz/Ro/vo
    HISTORY:
       2010-03-01 - Written - Bovy (NYU)
    """
    #For now assume pot == 'bar'
    (Rolr,alpha,phi,chi,t1) = potparams
    OmegaoOmegab= (Rolr**(1.-beta))/(1.+m.sqrt((1.+beta)/2.))
    #Final conditions
    u,v= UV
    v+= R**beta #Add circular velocity
    vR= u * OmegaoOmegab
    vT= v * OmegaoOmegab
    (vR,vT,R) = integrate_orbit((vR,vT,R),t=t,pot=pot,beta=beta,
                                potparams=potparams)
    vR/= OmegaoOmegab
    vT/= OmegaoOmegab
    return vRvTRToEL(vR,vT,R,beta)

def vRvTRToEL(vR,vT,R,beta):
    """
    NAME:
       vRvTRToEL
    PURPOSE:
       calculate the energy and angular momentum
    INPUT:
       vR - radial velocity
       vT - tangential velocity
       R - Galactocentric radius
    OUTPUT:
    HISTORY:
       2010-03-10 - Written - Bovy (NYU)
    """
    return (axipotential(R,beta)+0.5*vR**2.+0.5*vT**2.,vT*R)

def axipotential(R,beta=0.):
    """
    NAME:
       axipotential
    PURPOSE:
       return the axisymmetric potential at R/Ro
    INPUT:
       R - Galactocentric radius
       beta - rotation curve power-law
    OUTPUT:
       Pot(R)/vo**2.
    HISTORY:
       2010-03-01 - Written - Bovy (NYU)
    """
    if beta == 0.:
        return m.log(R)
    else: #non-flat rotation curve
        return R**(2.*beta)/2./beta

def integrate_orbit(vRvTR= (0.,1.,1.),t=-4.,pot='bar',beta=0.,
                    potparams=(0.9,0.01,25.*_degtorad,.8,None)):
    """
    NAME:
       integrate_orbit
    PURPOSE:
       integrate an orbit in a given potential
    INPUT:
       (vR,vT,R) - initial velocity in cylindrical coordinates +
                           initial R
       t - time to integrate for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-depedent potential
       beta - power-law index of rotation curve
       potparams - parameters for this potential
    OUTPUT:
       (vRF,vTF,RF)
    HISTORY:
       2010-03-01 - Written - Bovy (NYU)
    """
    #For now assume pot == 'bar'
    (Rolr,alpha,phi,chi,t1) = potparams
    OmegaoOmegab= (Rolr**(1.-beta))/(1.+m.sqrt((1.+beta)/2.))
    vR, vT, R= vRvTR
    vphi= vT/R #Angular velocity
    #Integrate the orbit
    Rb= chi*OmegaoOmegab**(1./(1.-beta)) #OmegabOmegao**-1. folded in
    if t1 == None:
        t1= m.fabs(t)/2.
    t1*= 2.*sc.pi
    tend= 2.*sc.pi*t
    args= (OmegaoOmegab,beta,alpha,Rb,t1,tend)
    time= [0.,tend]
    initCond= [R,vR,-phi,vphi]
    intOut= integrate.odeint(barEOM,initCond,time,args=args,
                             rtol=10.**-15.,mxstep=100000000)
    vRF= intOut[1,1]
    RF= intOut[1,0]
    vTF= RF*intOut[1,3]
    return (vRF,vTF,RF)

def barEOM(y,t,*args):
    """
    NAME:
       barEOM
    PURPOSE:
       function implementing the bar+rotation curve equation of motion
    INPUT:
       y - current position and velocity
       t - time t
       args - arguments = (h,OmegaoOmegab,beta,alpha,phi,Rb)
    OUTPUT:
       dy/dt
    HISTORY:
       2010-03-01 - Written - Bovy (NYU)
    """
    OmegaoOmegab,beta,alpha,Rb, t1, tend= args
    x= y[0]
    deltat= m.fabs(t-tend)
    if deltat < t1:
        xi= 2.*deltat/t1-1.
        smooth= (3./16.*xi**5.-5./8*xi**3.+15./16.*xi+.5)
    else: #bar is fully on
        smooth= 1.
    vx= y[1]
    phi= y[2]
    vphi= y[3]
    h= x**2.*vphi
    if x <= Rb:
        return [vx,h**2./x**3.-OmegaoOmegab**2.*x**(2.*beta-1.)
                -alpha*OmegaoOmegab**2.*m.cos(2.*(phi-t))*smooth*x**2./Rb**6.,
                vphi,
                2./3.*alpha*OmegaoOmegab**2.*smooth*m.sin(2.*(phi-t))*
                1./Rb**3.*((x/Rb)**3.-2)/x**2.-2.*vx*vphi/x]
    else: #outside of bar
        return [vx,h**2./x**3.-OmegaoOmegab**2.*x**(2.*beta-1.)
                -alpha*OmegaoOmegab**2.*m.cos(2.*(phi-t))*smooth/x**4.,
                vphi,
                -2./3.*alpha*OmegaoOmegab**2.*smooth*m.sin(2.*(phi-t))/x**5.
                -2.*vx*vphi/x]
            

if __name__ == '__main__':
    #Various tests
    print uvToELz((0,0))
    print uvToELz((-.1,-.1))

    """
    import sys
    sys.exit(-1)
    """

    import timeit
    
    s="""r=uvToELz((1.8*numpy.random.random()-.9,1.15*numpy.random.random()-.45))"""
    t = timeit.Timer(stmt=s,setup="from __main__ import uvToELz\nimport numpy")
    print "uvToELz: %.4f sec/pass" % (t.timeit(number=500)/500)

    if False:
        testN= 101
        print uvToELz_grid((-.5,.5,testN),(-.5,.5,testN))[testN/2,testN/2,:]
