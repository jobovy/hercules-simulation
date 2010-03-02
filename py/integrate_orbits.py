###############################################################################
#   integrate_orbits.py: module that integrates orbits in a given 
#                        non-axisymmetric, time-dependent potential
###############################################################################
import scipy as sc
import scipy.integrate as integrate
_degtorad= sc.pi/180.
def uvToELz(UV=(0.,0.),Rtheta=(1.,0.),t=-4.,pot='bar',
            potparams=(0.9,0.01,25.*_degtorad,.8,0.,2.)):
    """
    NAME:
       uvToELz
    PURPOSE:
       calculate initial (E.Lz) for final (u,v)
    INPUT:
       (u,v) - final radial and tangential velocity, divided by vcirc
       (R,theta) - final R and theta
       t - time to integrate backwards for 
           (interpretation depends on potential)
       pot - type of non-axisymmetric, time-depedent potential
       potparams - parameters for this potential
    OUTPUT:
       final (E,Lz)
    HISTORY:
       2010-03-01 - Written - Bovy (NYU)
    """
    #For now assume pot == 'bar'
    (Rolr,alpha,phi,chi,beta,t1) = potparams
    OmegaoOmegab= Rolr**(1.-beta)/(1.+sc.sqrt((1.+beta)/2.))
    #Final conditions
    u,v= UV
    v+= 1. #Add circular velocity
    R,theta = Rtheta
    x= R*sc.cos(theta)
    y= R*sc.sin(theta)
    vR= (u*x+v*y)/R * OmegaoOmegab
    vT= (v*x-u*y)/R * OmegaoOmegab
    (vR,vT,R) = integrate_orbit((vR,vT,R),t=t,pot=pot,potparams=potparams)
    return (axipotential(R,beta)+0.5*vR**2.+0.5*vT**2.,vT*R) #Assumes no perturbation at this time

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
        return sc.log(R)
    else: #non-flat rotation curve
        return R**(2.*beta)/2./beta

def integrate_orbit(vRvTR= (0.,1.,1.),t=-4.,pot='bar',
                    potparams=(0.9,0.01,25.*_degtorad,.8,0.,2.)):
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
       potparams - parameters for this potential
    OUTPUT:
       (vRF,vTF,RF)
    HISTORY:
       2010-03-01 - Written - Bovy (NYU)
    """
    #For now assume pot == 'bar'
    (Rolr,alpha,phi,chi,beta,t1) = potparams
    OmegaoOmegab= Rolr**(1.-beta)/(1.+sc.sqrt((1.+beta)/2.))
    vR, vT, R= vRvTR
    h= R*vT #specific angular momentum
    #Integrate the orbit
    Rb= chi*OmegaoOmegab**(1./(1.-beta))
    args= (h,OmegaoOmegab,beta,alpha,phi,Rb,t1)
    time= [0.,t]
    initCond= [R,vR]
    intOut= integrate.odeint(barEOM,initCond,time,args=args,rtol=10.**-15.)
    vRF= intOut[1,1]
    RF= intOut[1,0]
    vTF= h/RF
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
    h, OmegaoOmegab,beta,alpha,phi,Rb, t1= args
    x= y[0]
    xi= 2.*t/t1-1.
    if x <= Rb:
        barsign= 1.
    else: #outside of bar
        barsign= -1.
    return [y[1],h**2./x**3.-OmegaoOmegab**2.*x**(2.*beta-1.)+barsign*alpha*OmegaoOmegab**2.*sc.cos(2.*(phi-t))*(3./16.*xi**5.-5./8*xi**3.+15./16.*xi+.5)]
            

if __name__ == '__main__':
    print uvToELz()
