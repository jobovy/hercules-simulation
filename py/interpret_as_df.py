###############################################################################
#   interpret_as_df.py: module that interprets (E,Lz) pairs in terms of a 
#                        distribution function
###############################################################################
import scipy as sc
def interpret_as_df(EL,beta=0.,xD=0.33,xS=1.,Sro=0.2,type='dehnen'):
    """
    NAME:
       interpret_as_df
    PURPOSE:
       interpret (E,L) as a DF (probability)
    INPUT:
       EL - (E,L) - Energy and angular momentum (E/vo^2,L/Ro/vo)
       beta - rotation curve power-law
       xD - disk surface mass scalelength / Ro
       xS - disk velocity dispersion scalelength / Ro
       Sro - disk velocity dispersion at Ro (/vo)
       type - type of DF ('dehnen' or 'shu')
    OUTPUT:
       DF(E,L)
    HISTORY:
       2010-03-08 - Written - Bovy (NYU)
    """
    if type == 'dehnen':
        return _interpret_as_df_dehnen(EL,beta,xD,xS,Sro)

def _interpret_as_df_dehnen(EL,beta,xD,xS,Sro):
    E, L= EL
    #Calculate Re,LE, OmegaE
    if beta == 0.:
        xE= sc.exp(E-.5)
        LE= xE
        OmegaE= 1./xE
    else: #non-flat rotation curve
        xE= (2.*E/(1.+1./beta))**(1./2./beta)
        LE= xE**(beta+1.)
        OmegaE= xE*(beta-1.)
    SRE2= Sro**2.*sc.exp(-2.*(xE-1.)/xS)
    return sc.exp(-xE/xD)/SRE2*sc.exp(OmegaE*(L-LE)/SRE2)

if __name__ == '__main__':
    print interpret_as_df((0.54636764432720319, 0.99999999999999989),type='dehnen')

    from integrate_orbits import uvToELz
    """
    step=0.1
    EL= uvToELz((0.,0.))
    print "step= ", step
    print "(0,0):",EL, interpret_as_df(EL)

    EL= uvToELz((step,0.))
    print "(1,0):",EL, interpret_as_df(EL)

    EL= uvToELz((0.,step))
    print "(0,1):",EL, interpret_as_df(EL)

    EL= uvToELz((-step,0.))
    print "(-1,0):",EL, interpret_as_df(EL)

    EL= uvToELz((0.,-step))
    print "(0,-1):",EL, interpret_as_df(EL)
    """

    _degtorad= sc.pi/180.
    print interpret_as_df(uvToELz((0,0),potparams=(1.,0.01,25.*_degtorad,.8,None)))
    print interpret_as_df(uvToELz((-.1,-.1),potparams=(1.,0.01,25.*_degtorad,.8,None)))


