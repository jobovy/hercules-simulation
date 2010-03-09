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
        xE= sc.exp(2.*E-1.)
        LE= xE
        OmegaE= 1./xE
    else: #non-flat rotation curve
        xE= (2.*E/(1.+1./2./beta))**(1./2./beta)
        LE= xE**(beta+1.)
        xE= xE*(beta-1.)
    SRE2= Sro**2.*sc.exp(-2.*(xE-1.)/xS)
    return sc.exp(-xE/xD)/SRE2*sc.exp(OmegaE*(L-LE)/SRE2)

if __name__ == '__main__':
    print interpret_as_df((0.54636764432720319, 0.99999999999999989),type='dehnen')
