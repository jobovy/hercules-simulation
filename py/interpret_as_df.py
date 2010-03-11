###############################################################################
#   interpret_as_df.py: module that interprets (E,Lz) pairs in terms of a 
#                        distribution function
###############################################################################
import os, os.path
import cPickle as pickle
import scipy as sc
import scipy.integrate as integrate
from integrate_orbits import vRvTRToEL
class distF:
    """Class that represents a DF"""
    def __init__(self,dftype='dehnen',dfparams=(0.33,1.0,0.2),beta=0.,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           Initialize a DF
        INPUT:
           dftype= 'dehnen' or 'corrected-dehnen'
           dfparams - parameters of the df (xD,xS,Sro)
                      xD - disk surface mass scalelength / Ro
                      xS - disk velocity dispersion scalelength / Ro
                      Sro - disk velocity dispersion at Ro (/vo)
           beta - power-law index of the rotation curve
           + DFcorrection kwargs
        OUTPUT:
        HISTORY:
            2010-03-10 - Written - Bovy (NYU)
        """
        self._dftype= dftype
        self._dfparams= dfparams
        self._beta= beta
        if dftype == 'corrected-dehnen':
            #Load corrections
            corr= DFcorrection(**kwargs)
        return None

    def eval(self,E,L):
        """
        NAME:
           eval
        PURPOSE:
           evaluate the distribution function
        INPUT:
           E - energy (/vo^2)
           L - angular momentun (/ro/vo)
        OUTPUT:
           DF(E,L)
        HISTORY:
           2010-03-10 - Written - Bovy (NYU)
        """
        if self._dftype == 'dehnen':
            return self._eval_dehnen(E,L)
        elif self._dftype == 'corrected-dehnen':
            return self._eval_corrected_dehnen(E,L)
    
    def _eval_dehnen(self,E,L):
        """Internal function that evaluates the uncorrected Dehnen DF"""
        #Calculate Re,LE, OmegaE
        if self._beta == 0.:
            xE= sc.exp(E-.5)
            LE= xE
            OmegaE= 1./xE
        else: #non-flat rotation curve
            xE= (2.*E/(1.+1./self._beta))**(1./2./self._beta)
            LE= xE**(self._beta+1.)
            OmegaE= xE**(self._beta-1.)
        SRE2= self._dfparams[2]**2.*sc.exp(-2.*(xE-1.)/self._dfparams[1])
        return 1./2./sc.pi*sc.sqrt(2./(1.+self._beta))*sc.exp(-xE/self._dfparams[0])/SRE2*sc.exp(OmegaE*(L-LE)/SRE2)
        
    def _eval_corrected_dehnen(self,E,L):
        """Internal function that evaluates the corrected Dehnen DF"""
        #Calculate Re,LE, OmegaE
        if self._beta == 0.:
            xE= sc.exp(E-.5)
            LE= xE
            OmegaE= 1./xE
        else: #non-flat rotation curve
            xE= (2.*E/(1.+1./self._beta))**(1./2./self._beta)
            LE= xE**(self._beta+1.)
            OmegaE= xE*(self._beta-1.)
        correction= corr.correct(xE)
        SRE2= self._dfparams[2]**2.*sc.exp(-2.*(xE-1.)/self._dfparams[1])*correction[1]
        return 1./2./sc.pi*sc.sqrt(2./(1.+self._beta))*sc.exp(-xE/self._dfparams[0])/SRE2*sc.exp(OmegaE*(L-LE)/SRE2)*correction[0]
        
    def _calc_surfacemass(self,R):
        """Internal function that calculates the surface mass for a given DF at R"""
        bound= 1.
        #BOVY: Normalize by (0.,0.) value first?
        return integrate.dblquad(_surfaceIntegrand,-bound,bound,
                                 lambda x: -bound, lambda x: bound,
                                 (R,self))

def _surfaceIntegrand(vR,vT,R,df):
    """Internal function that is the integrand for the surface mass integration"""
    E,L= vRvTRToEL(vR,vT,R,df._beta)
    return df.eval(E,L)



class DFcorrection:
    """Class that contains the corrections necessary to reach
    exponential profiles"""
    def __init__(self,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize the corrections: set them, load them, or calculate
           and save them
        OPTIONAL INPUTS:
           corrections - if Set, these are the corrections and they should
                         be used as such
           npoints - number of points from 0 to Rmax
           rmax - correct up to this radius (/ro)
           savedir - save the corrections in this directory
           dfparams - parameters of the df: xD, xS, Sro
           beta - power-law index of the rotation curve (when calculating)
           dftype - 'dehnen'
           niter - number of iterations to perform to calculate the corrections
        OUTPUT:
        HISTORY:
           2010-03-10 - Written - Bovy (NYU)
        """
        if not kwargs.has_key('dfparams'):
            self._dfparams= (0.33,1.,0.2)
        else:
            self._dfparams= kwargs['dfparams']
        if not kwargs.has_key('rmax'):
            self._rmax= 10.*self._dfparams[0]
        else:
            self._rmax= kwargs['rmax']
        if not kwargs.has_key('niter'):
            self._niter= 4
        else:
            self._niter= kwargs['niter']
        if not kwargs.has_key('npoints'):
            if kwargs.has_key('corrections'):
                self._npoints= kwargs['corrections'].shape[0]
            else:
                self._npoints= 101
        else:
            self._npoints= kwargs['npoints']
        if kwargs.has_key('dftype'):
            self._dftype= kwargs['dftype']
        else:
            self._dftype= 'dehnen'
        if kwargs.has_key('beta'):
            self._beta= kwargs['beta']
        else:
            self._beta= 0.
        if kwargs.has_key('corrections'):
            self._corrections= kwargs['corrections']
            if not len(self._corrections) == self._npoints:
                raise DFcorrectionError("Number of corrections has to be equal to the number of points npoints")
            return None
        if kwargs.has_key('savedir'):
            self._savedir= kwargs['savedir']
        else:
            self._savedir= '.'
        self._savefilename= os.path.join(self._savedir,'dfcorrection_'+
                                         self._dftype+
                                         '_%4.2f_%4.2f_%4.2f_%i_%4.2f_%i.sav'
                                         % (self._dfparams[0],self._dfparams[1],
                                            self._dfparams[2],self._npoints,
                                            self._rmax,self._niter))
        if os.path.exists(self._savefilename):
            self._corrections= pickle.load(self._savefilename)['corrections']
            return None
        #Calculate the corrections
        self._corrections= self._calc_corrections()
        return None

    def correct(self,R):
        """
        NAME:
           correct
        PURPOSE:
           calculate the correction in Sigma and sigma2 at R
        INPUT:
           R - Galactocentric radius(/ro)
        OUTPUT:
           [Sigma correction, sigma2 correction]
        HISTORY:
           2010-03-10 - Written - Bovy (NYU)
        """
        if R > self._rmax:
            return [1.,1.]
        nearR= round(R/self._rmax) #BOVY: INTERPOLATE?
        return self._corrections[nearR,:]

    def _calc_corrections(self):
        """Internal function that calculates the corrections"""
        self._uncorrdf= distF(dfparams=self._dfparams,dftype=self._dftype,beta=self._beta)
        corrections= sc.ones((self._npoints,2))
        rs= sc.linspace(0,self._rmax,self._npoints)
        for ii in range(self._niter):
            currentDF= distF(dftype='corrected-'+self.dftype,dfparams=self.dfparams,
                             beta=self._beta,corrections=corrections)
            this_surface= currentDF._calc_surfacemass(rs[ii])
            corrections[ii,0]*= this_surface/self._surfacemass(rs[ii])
            thisSigma2= currentDF._calc_sigma2surface(rs[ii])/this_surface
            corrections[ii,1]*= thisSigma2/self._sigma2(rs[ii])
        #Save
        saveDict= {'corrections': corrections}
        savefile= open(self._savefilename,'w')
        pickle.dump(saveDict,savefile)
        savefile.close()
    
    def _calc_sigma2surface(self,R):
        """Internal function that calculates sigma^2 Sigma for a given DF at R"""
        return None

    def _surfacemass(self,R):
        """Internal function that gives the surface mass at R"""
        return sc.exp(-R/self._dfparams[0])

    def _sigma2(self,R):
        """Internal function that gives the R-velocity variance at R"""
        return self_dfparams[2]**2.*sc.exp(-2.*(R-1.)/self._dfparams[1])

class DFcorrectionError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
        

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


