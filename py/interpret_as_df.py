###############################################################################
#   interpret_as_df.py: module that interprets (E,Lz) pairs in terms of a 
#                        distribution function
#
#   ToDo:
#      Look at BOVY: e.g., INTERPOLATE 2nd
#      Is it faster to split all functions between correct and don't correct?
#      upgrade _calc_surfacemass, _calc_sigma2surfacemass, and _calc_sigma2
#        to 'public' procedures; propagate to test_interpret_as_df.py
#      Allow more general surfacemass and sigma2 functions be specified
###############################################################################
_EPSREL=10.**-14.
_NSIGMA= 4.
import copy
import os, os.path
import cPickle as pickle
import math as m
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
           + DFcorrection kwargs (except for those already specified)
        OUTPUT:
        HISTORY:
            2010-03-10 - Written - Bovy (NYU)
        """
        self._dftype= dftype
        self._dfparams= dfparams
        self._beta= beta
        if dftype == 'corrected-dehnen':
            #Load corrections
            self._corr= DFcorrection(dfparams=dfparams,dftype='dehnen',
                                     beta=beta,**kwargs)
        return None

    def eval(self,E,L,log=False):
        """
        NAME:
           eval
        PURPOSE:
           evaluate the distribution function
        INPUT:
           E - energy (/vo^2)
           L - angular momentun (/ro/vo)
           log - return the log (default: false)
        OUTPUT:
           DF(E,L)
        HISTORY:
           2010-03-10 - Written - Bovy (NYU)
        """
        if self._dftype == 'dehnen':
            return self._eval_dehnen(E,L,log)
        elif self._dftype == 'corrected-dehnen':
            return self._eval_dehnen(E,L,log,correct=True)
    
    def _eval_dehnen(self,E,L,log,correct=False):
        """Internal function that evaluates the (corrected) Dehnen DF"""
        #Calculate Re,LE, OmegaE
        if self._beta == 0.:
            xE= sc.exp(E-.5)
            LE= xE
            OmegaE= 1./xE
        else: #non-flat rotation curve
            xE= (2.*E/(1.+1./self._beta))**(1./2./self._beta)
            LE= xE**(self._beta+1.)
            OmegaE= xE*(self._beta-1.)
        if correct:
            correction= self._corr.correct(xE)
        else:
            correction= sc.ones(2)
        SRE2= self._eval_SR2(xE)*correction[1]
        return 1./2./sc.pi*sc.sqrt(2./(1.+self._beta))*sc.exp(-xE/self._dfparams[0])/SRE2*sc.exp(OmegaE*(L-LE)/SRE2)*correction[0]
        
    def _eval_SR2(self,R,log=False):
        """Internal function that evaluates sigmaR^2 for an exponential profile"""
        if log:
            return 2.*sc.log(self._dfparams[2])-2.*(R-1.)/self._dfparams[1]
        else:
            return self._dfparams[2]**2.*sc.exp(-2.*(R-1.)/self._dfparams[1])

    def _eval_surfacemass(self,R,log=False):
        """Internal function that evaluates Sigma(R) for an exponential profile"""
        if log:
            return -R/self._dfparams[0]
        else:
            return sc.exp(-R/self._dfparams[0])

    def _calc_surfacemass(self,R,romberg=False,nsigma=None,relative=False):
        """Internal function that calculates the surface mass for a given DF at R"""
        if nsigma == None:
            nsigma= _NSIGMA
        logSigmaR= self._eval_surfacemass(R,log=True)
        sigmaR2= self._eval_SR2(R)
        sigmaR1= sc.sqrt(sigmaR2)
        logsigmaR2= sc.log(sigmaR2)
        gamma= sc.sqrt(2./(1.+self._beta))
        if relative:
            norm= 1.
        else:
            norm= sc.exp(logSigmaR)
        if romberg:
            return bovy_dblquad(_surfaceIntegrand,
                                gamma*R**self._beta/sigmaR1-nsigma,
                                gamma*R**self._beta/sigmaR1+nsigma,
                                lambda x: 0., lambda x: nsigma,
                                [R,self,logSigmaR,logsigmaR2,sigmaR1,gamma],
                                tol=10.**-8)/sc.pi*norm
        else:
            return integrate.dblquad(_surfaceIntegrand,
                                     gamma*R**self._beta/sigmaR1-nsigma,
                                     gamma*R**self._beta/sigmaR1+nsigma,
                                     lambda x: 0., lambda x: nsigma,
                                     (R,self,logSigmaR,logsigmaR2,sigmaR1,
                                      gamma),
                                     epsrel=_EPSREL)[0]/sc.pi*norm

    def _eval_surfaceIntegrand(self,E,L,logSigmaR,logsigmaR2):
        """Internal function that has the normalized DF for the surface mass integral"""
        if self._dftype == 'dehnen' or self._dftype == 'corrected-dehnen':
            return self._eval_surfaceIntegrand_dehnen(E,L,logSigmaR,logsigmaR2)
        
    def _eval_surfaceIntegrand_dehnen(self,E,L,logSigmaR,logsigmaR2):
        """Internal function that has the normalized DF for the surface mass integral for the uncorrected Dehnen DF"""
        #Calculate Re,LE, OmegaE
        if self._beta == 0.:
            xE= sc.exp(E-.5)
            logOLLE= sc.log(L/xE-1.)
            #LE= xE
            #OmegaE= 1./xE
        else: #non-flat rotation curve
            xE= (2.*E/(1.+1./self._beta))**(1./2./self._beta)
            logOLLE= self._beta*sc.log(xE)+sc.log(L/xE-xE**self._beta)
        if hasattr(self,'_corr'):
            correction= self._corr.correct(xE)
        else:
            correction= sc.ones(2)
        SRE2= self._eval_SR2(xE,log=True)+sc.log(correction[1])
        return sc.exp(logsigmaR2-SRE2+self._eval_surfacemass(xE,log=True)-logSigmaR+sc.exp(logOLLE-SRE2))*correction[0]

    def _calc_sigma2surfacemass(self,R,romberg=False,nsigma=None,
                                relative=False):
        """Internal function that calculates the velocity variance * surface 
        mass for a given DF at R"""
        if nsigma == None:
            nsigma= _NSIGMA
        logSigmaR= self._eval_surfacemass(R,log=True)
        sigmaR2= self._eval_SR2(R)
        sigmaR1= sc.sqrt(sigmaR2)
        logsigmaR2= sc.log(sigmaR2)
        gamma= sc.sqrt(2./(1.+self._beta))
        if relative:
            norm= 1.
        else:
            norm= sc.exp(logSigmaR+logsigmaR2)
        if romberg:
            return bovy_dblquad(_sigma2surfaceIntegrand,
                                gamma*R**self._beta/sigmaR1-nsigma,
                                gamma*R**self._beta/sigmaR1+nsigma,
                                lambda x: 0., lambda x: nsigma,
                                [R,self,logSigmaR,logsigmaR2,sigmaR1,gamma],
                                tol=10.**-8)/sc.pi*norm
        else:
            return integrate.dblquad(_sigma2surfaceIntegrand,
                                     gamma*R**self._beta/sigmaR1-nsigma,
                                     gamma*R**self._beta/sigmaR1+nsigma,
                                     lambda x: 0., lambda x: nsigma,
                                     (R,self,logSigmaR,logsigmaR2,sigmaR1,
                                      gamma),
                                     epsrel=_EPSREL)[0]/sc.pi*norm

    def _calc_sigma2(self,R,romberg=False,nsigma=None):
        """Shortcut to calculate sigma2(R)"""
        return self._calc_sigma2surfacemass(R,romberg,nsigma)/self._calc_surfacemass(R,romberg,nsigma)

def _surfaceIntegrand(vR,vT,R,df,logSigmaR,logsigmaR2,sigmaR1,gamma):
    """Internal function that is the integrand for the surface mass integration"""
    E,L= _vRpvTpRToEL(vR,vT,R,df._beta,sigmaR1,gamma)
    return df._eval_surfaceIntegrand(E,L,logSigmaR,logsigmaR2)

def _sigma2surfaceIntegrand(vR,vT,R,df,logSigmaR,logsigmaR2,sigmaR1,gamma):
    """Internal function that is the integrand for the sigma-squared times
    surface mass integration"""
    E,L= _vRpvTpRToEL(vR,vT,R,df._beta,sigmaR1,gamma)
    return vR**2.*df._eval_surfaceIntegrand(E,L,logSigmaR,logsigmaR2)

def _vRpvTpRToEL(vR,vT,R,beta,sigmaR1,gamma):
    """Internal function that calculates E and L given velocities normalized by the velocity dispersion"""
    vR*= sigmaR1
    vT*= sigmaR1/gamma
    return vRvTRToEL(vR,vT,R,beta)

def _oned_intFunc(x,twodfunc,gfun,hfun,tol,args):
    """Internal function for bovy_dblquad"""
    thisargs= copy.deepcopy(args)
    thisargs.insert(0,x)
    return integrate.romberg(twodfunc,gfun(x),hfun(x),args=thisargs,tol=tol)

def bovy_dblquad(func, a, b, gfun, hfun, args=(), tol=1.48e-08):
    """
    NAME:
       bovy_dblquad
    PURPOSE:
       like scipy.integrate's dblquad, but using Romberg integration for the one-d integrals and using tol
    INPUT:
       same as scipy.integrate.dblquad except for tol and epsrel,epsabs
    OUTPUT:
       value
    HISTORY:
       2010-03-11 - Written - Bpvy (NYU)
    """
    return integrate.romberg(_oned_intFunc,a,b,args=(func,gfun,hfun,tol,args),tol=tol)



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
            savefile= open(self._savefilename,'r')
            self._corrections= pickle.load(savefile)
            savefile.close()
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
        #nearR= round(R/self._rmax*(self._npoints-1.)) #BOVY: INTERPOLATE?
        #return self._corrections[nearR,:]
        thisR= R/self._rmax*(self._npoints-1.)
        lowR= m.floor(thisR)
        highR= m.ceil(thisR)
        #return self._corrections[lowR,:]+(thisR-lowR)/(highR-lowR)*(self._corrections[highR,:]-self._corrections[lowR,:]) #BOVY: CLEAN UP, interpolate log?
        return sc.exp(sc.log(self._corrections[lowR,:])+(thisR-lowR)/(highR-lowR)*(sc.log(self._corrections[highR,:])-sc.log(self._corrections[lowR,:]))) #BOVY: CLEAN UP, interpolate log?
    #BOVY: SECOND ORDER INTERPOLATE?

    def _calc_corrections(self):
        """Internal function that calculates the corrections"""     
        searchIter= self._niter-1
        while searchIter > 0:
            trySavefilename= os.path.join(self._savedir,'dfcorrection_'+
                                          self._dftype+
                                          '_%4.2f_%4.2f_%4.2f_%i_%4.2f_%i.sav'
                                          % (self._dfparams[0],
                                             self._dfparams[1],
                                             self._dfparams[2],self._npoints,
                                             self._rmax,searchIter))
            if os.path.exists(trySavefilename):
                trySavefile= open(trySavefilename,'r')
                corrections= pickle.load(trySavefile)
                trySavefile.close()
                break
            else:
                searchIter-= 1
        if searchIter == 0:
            corrections= sc.ones((self._npoints,2))
        rs= sc.linspace(1./10**8.,self._rmax,self._npoints)
        for ii in range(searchIter,self._niter):
            currentDF= distF(dftype='corrected-'+self._dftype,
                             dfparams=self._dfparams,
                             beta=self._beta,corrections=corrections)
            newcorrections= sc.zeros((self._npoints,2))
            for jj in range(self._npoints):
                thisSurface= currentDF._calc_surfacemass(rs[jj])
                newcorrections[jj,0]= currentDF._eval_surfacemass(rs[jj])/thisSurface
                newcorrections[jj,1]= currentDF._eval_SR2(rs[jj])*thisSurface/currentDF._calc_sigma2surfacemass(rs[jj])
                print jj, newcorrections[jj,:]
            corrections*= newcorrections
        #Save
        savefile= open(self._savefilename,'w')
        pickle.dump(corrections,savefile)
        savefile.close()
        return corrections
    
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


