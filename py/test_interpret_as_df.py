###############################################################################
#   test_interpret_as_df.py: module that tests the interpret_as_df module
# ToDo: Make it work with new dehnenDF class
###############################################################################
import interpret_as_df as df
import scipy as sc
import bovy_plot as plot
import os, os.path
def test_calc_surfacemass(baseplotfilename,format='png',ngrid=101):
    """
    NAME:
       test_calc_surfacemass
    PURPOSE:
       test the calculation of the surface mass density
    INPUT:
       baseplotfilename - basefilename for plots (we add _1.png etc.)
    OPTIONAL INPUT:
       format - specify a different format from 'png'
       ngrid - number of points to sample the surface mass density on
    OUTPUT:
       outputs plots
    HISTORY:
       2010-03-23 - Written - Bovy (NYU)
    """
    #Test 1: beta= 0., sigma_0= 0.5, Rsigma= 3 x Rscale, Rscale= Ro/3.
    if not os.path.exists(basefilename+'_1.png'):
        distFunc= df.dehnenDF(profileParams=(1./3.,1.,0.5*sc.exp(-1.)),
                              beta=0.)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= sc.log(distFunc.surfacemass(xs[ii]))+xs[ii]*3.
        plot.bovy_print(fig_width=5,fig_height=2.5)
        plot.bovy_plot(xs*3.,sigma,'k',xlabel=r'R/R_s',
                       ylabel=r'\ln \Sigma_f - \ln \Sigma',xrange=[0,10],
                       yrange=[-0.25,0.25])
        plot.bovy_text(r'$\beta = 0.0\,,\quad \sigma_0 = 0.5 v_0\,,\quad R_{\sigma} = 3 R_s$',
                       bottom_right=True)
        plot.bovy_end_print(basefilename+'_1.'+format)
        print "Wrote "+basefilename+'_1.'+format
    else:
        print basefilename+'_1.png exists'
        print 'Move or delete this file to execute the test'


    #Test 2: beta= 0.2, sigma_0= 0.5, Rsigma= 3 x Rscale, Rscale= Ro/3.
    if not os.path.exists(basefilename+'_2.png'):
        distFunc= df.dehnenDF(profileParams=(1./3.,1.,0.5*sc.exp(-1.)),
                        beta=0.2)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= sc.log(distFunc.surfacemass(xs[ii]))+xs[ii]*3.
        plot.bovy_print(fig_width=5,fig_height=2.5)
        plot.bovy_plot(xs*3.,sigma,'k',xlabel=r'R/R_s',
                       ylabel=r'\ln \Sigma_f - \ln \Sigma',xrange=[0,10],
                       yrange=[-0.12,0.12])
        plot.bovy_text(r'$\beta = 0.2\,,\quad \sigma_0 = 0.5 v_0\,,\quad R_{\sigma} = 3 R_s$',
                       bottom_right=True)
        plot.bovy_end_print(basefilename+'_2.'+format)
        print "Wrote "+basefilename+'_2.'+format
    else:
        print basefilename+'_2.png exists'
        print 'Move or delete this file to execute the test'


    #Test 3: beta= -0.2, sigma_0= 0.5, Rsigma= 3 x Rscale, Rscale= Ro/3.
    if not os.path.exists(basefilename+'_3.png'):
        distFunc= df.dehnenDF(profileParams=(1./3.,1.,0.5*sc.exp(-1.)),
                        beta=-0.2)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= sc.log(distFunc.surfacemass(xs[ii]))+xs[ii]*3.
        plot.bovy_print(fig_width=5,fig_height=2.5)
        plot.bovy_plot(xs*3.,sigma,'k',xlabel=r'R/R_s',
                       ylabel=r'\ln \Sigma_f - \ln \Sigma',xrange=[0,10],
                       yrange=[-0.45,0.45])
        plot.bovy_text(r'$\beta = -0.2\,,\quad \sigma_0 = 0.5 v_0\,,\quad R_{\sigma} = 3 R_s$',
                       bottom_right=True)
        plot.bovy_end_print(basefilename+'_3.'+format)
        print "Wrote "+basefilename+'_3.'+format
    else:
        print basefilename+'_3.png exists'
        print 'Move or delete this file to execute the test'


    #Test 4: beta= 0., sigma_0= 1, Rsigma= 3 x Rscale, Rscale= Ro/3.
    if not os.path.exists(basefilename+'_4.png'):
        distFunc= df.dehnenDF(profileParams=(1./3.,1.,sc.exp(-1.)),
                        beta=0.)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= sc.log(distFunc.surfacemass(xs[ii]))+xs[ii]*3.
        plot.bovy_print(fig_width=5,fig_height=2.5)
        plot.bovy_plot(xs*3.,sigma,'k',xlabel=r'R/R_s',
                       ylabel=r'\ln \Sigma_f - \ln \Sigma',xrange=[0,10],
                       yrange=[-0.95,0.95])
        plot.bovy_text(r'$\beta = 0.0\,,\quad \sigma_0 = v_0\,,\quad R_{\sigma} = 3 R_s$',
                       bottom_right=True)
        plot.bovy_end_print(basefilename+'_4.'+format)
        print "Wrote "+basefilename+'_4.'+format
    else:
        print basefilename+'_4.png exists'
        print 'Move or delete this file to execute the test'



def test_calc_sigma(baseplotfilename,format='png',ngrid=101):
    """
    NAME:
       test_calc_sigma
    PURPOSE:
       test the calculation of the radial velocity dispersion
    INPUT:
       baseplotfilename - basefilename for plots (we add _1.png etc.)
    OPTIONAL INPUT:
       format - specify a different format from 'png'
       ngrid - number of points to sample the surface mass density on
    OUTPUT:
       outputs plots
    HISTORY:
       2010-03-23 - Written - Bovy (NYU)
    """
    #Test 1: beta= 0., sigma_0= 0.5, Rsigma= 3 x Rscale, Rscale= Ro/3.
    if not os.path.exists(basefilename+'_1.png'):
        distFunc= df.dehnenDF(profileParams=(1./3.,1.,0.5*sc.exp(-1.)),
                        beta=0.)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= 0.5*sc.log(distFunc.sigma2(xs[ii]))+xs[ii]-sc.log(0.5)
        plot.bovy_print(fig_width=5,fig_height=2.5)
        plot.bovy_plot(xs*3.,sigma,'k',xlabel=r'R/R_s',
                       ylabel=r'\ln \sigma_f - \ln \sigma',xrange=[0,10],
                       yrange=[-0.25,0.25])
        plot.bovy_text(r'$\beta = 0.0\,,\quad \sigma_0 = 0.5 v_0\,,\quad R_{\sigma} = 3 R_s$',
                       bottom_right=True)
        plot.bovy_end_print(basefilename+'_1.'+format)
        print "Wrote "+basefilename+'_1.'+format
    else:
        print basefilename+'_1.png exists'
        print 'Move or delete this file to execute the test'


    #Test 2: beta= 0.2, sigma_0= 0.5, Rsigma= 3 x Rscale, Rscale= Ro/3.
    if not os.path.exists(basefilename+'_2.png'):
        distFunc= df.dehnenDF(profileParams=(1./3.,1.,0.5*sc.exp(-1.)),
                        beta=0.2)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= 0.5*sc.log(distFunc.sigma2(xs[ii]))+xs[ii]-sc.log(0.5)
        plot.bovy_print(fig_width=5,fig_height=2.5)
        plot.bovy_plot(xs*3.,sigma,'k',xlabel=r'R/R_s',
                       ylabel=r'\ln \sigma_f - \ln \sigma',xrange=[0,10],
                       yrange=[-0.12,0.12])
        plot.bovy_text(r'$\beta = 0.2\,,\quad \sigma_0 = 0.5 v_0\,,\quad R_{\sigma} = 3 R_s$',
                       bottom_right=True)
        plot.bovy_end_print(basefilename+'_2.'+format)
        print "Wrote "+basefilename+'_2.'+format
    else:
        print basefilename+'_2.png exists'
        print 'Move or delete this file to execute the test'


    #Test 3: beta= -0.2, sigma_0= 0.5, Rsigma= 3 x Rscale, Rscale= Ro/3.
    if not os.path.exists(basefilename+'_3.png'):
        distFunc= df.dehnenDF(profileParams=(1./3.,1.,0.5*sc.exp(-1.)),
                        beta=-0.2)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= 0.5*sc.log(distFunc.sigma2(xs[ii]))+xs[ii]-sc.log(0.5)
        plot.bovy_print(fig_width=5,fig_height=2.5)
        plot.bovy_plot(xs*3.,sigma,'k',xlabel=r'R/R_s',
                       ylabel=r'\ln \sigma_f - \ln \sigma',xrange=[0,10],
                       yrange=[-0.45,0.45])
        plot.bovy_text(r'$\beta = -0.2\,,\quad \sigma_0 = 0.5 v_0\,,\quad R_{\sigma} = 3 R_s$',
                       bottom_right=True)
        plot.bovy_end_print(basefilename+'_3.'+format)
        print "Wrote "+basefilename+'_3.'+format
    else:
        print basefilename+'_3.png exists'
        print 'Move or delete this file to execute the test'


    #Test 4: beta= 0., sigma_0= 1, Rsigma= 3 x Rscale, Rscale= Ro/3.
    if not os.path.exists(basefilename+'_4.png'):
        distFunc= df.dehnenDF(profileParams=(1./3.,1.,sc.exp(-1.)),
                        beta=0.)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= 0.5*sc.log(distFunc.sigma2(xs[ii]))+xs[ii]
        plot.bovy_print(fig_width=5,fig_height=2.5)
        plot.bovy_plot(xs*3.,sigma,'k',xlabel=r'R/R_s',
                       ylabel=r'\ln \sigma_f - \ln \sigma',xrange=[0,10],
                       yrange=[-0.95,0.95])
        plot.bovy_text(r'$\beta = 0.0\,,\quad \sigma_0 = v_0\,,\quad R_{\sigma} = 3 R_s$',
                       bottom_right=True)
        plot.bovy_end_print(basefilename+'_4.'+format)
        print "Wrote "+basefilename+'_4.'+format
    else:
        print basefilename+'_4.png exists'
        print 'Move or delete this file to execute the test'

def test_surfacemass_corrections(baseplotfilename,format='png'):
    """
    NAME:
       test_surfacemass_corrections
    PURPOSE:
       test the calculation of the surface mass density corrections
    INPUT:
       baseplotfilename - basefilename for plots
    OPTIONAL INPUT:
       format - specify a different format from 'png'
    OUTPUT:
       outputs plots
    HISTORY:
       2010-03-29 - Written - Bovy (NYU)
    """
    savedir='../corrections/'
    #Test 1: beta= 0., sigma_R0= 0.2, Rsigma= 3 x Rscale, Rscale= Ro/3.
    if not os.path.exists(basefilename+'_sigma0_0.5.png'):
        profileParams= (1./3.,1.,0.2)
        xs= sc.linspace(0,15,151)
        df1= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=1)
        df2= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=2)
        df3= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=3)
        df4= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=4)
        df5= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=5)
        df9= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=9)
        df10= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=10)
        df14= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=14)
        df15= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=15)
        df19= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=19)
        df20= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=20)
        plot.bovy_print(fig_width=10,fig_height=4.5)
        plot.bovy_plot(xs,-sc.log(df1._corr._corrections[:,0]),
                       'k',xlabel=r'R/R_s',
                       ylabel=r'\ln \Sigma_f - \ln \Sigma',xrange=[0,15],
                       yrange=[-0.015,0.015])
        plot.bovy_plot(xs,-sc.log(df2._corr._corrections[:,0]/df1._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df3._corr._corrections[:,0]/df2._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df4._corr._corrections[:,0]/df3._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df5._corr._corrections[:,0]/df4._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df10._corr._corrections[:,0]/df9._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df15._corr._corrections[:,0]/df14._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df20._corr._corrections[:,0]/df19._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_text(r'$\beta = 0.0\,,\quad \sigma_{R}(3\,R_s) = 0.2 v_0\,,\quad R_{\sigma} = 3 R_s$',
                       bottom_right=True)
        plot.bovy_end_print(basefilename+'_sigma0_0.5.'+format)
        print "Wrote "+basefilename+'_sigma0_0.5.'+format
    else:
        print basefilename+'_sigma0_0.5'+format+' exists'
        print 'Move or delete this file to execute the test'

    #Test 2: beta= 0., sigma_0= 1, Rsigma= 3 x Rscale, Rscale= Ro/3.
    if not os.path.exists(basefilename+'_sigma0_1.0.png'):
        profileParams= (1./3.,1.,sc.exp(-1.))
        xs= sc.linspace(0,15,151)
        df1= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=1)
        df2= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=2)
        df3= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=3)
        df4= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=4)
        df5= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=5)
        df9= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=9)
        df10= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=10)
        df14= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=14)
        df15= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=15)
        df19= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=19)
        df20= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=20)
        plot.bovy_print(fig_width=10,fig_height=4.5)
        plot.bovy_plot(xs,-sc.log(df1._corr._corrections[:,0]),
                       'k',xlabel=r'R/R_s',
                       ylabel=r'\ln \Sigma_f - \ln \Sigma',xrange=[0,15],
                       yrange=[-0.015,0.015])
        plot.bovy_plot(xs,-sc.log(df2._corr._corrections[:,0]/df1._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df3._corr._corrections[:,0]/df2._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df4._corr._corrections[:,0]/df3._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df5._corr._corrections[:,0]/df4._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df10._corr._corrections[:,0]/df9._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df15._corr._corrections[:,0]/df14._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df20._corr._corrections[:,0]/df19._corr._corrections[:,0]),
                       overplot=True)
        plot.bovy_text(r'$\beta = 0.0\,,\quad \sigma_0 = v_0\,,\quad R_{\sigma} = 3 R_s$',
                       bottom_right=True)
        plot.bovy_end_print(basefilename+'_sigma0_1.0.'+format)
        print "Wrote "+basefilename+'_sigma0_1.0.'+format
    else:
        print basefilename+'_sigma0_1.0'+format+' exists'
        print 'Move or delete this file to execute the test'


def test_sigma_corrections(baseplotfilename,format='png'):
    """
    NAME:
       test_sigma_corrections
    PURPOSE:
       test the calculation of the sigma_R corrections
    INPUT:
       baseplotfilename - basefilename for plots
    OPTIONAL INPUT:
       format - specify a different format from 'png'
    OUTPUT:
       outputs plots
    HISTORY:
       2010-03-29 - Written - Bovy (NYU)
    """
    savedir='../corrections/'
    #Test 1: beta= 0., sigma_R0= 0.2, Rsigma= 3 x Rscale, Rscale= Ro/3.
    if not os.path.exists(basefilename+'_sigma0_0.5.png'):
        profileParams= (1./3.,1.,0.2)
        xs= sc.linspace(0,15,151)
        df1= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=1)
        df2= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=2)
        df3= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=3)
        df4= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=4)
        df5= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=5)
        df9= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=9)
        df10= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=10)
        df14= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=14)
        df15= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=15)
        df19= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=19)
        df20= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=20)
        plot.bovy_print(fig_width=10,fig_height=4.5)
        plot.bovy_plot(xs,-sc.log(df1._corr._corrections[:,1])/2.,
                       'k',xlabel=r'R/R_s',
                       ylabel=r'\ln \Sigma_f - \ln \Sigma',xrange=[0,15],
                       yrange=[-0.015,0.015])
        plot.bovy_plot(xs,-sc.log(df2._corr._corrections[:,1]/df1._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df3._corr._corrections[:,1]/df2._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df4._corr._corrections[:,1]/df3._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df5._corr._corrections[:,1]/df4._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df10._corr._corrections[:,1]/df9._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df15._corr._corrections[:,1]/df14._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df20._corr._corrections[:,1]/df19._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_text(r'$\beta = 0.0\,,\quad \sigma_{R}(3\,R_s) = 0.2 v_0\,,\quad R_{\sigma} = 3 R_s$',
                       bottom_right=True)
        plot.bovy_end_print(basefilename+'_sigma0_0.5.'+format)
        print "Wrote "+basefilename+'_sigma0_0.5.'+format
    else:
        print basefilename+'_sigma0_0.5'+format+' exists'
        print 'Move or delete this file to execute the test'

    #Test 2: beta= 0., sigma_0= 1, Rsigma= 3 x Rscale, Rscale= Ro/3.
    if not os.path.exists(basefilename+'_sigma0_1.0.png'):
        profileParams= (1./3.,1.,sc.exp(-1.))
        xs= sc.linspace(0,15,151)
        df1= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=1)
        df2= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=2)
        df3= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=3)
        df4= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=4)
        df5= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                        beta=0.,niter=5)
        df9= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=9)
        df10= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=10)
        df14= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=14)
        df15= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=15)
        df19= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=19)
        df20= df.dehnenDF(profileParams=profileParams,savedir=savedir,
                         beta=0.,niter=20)
        plot.bovy_print(fig_width=10,fig_height=4.5)
        plot.bovy_plot(xs,-sc.log(df1._corr._corrections[:,1])/2.,
                       'k',xlabel=r'R/R_s',
                       ylabel=r'\ln \Sigma_f - \ln \Sigma',xrange=[0,15],
                       yrange=[-0.015,0.015])
        plot.bovy_plot(xs,-sc.log(df2._corr._corrections[:,1]/df1._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df3._corr._corrections[:,1]/df2._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df4._corr._corrections[:,1]/df3._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df5._corr._corrections[:,1]/df4._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df10._corr._corrections[:,1]/df9._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df15._corr._corrections[:,1]/df14._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_plot(xs,-sc.log(df20._corr._corrections[:,1]/df19._corr._corrections[:,1])/2.,
                       overplot=True)
        plot.bovy_text(r'$\beta = 0.0\,,\quad \sigma_0 = v_0\,,\quad R_{\sigma} = 3 R_s$',
                       bottom_right=True)
        plot.bovy_end_print(basefilename+'_sigma0_1.0.'+format)
        print "Wrote "+basefilename+'_sigma0_1.0.'+format
    else:
        print basefilename+'_sigma0_1.0'+format+' exists'
        print 'Move or delete this file to execute the test'



if __name__ == '__main__':
    #Test surface-mass
    print "Testing surface-mass calculation ..."
    basefilename= '../dfTest/testSurfacemass'
    test_calc_surfacemass(basefilename,ngrid=101)

    #Test sigma
    print "Testing sigma calculation ..."
    basefilename= '../dfTest/testSigma'
    test_calc_sigma(basefilename,ngrid=101)

    #Test surfacemass corrections
    print "Testing surface-mass corrections ..."
    basefilename= '../dfTest/testSurfacemassCorrections'
    test_surfacemass_corrections(basefilename)

    #Test surfacemass corrections
    print "Testing sigma corrections ..."
    basefilename= '../dfTest/testSigmaCorrections'
    test_sigma_corrections(basefilename)
