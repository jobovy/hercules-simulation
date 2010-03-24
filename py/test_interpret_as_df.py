###############################################################################
#   test_interpret_as_df.py: module that tests the interpret_as_df module
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
        distF= df.distF(dftype='dehnen',dfparams=(1./3.,1.,0.5*sc.exp(-1.)),
                        beta=0.)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= sc.log(distF._calc_surfacemass(xs[ii]))+xs[ii]*3.
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
        distF= df.distF(dftype='dehnen',dfparams=(1./3.,1.,0.5*sc.exp(-1.)),
                        beta=0.2)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= sc.log(distF._calc_surfacemass(xs[ii]))+xs[ii]*3.
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
        distF= df.distF(dftype='dehnen',dfparams=(1./3.,1.,0.5*sc.exp(-1.)),
                        beta=-0.2)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= sc.log(distF._calc_surfacemass(xs[ii]))+xs[ii]*3.
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
        distF= df.distF(dftype='dehnen',dfparams=(1./3.,1.,sc.exp(-1.)),
                        beta=0.)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= sc.log(distF._calc_surfacemass(xs[ii]))+xs[ii]*3.
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
        distF= df.distF(dftype='dehnen',dfparams=(1./3.,1.,0.5*sc.exp(-1.)),
                        beta=0.)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= 0.5*sc.log(distF._calc_sigma2(xs[ii]))+xs[ii]-sc.log(0.5)
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
        distF= df.distF(dftype='dehnen',dfparams=(1./3.,1.,0.5*sc.exp(-1.)),
                        beta=0.2)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= 0.5*sc.log(distF._calc_sigma2(xs[ii]))+xs[ii]-sc.log(0.5)
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
        distF= df.distF(dftype='dehnen',dfparams=(1./3.,1.,0.5*sc.exp(-1.)),
                        beta=-0.2)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= 0.5*sc.log(distF._calc_sigma2(xs[ii]))+xs[ii]-sc.log(0.5)
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
        distF= df.distF(dftype='dehnen',dfparams=(1./3.,1.,sc.exp(-1.)),
                        beta=0.)
        xs= sc.linspace(0.00001,10./3.,ngrid)
        sigma= sc.zeros(ngrid)
        for ii in range(ngrid):
            sigma[ii]= 0.5*sc.log(distF._calc_sigma2(xs[ii]))+xs[ii]
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




if __name__ == '__main__':
    #Test surface-mass
    print "Testing surface-mass calculation ..."
    basefilename= '../dfTest/testSurfacemass'
    test_calc_surfacemass(basefilename,ngrid=101)

    #Test sigma
    print "Testing sigma calculation ..."
    basefilename= '../dfTest/testSigma'
    test_calc_sigma(basefilename,ngrid=101)
