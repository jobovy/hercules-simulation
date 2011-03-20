###############################################################################
#   apogee_figures.py: make nice figures for the bar prediction in APOGEE
###############################################################################
import os, os.path
import sys
import cPickle as pickle
import math as m
import scipy as sc
from scipy import signal
from optparse import OptionParser
import galpy.util.bovy_plot as plot
from matplotlib import pyplot
import matplotlib.ticker as ticker
import galpy.util.bovy_plot as bovy_plot
#from calc_veldist_2d import calc_veldist_2d
from calc_veldist_1d import predictVlos, predictVlosConvolve
_degtorad= m.pi/180.
_radtodeg= 180./m.pi
_DEFAULTL= 235.
_VLOSGRID=26
_DGRID= 21
def apogee_figures(plotfilename,savefilename=None,bar_angle=25.,dt=None,
                   l=None,rolr=None,bar_strength=None,slope=None,
                   vlosgrid=201,dgrid=101,
                   conditional=False):
    """
    NAME:
       apogee_figures
    PURPOSE:
       make a vlos-d plot for APOGEE
    INPUT:
       plotfilename - name of the file the figure will be saved to
       savefilename - name of the file the velocity distributions will
                      be saved to
       bar_angle= angle between the GC-Sun lin and the bar major axis
       dt= - time to integrate for (in bar-periods)
       l= - Galactic longitude
       rolr= radius of outer lindblad radius
       bar_strength= strength of the bar
       slop= slope of the rotation curve (power-law index)
       conditional= if True, normalize each velocity distribution independently
    OUTPUT:
       saves velocity distributions to pickle save file and produces plot
    HISTORY:
       2011-03-19 - Written - Bovy (NYU)
    """
    #Grid
    vlosgrid=_VLOSGRID
    dgrid= _DGRID
    vloslinspace= (-.9,.9,vlosgrid)
    vloss= sc.linspace(*vloslinspace)
    dlinspace= (0.0001,10./8.,dgrid)
    ds= sc.linspace(*dlinspace)

    #Set up parameters
    potparams= (rolr,bar_strength,bar_angle*_degtorad,.8,None)

    if os.path.exists(savefilename):
        savefile= open(savefilename,"rb")
        vlosds= pickle.load(savefile)
        dd= pickle.load(savefile)
        savefile.close()
    else:
        vlosds= []
        dd= 0

    while dd < dgrid:
        print "Working on %i / %i ..." % (dd+1,dgrid)
        #Calculate vlos for this distance
        if dt is None:
            vlosd= predictVlos(vloslinspace,
                               l=(l*_degtorad),
                               d=ds[dd],
                               distCoord='Sun',
                               pot='bar',beta=slope,
                               potparams=potparams)
        else:
            vlosd= predictVlos(vloslinspace,
                               l=(l*_degtorad),
                               d=ds[dd],
                               distCoord='Sun',
                               pot='bar',beta=slope,
                               potparams=potparams,t=dt)
        vlosds.append(vlosd)
        dd+= 1
        #Save
        savefile= open(savefilename,"wb")
        pickle.dump(vlosds,savefile)
        pickle.dump(dd,savefile)
        savefile.close()

    #Plot
    if conditional:
        newvlosds= []
        for vlosd in vlosds:
            newvlosds.append(vlosd/(sc.nansum(vlosd)*(vloss[1]-vloss[0])))
        vlosds= newvlosds
    bovy_plot.bovy_print()
    bovy_plot.bovy_dens2d(sc.array(vlosds).T,
                          origin='lower',cmap='gist_yarg',aspect=0.66,
                          xlabel=r'$d\ /\ R_0$',
                          ylabel=r'$(v_{\mathrm{los}} - \vec{v}_c \cdot \vec{d})\ /\ v_0$',
                          yrange=sc.array([vloslinspace[0],vloslinspace[1]]),
                          xrange=sc.array([dlinspace[0],dlinspace[1]]),
                          contours=False,cntrmass=False)
    if bar_strength == 0.:
        bovy_plot.bovy_text(r'$l = %i^\circ$' % int(l),
                            top_right=True)
    else:
        bovy_plot.bovy_text(r'$l = %i^\circ$' % int(l) +'\n'+
                            r'$R_{\mathrm{OLR}} = %3.1f$' % rolr +'\n'+
                            r'$\alpha = %5.3f$' % bar_strength+'\n'+
                            r'$\phi_{\mathrm{bar}} = %i^\circ$'
                            % int(bar_angle),
                            top_right=True)
    bovy_plot.bovy_end_print(plotfilename)


def get_options():
    usage = "usage: %prog [options] <plotfilename>\n\nplotfilename= name of the file that the figure will be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("-s","--savefile",dest="savefilename",default=None,
                      help="Name of the file the solution prediction will be saved to")
#    parser.add_option("-c", "--col",dest="col",type='int',
#                      default=None,
#                      help="Only calculate one column of los")
    parser.add_option("-l",dest="l",type='float',
                      default=_DEFAULTL,
                      help="Galactic longitude to consider in degree")
    parser.add_option("--convolve", dest="convolve",
                      default=False,action="store_true", 
                      help="Make plot showing influence of distance uncertainties")
    parser.add_option("--vrconvolve", dest="vrconvolve",
                      default=False,action="store_true", 
                      help="Make plot showing influence of los velocity uncertainties")
#    parser.add_option("--df", dest="df",
#                      default=False,action="store_true", 
#                      help="Make plot showing influence of distribution function assumptions")
    parser.add_option("--rolr", dest="rolr",type='float',
                      default=0.9,
                      help="R_OLR")
    parser.add_option("--barstrength", dest="barstrength",type='float',
                      default=0.01,
                      help="Bar strength")
    parser.add_option("--barangle", dest="barangle",type='float',
                      default=25.,
                      help="Bar angle [deg]")
    parser.add_option("--dt", dest="dt",type=float,
                      default=None,
                      help="Time to integrate for")
    parser.add_option("--slope", dest="slope",type='float',
                      default=0.,
                      help="slope of the rotation curve (power-law)")
    parser.add_option("--conditional", dest="conditional",
                      default=False,action="store_true", 
                      help="Normalize all of the velocity distributions indepently")

    return parser

if __name__ == '__main__':
    parser= get_options()
    (options,args)= parser.parse_args()
    if len(args) == 0 or options.savefilename is None:
        parser.print_help()
        import sys
        sys.exit(-1)
    if options.dt == 'None' or options.dt is None:
        dt= None
    else:
        dt= float(options.dt)
    apogee_figures(args[0],savefilename=options.savefilename,
                   dt=dt,bar_angle=options.barangle,
                   l=options.l,rolr=options.rolr,
                   bar_strength=options.barstrength,slope=options.slope,
                   conditional=options.conditional)
    
