from pylab import *
from numpy import *
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist import SubplotHost
from mpl_toolkits.axisartist import GridHelperCurveLinear
from mpl_toolkits.axisartist import ParasiteAxesAuxTrans

# astronomical unit
AU = 149597870700 # m

toSave = False # whether to save the plots
resultsDir = '' # where to save the plots

# plot Pt(A) for GSR
def plotGsrAPt(APtIn, APtOut, APtFit):
    figure()
    plot(APtFit['x'], APtFit['y']*1e9, 'k') # fitted Pt(A)
    plot(APtIn['x'], APtIn['y']*1e9, 'go') # inward Pt(A)
    plot(APtOut['x'], APtOut['y']*1e9, 'ro') # outward Pt(A)
    # set the labels
    xlabel('A [Tm]')
    ylabel('Pt [nPa]')
    # switch on the grid
    grid(True)
    # save the plots in png and eps
    if toSave:
        savefig(resultsDir+'/eps/gsr_APt.eps', format='eps')
        savefig(resultsDir+'/png/gsr_APt.png', format='png')

# plot Bz(A) for GSR
def plotGsrABz(ABz, ABzFit):
    figure()
    plot(ABzFit['x'], ABzFit['y']*1e9, 'k') # fitted Bz(A)
    plot(ABz['x'], ABz['y']*1e9, 'ro') # Bz(A)
    # labels
    xlabel('A [Tm]')
    ylabel('Bz [nT]')
    # grid
    grid(True)
    # save
    if toSave:
        savefig(resultsDir+'/eps/gsr_ABz.eps', format='eps')
        savefig(resultsDir+'/png/gsr_ABz.png', format='png')

# plot dPt/dA for GSR
def plotGsrAdPt(AdPt):
    figure()
    plot(AdPt['x'], AdPt['y']*1e9, 'k') # dPt/dA(A)
    # labels
    xlabel('A [Tm]')
    ylabel('dPt/dA [nPa/Tm]')
    # grid
    grid(True)
    # save
    if toSave:
        savefig(resultsDir+'/eps/gsr_AdPt.eps', format='eps')
        savefig(resultsDir+'/png/gsr_AdPt.png', format='png')

# plot Bz map for GSR
def plotGsrBzMap(X, Y, Axy, Bz, Ab):
    # draw negative contour lines in solid (dashed by default)
    rcParams['contour.negative_linestyle'] = 'solid'
    figure()
    # filled contour plot of Bz(x,y)
    cc = contourf(X/AU, Y/AU, transpose(reshape(Bz*1e9, (X.size,-1))), 600)
    # remove gaps between the colored areas
    for c in cc.collections:
        c.set_antialiased(False)
    # line contour plot of A(x,y)
    contour(X/AU, Y/AU, transpose(reshape(Axy, (X.size,-1))), 40, colors='k')
    # border line of the flux rope (Ab)
    contour(X/AU, Y/AU, transpose(reshape(Axy, (X.size,-1))), levels=[Ab], colors='w', linewidths=5)

    # determine the central part of the flux rope
    cp = contour(X/AU, Y/AU, transpose(reshape(Axy, (X.size,-1))), levels=[min(Axy)], colors='w', linewidths=5)
    p = cp.collections[0].get_paths()[0]
    v = p.vertices
    # plot the central point of the flux rope (Ac)
    plot(mean(v[:,0]), mean(v[:,1]), '.w', markersize=12)

    # tighten the axes
    axis('tight')

    # labels
    xlabel('X [AU]')
    ylabel('Y [AU]')

    # add a colorbar
    colorbar(cc, pad=0.07)

    # save
    if toSave:
        savefig(resultsDir+'/eps/gsr_BzMap.eps', format='eps')
        savefig(resultsDir+'/png/gsr_BzMap.png', format='png')

# plot residue map for GSR
def plotGsrResidue(theta, phi, residue, optTheta, optPhi, mvabTheta=None, mvabPhi=None, mvubTheta=None, mvubPhi=None):
    fig = figure()
    fig.clf()

    # some matplotlib setup stuff which I don't fully understand but it works
    tr = Affine2D().scale(pi/180., 1.) + PolarAxes.PolarTransform()
    extreme_finder = angle_helper.ExtremeFinderCycle(20, 20,
                                                     lon_cycle = 360,
                                                     lat_cycle = None,
                                                     lon_minmax = None,
                                                     lat_minmax = (0, inf),
                                                     )
    grid_locator1 = angle_helper.LocatorDMS(12)
    tick_formatter1 = angle_helper.FormatterDMS()
    grid_helper = GridHelperCurveLinear(tr,
                                        extreme_finder=extreme_finder,
                                        grid_locator1=grid_locator1,
                                        tick_formatter1=tick_formatter1
                                        )
    ax1 = SubplotHost(fig, 1, 1, 1, grid_helper=grid_helper)
    fig.add_subplot(ax1)
    ax1.axis["right"].major_ticklabels.set_visible(True)
    ax1.axis["top"].major_ticklabels.set_visible(True)
    ax1.axis["right"].get_helper().nth_coord_ticks=0
    ax1.axis["bottom"].get_helper().nth_coord_ticks=1

    # draw the filled contoured map in polar coordinates
    cc = ax1.contourf(transpose(mat(theta))*mat(cos(phi*pi/180)), transpose(mat(theta))*mat(sin(phi*pi/180)), 1/transpose(reshape(residue, (phi.size,-1))), 100)
    # remove gaps between the contour lines
    for c in cc.collections:
        c.set_antialiased(False)

    # show the MVAB direction
    if mvabTheta is not None and mvabPhi is not None:
        ax1.plot(mvabTheta*cos(mvabPhi*pi/180), mvabTheta*sin(mvabPhi*pi/180), 'sk', markersize=8)

    # show the MVUB direction
    if mvubTheta is not None and mvubPhi is not None:
        ax1.plot(mvubTheta*cos(mvubPhi*pi/180), mvubTheta*sin(mvubPhi*pi/180), 'dk', markersize=8)

    # show the optimal direction
    ax1.plot(optTheta*cos(optPhi*pi/180), optTheta*sin(optPhi*pi/180), '.k', markersize=15)

    # aspect and initial axes limits
    ax1.set_aspect(1.)
    ax1.set_xlim(-90, 90)
    ax1.set_ylim(-90, 90)

    # add grid
    ax1.grid(True)

    # add colobar
    colorbar(cc, pad=0.07)

    # save
    if toSave:
        savefig(resultsDir+'/eps/gsr_ResidualMap.eps', format='eps')
        savefig(resultsDir+'/png/gsr_ResidualMap.png', format='png')

# plot B rotation for MVA
def plotMvaBrot(Bx, By):
    figure()
    plot(Bx*1e9, By*1e9, '.k') # (Bx,By)
    # labels
    xlabel('Bx [nT]')
    ylabel('By [nT]')
    # grid
    grid(True)
    # save
    if toSave:
        savefig(resultsDir+'/eps/mva_Brot.eps', format='eps')
        savefig(resultsDir+'/png/mva_Brot.png', format='png')

# show all pending plots
def showPlots():
    show()

# initialize the plotter
def initPlotter(newToSave, newResultsDir):
    global toSave
    global resultsDir
    toSave = newToSave
    resultsDir = newResultsDir

