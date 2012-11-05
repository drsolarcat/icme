from pylab import *
from numpy import *
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist import SubplotHost
from mpl_toolkits.axisartist import GridHelperCurveLinear
from mpl_toolkits.axisartist import ParasiteAxesAuxTrans
from datetime import datetime

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

# plot Pt(A) for GSR
def plotGsrAPtFull(APt):
    figure()
    plot(APt['x'], APt['y']*1e9, 'k') # full Pt(A)
    # set the labels
    xlabel('A [Tm]')
    ylabel('Pt [nPa]')
    # switch on the grid
    grid(True)
    # save the plots in png and eps
    if toSave:
        savefig(resultsDir+'/eps/gsr_APt_full.eps', format='eps')
        savefig(resultsDir+'/png/gsr_APt_full.png', format='png')

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
def plotGsrBzMap(X, Y, Axy, Bz, Ab, Aa, Bx, By, xx, yy, zz):
    # draw negative contour lines in solid (dashed by default)
    rcParams['contour.negative_linestyle'] = 'solid'
    figure()
    # filled contour plot of Bz(x,y)
    cc = contourf(X/AU, Y/AU, transpose(reshape(Bz*1e9, (X.size,-1))), int((max(Bz)-min(Bz))*1e9/0.1))
    #cc.set_cmap('autumn')
    # remove gaps between the colored areas
    for c in cc.collections:
        c.set_antialiased(False)
    # line contour plot of A(x,y)
    contour(X/AU, Y/AU, transpose(reshape(Axy, (X.size,-1))), 40, colors='k')
    # border line of the flux rope (Ab)
    contour(X/AU, Y/AU, transpose(reshape(Axy, (X.size,-1))), levels=[Ab], colors='w', linewidths=5)

    # determine the central part of the flux rope
    cp = contour(X/AU, Y/AU, transpose(reshape(Axy, (X.size,-1))), levels=[0.99*Aa], colors='w', linewidths=0)
    p = cp.collections[0].get_paths()[0]
    v = p.vertices
    # plot the central point of the flux rope (Ac)
    plot(mean(v[:,0]), mean(v[:,1]), '.w', markersize=12)

    # plot (Bx,By) quiver plot
    quiver(X/AU, zeros(X.size), -Bx*1e9, By*1e9, units='xy', angles='xy')

    # axes projection
    dx = (max(X)-min(X))/AU/10;
    dy = (max(Y)-min(Y))/AU/10;
    xxx=[dot([1,0,0],xx), dot([1,0,0],yy), dot([1,0,0],zz)];
    yyy=[dot([0,1,0],xx), dot([0,1,0],yy), dot([0,1,0],zz)];
    zzz=[dot([0,0,1],xx), dot([0,0,1],yy), dot([0,0,1],zz)];
    plot((1+array([0,xxx[0]]))*dx, (max(Y)/AU/dy-1+array([0,xxx[1]]))*dy, '-c', linewidth=3)
    plot((1+array([0,yyy[0]]))*dx, (max(Y)/AU/dy-1+array([0,yyy[1]]))*dy, '-m', linewidth=3)
    plot((1+array([0,zzz[0]]))*dx, (max(Y)/AU/dy-1+array([0,zzz[1]]))*dy, '-y', linewidth=3)

    # tighten the axes
    axis('tight')

    # labels
    xlabel('X [AU]')
    ylabel('Y [AU]')

    # add a colorbar
    cb = colorbar(cc, pad=0.07)
    cb.locator = MaxNLocator(14)
    cb.update_ticks()
    cb.set_label('Bz [nT]')

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
    cb = colorbar(cc, pad=0.07)
    cb.locator = MaxNLocator(14)
    cb.update_ticks()
    cb.set_label(r"$\tilde{\mathcal{R}}$")

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

def plotBrot(year, month, day, hour, minute, second,
             theta, phi,
             year1fr, month1fr, day1fr, hour1fr, minute1fr, second1fr,
             year2fr, month2fr, day2fr, hour2fr, minute2fr, second2fr):

    dates = []

    for i in xrange(len(year)):
        dates.append(datetime(year[i], month[i], day[i],
                              hour[i], minute[i], second[i]))

    date1fr = datetime(year1fr, month1fr, day1fr, hour1fr, minute1fr, second1fr)
    date2fr = datetime(year2fr, month2fr, day2fr, hour2fr, minute2fr, second2fr)

#    major = HourLocator([0, 12])
    major = DayLocator()
    minor = HourLocator()
    majorFormat = DateFormatter('%Y-%m-%d')

    figure()
    subplots_adjust(hspace=0.001)

    labelx = -0.08

    ax1 = subplot(211)
    ax1.plot(dates, theta, 'k')
#    ax1.axvline(date1fr, color='g', linestyle='dashed')
#    ax1.axvline(date2fr, color='g', linestyle='dashed')
    ax1.xaxis.set_major_locator(major)
    ax1.xaxis.set_major_formatter(majorFormat)
    ax1.xaxis.set_minor_locator(minor)
    ax1.yaxis.set_major_locator(MaxNLocator(nbins=4, symmetric=True))
    ax1.yaxis.set_minor_locator(MultipleLocator(1))
    ax1.set_ylabel('theta')
    ax1.yaxis.set_label_coords(labelx, 0.5)

    ax2 = subplot(212)
    ax2.plot(dates, phi, 'k')
#    ax2.axvline(date1fr, color='g', linestyle='dashed')
#    ax2.axvline(date2fr, color='g', linestyle='dashed')
    ax2.xaxis.set_major_locator(major)
    ax2.xaxis.set_major_formatter(majorFormat)
    ax2.xaxis.set_minor_locator(minor)
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, symmetric=True))
    ax2.yaxis.set_minor_locator(MultipleLocator(1))
    ax2.set_ylabel('phi')
    ax2.yaxis.set_label_coords(labelx, 0.5)

    # switch off the labels on upper plots
    xticklabels = ax1.get_xticklabels()
    setp(xticklabels, visible=False)

# plot in-situ data
def plotData(year, month, day, hour, minute, second,
             B, Bx, By, Bz, Vp, Vx, Vy, Vz, Pth, Np, Tp, Vth, beta,
             year1mc, month1mc, day1mc, hour1mc, minute1mc, second1mc,
             year2mc, month2mc, day2mc, hour2mc, minute2mc, second2mc,
             year1fr, month1fr, day1fr, hour1fr, minute1fr, second1fr,
             year2fr, month2fr, day2fr, hour2fr, minute2fr, second2fr):

    dates = []

    for i in xrange(len(year)):
        dates.append(datetime(year[i], month[i], day[i],
                              hour[i], minute[i], second[i]))

    date1mc = datetime(year1mc, month1mc, day1mc, hour1mc, minute1mc, second1mc)
    date2mc = datetime(year2mc, month2mc, day2mc, hour2mc, minute2mc, second2mc)
    date1fr = datetime(year1fr, month1fr, day1fr, hour1fr, minute1fr, second1fr)
    date2fr = datetime(year2fr, month2fr, day2fr, hour2fr, minute2fr, second2fr)

#    major = HourLocator([0, 12])
    major = DayLocator()
    minor = HourLocator()
    majorFormat = DateFormatter('%Y-%m-%d')

    figure(figsize=[8,10])
    subplots_adjust(hspace=0.001)

    labelx = -0.08

    # magnetic field # nT
    ax1 = subplot(711)
    ax1.plot(dates, B*1e9, 'k')
#    ax1.plot(dates, Bx*1e9, 'r')
#    ax1.plot(dates, By*1e9, 'g')
#    ax1.plot(dates, Bz*1e9, 'b')
    ax1.axvline(date1mc, color='r')
    ax1.axvline(date2mc, color='r')
    ax1.axvline(date1fr, color='g', linestyle='dashed')
    ax1.axvline(date2fr, color='g', linestyle='dashed')
    ax1.xaxis.set_major_locator(major)
    ax1.xaxis.set_major_formatter(majorFormat)
    ax1.xaxis.set_minor_locator(minor)
    ax1.yaxis.set_major_locator(MaxNLocator(nbins=4, symmetric=True))
    ax1.yaxis.set_minor_locator(MultipleLocator(1))
    ax1.set_ylabel('B [nT]')
    ax1.yaxis.set_label_coords(labelx, 0.5)

    # proton bilk speed # km/s
    ax2 = subplot(712)
#    ax2.plot(dates, Vp*1e-3-mean(Vp)*1e-3, 'k')
    ax2.plot(dates, Vp*1e-3-450, 'k')
#    ax2.plot(dates, Vx*1e-3-mean(Vx)*1e-3, 'r')
    ax2.plot(dates, Vx*1e-3+450, 'r')
#    ax2.plot(dates, Vy*1e-3-mean(Vy)*1e-3, 'g')
    ax2.plot(dates, Vy*1e-3, 'g')
#    ax2.plot(dates, Vz*1e-3-mean(Vz)*1e-3, 'b')
    ax2.plot(dates, Vz*1e-3, 'b')
    ax2.axvline(date1mc, color='r')
    ax2.axvline(date2mc, color='r')
    ax2.axvline(date1fr, color='g', linestyle='dashed')
    ax2.axvline(date2fr, color='g', linestyle='dashed')
    ax2.xaxis.set_major_locator(major)
    ax2.xaxis.set_major_formatter(majorFormat)
    ax2.xaxis.set_minor_locator(minor)
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='upper'))
#    ax2.yaxis.set_minor_locator(MultipleLocator(50))
    ax2.set_ylabel('Vp [km/s]')
    ax2.yaxis.set_label_coords(labelx, 0.5)

    # plasma pressure # nPa
    ax3 = subplot(713)
    ax3.plot(dates, Pth*1e9, 'k')
    ax3.axvline(date1mc, color='r')
    ax3.axvline(date2mc, color='r')
    ax3.axvline(date1fr, color='g', linestyle='dashed')
    ax3.axvline(date2fr, color='g', linestyle='dashed')
    ax3.xaxis.set_major_locator(major)
    ax3.xaxis.set_major_formatter(majorFormat)
    ax3.xaxis.set_minor_locator(minor)
    ax3.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1,2,4], prune='upper'))
    ax3.set_ylabel('Pth [nPa]')
    ax3.yaxis.set_label_coords(labelx, 0.5)

    # proton density # cm^-3
    ax4 = subplot(714)
    ax4.plot(dates, Np*1e-6, 'k')
    ax4.axvline(date1mc, color='r')
    ax4.axvline(date2mc, color='r')
    ax4.axvline(date1fr, color='g', linestyle='dashed')
    ax4.axvline(date2fr, color='g', linestyle='dashed')
    ax4.xaxis.set_major_locator(major)
    ax4.xaxis.set_major_formatter(majorFormat)
    ax4.xaxis.set_minor_locator(minor)
    ax4.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))
    ax4.set_ylabel('Np [cm^-3]')
    ax4.yaxis.set_label_coords(labelx, 0.5)

    # proton temperature # K
    ax5 = subplot(715)
    ax5.semilogy(dates, Tp, 'k')
    ax5.axvline(date1mc, color='r')
    ax5.axvline(date2mc, color='r')
    ax5.axvline(date1fr, color='g', linestyle='dashed')
    ax5.axvline(date2fr, color='g', linestyle='dashed')
    ax5.xaxis.set_major_locator(major)
    ax5.xaxis.set_major_formatter(majorFormat)
    ax5.xaxis.set_minor_locator(minor)
    ax5.set_ylabel('Tp [K]')
    ax5.yaxis.set_label_coords(labelx, 0.5)

    # thermal speed # km/s
    ax6 = subplot(716)
    ax6.plot(dates, Vth*1e-3, 'k')
    ax6.axvline(date1mc, color='r')
    ax6.axvline(date2mc, color='r')
    ax6.axvline(date1fr, color='g', linestyle='dashed')
    ax6.axvline(date2fr, color='g', linestyle='dashed')
    ax6.xaxis.set_major_locator(major)
    ax6.xaxis.set_major_formatter(majorFormat)
    ax6.xaxis.set_minor_locator(minor)
    ax6.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))
    ax6.yaxis.set_minor_locator(MultipleLocator(5))
    ax6.set_ylabel('Vth [km/s]')
    ax6.yaxis.set_label_coords(labelx, 0.5)

    # plasma beta
    ax7 = subplot(717)
    ax7.semilogy(dates, beta, 'k')
    ax7.axvline(date1mc, color='r')
    ax7.axvline(date2mc, color='r')
    ax7.axvline(date1fr, color='g', linestyle='dashed')
    ax7.axvline(date2fr, color='g', linestyle='dashed')
    ax7.xaxis.set_major_locator(major)
    ax7.xaxis.set_major_formatter(majorFormat)
    ax7.xaxis.set_minor_locator(minor)
    ax7.set_ylabel('beta')
    ax7.yaxis.set_label_coords(labelx, 0.5)

    # switch off the labels on upper plots
    xticklabels = ax1.get_xticklabels()+ax2.get_xticklabels()+ \
                  ax3.get_xticklabels()+ax4.get_xticklabels()+ \
                  ax5.get_xticklabels()+ax6.get_xticklabels()
    setp(xticklabels, visible=False)

    # save
    if toSave:
        savefig(resultsDir+'/eps/data.eps', format='eps')
        savefig(resultsDir+'/png/data.png', format='png')

# plot simple 1D data
def plotData1D(dataY):
    figure()
    plot(dataY, '.k')
    # grid
    grid(True)

# show all pending plots
def showPlots():
    show()

# initialize the plotter
def initPlotter(newToSave, newResultsDir):
    global toSave
    global resultsDir
    toSave = newToSave
    resultsDir = newResultsDir

