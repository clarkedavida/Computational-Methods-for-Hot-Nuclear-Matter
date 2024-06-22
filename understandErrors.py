# 
# understandErrors.py                                                               
# 
# D. Clarke 
# 
# Are jackknife errors very different from bootstrap errors? (Results seem to
# be more or less comparable as long as the time series is binned first.) 
# 

from latqcdtools.base.readWrite import readTable
from latqcdtools.statistics.jackknife import jackknife
from latqcdtools.statistics.bootstr import bootstr
from latqcdtools.statistics.statistics import std_mean, binSeries
from latqcdtools.base.plotting import plt, plot_dots, set_params
import latqcdtools.base.logger as logger

P, plaq = readTable('Nt16/Nt16_Ns48b65457.txt')

jack = []
boot = []

logger.info('Nmeas=',len(P))

NBS = [300,600,900,3000]
NBINS = [16,32,64,128,192,256,512]

for Nbootstrap in NBS:
    Pm, Pe = bootstr(std_mean,binSeries(P,128),Nbootstrap)
    boot.append(Pe)

for nbin in NBINS:
    Pm, Pe = jackknife(std_mean,P,nbin)
    jack.append(Pe)

fig, ax1 = plt.subplots()


# Plot data on the primary axis
plot_dots(NBS, boot, label='bootstrap',ax=ax1)
set_params(xlabel=r'$N_{\rm bootstrap}$',ax=ax1,legendpos=9)
ax2 = ax1.twiny()
plot_dots(NBINS, jack, ax=ax2,label='jackknife')
set_params(xlabel=r'$N_{\rm bins}$',ax=ax2,legendpos=7,xlabelpos=(.9,.9))
ax2.xaxis.set_ticks_position('top')
ax2.xaxis.set_label_position('top')


# Display the plot
plt.show()
