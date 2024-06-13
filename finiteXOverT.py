# plot <X> over T

from latqcdtools.base.plotting import plt, getColorGradient, plot_dots, plot_lines, set_params
from latqcdtools.base.readWrite import readTable, writeTable
from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.math.spline import getSpline
import numpy as np
from latqcdtools.math.num_deriv import diff_deriv
from latqcdtools.statistics.bootstr import bootstr_from_gauss
from latqcdtools.base.printErrorBars import get_err_str
import latqcdtools.base.logger as logger
from scipy.stats import linregress

Nt = 12
Nss = [18, 24, 36, 48, 60]

Tcs = []
Tces = []

Ntcolors = getColorGradient(len(Nss))

icolor = 0
max_X = []
max_X_err = []
for Ns in Nss:
    data = readTable('Nt' + str(Nt) + '_Ns' + str(Ns) + '.txt')


    beta = data[0]
    X, Xe = data[3], data[4]
    max_X.append(X[4])
    max_X_err.append(Xe[4])



    T = []
    for b in beta:
        lp = latticeParams(Ns, Nt, b, scaleType='r0')
        T.append(lp.getT())
    tt = np.linspace(T[0], T[-1], 1001)


    def getTc(pm) -> float:
        """ Use a natural spline to find inflection point (Tc)

        Args:
            pm (array): Polyakov loop means

        Returns:
            float: Tc
        """
        spl = getSpline(T, pm, natural=True)
        dPdT = diff_deriv(tt, spl)
        maxIndex = np.argmax(dPdT)
        return tt[maxIndex]


    Tc, Tce = bootstr_from_gauss(getTc, data[1], data[2], 1000)

    # For each Nt, plot the spline, data, and estimated Tc band.
    Pspl = getSpline(T, X, natural=True)
    Tspl = np.linspace(T[0], T[-1], 1001)
    plot_lines(T, X, Xe, color=Ntcolors[icolor],label='$N_s=' + str(Ns) + '$')
    Tcs.append(Tc)
    Tces.append(Tce)

    logger.info(Ns, ":", get_err_str(Tc, Tce))
    icolor += 1


# Finish the plot
set_params(xlabel='$T$ [MeV]', legendpos=6, ylabel='$\\langle\\chi\\rangle$')
plt.savefig('Chi_overT.png')

print(max_X)
print(max_X_err)


