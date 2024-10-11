# 
# splineTc.py 
# 
# D. Clarke
#
# Find Tc using Polyakov loop inflection point. 
# 

from latqcdtools.base.plotting import plt, getColorGradient, plot_dots, plot_lines, set_params 
from latqcdtools.base.readWrite import readTable, writeTable
from latqcdtools.physics.lattice_params import latticeParams 
from latqcdtools.math.spline import getSpline 
import numpy as np
from latqcdtools.math.num_deriv import diff_deriv
from latqcdtools.statistics.bootstr import bootstr_from_gauss
from latqcdtools.base.printErrorBars import get_err_str
import latqcdtools.base.logger as logger
from latqcdtools.statistics.statistics import dev_by_dist

# make a list of the Ntau i want to look at
Nt  = 12
Nss = [18, 24, 36, 48, 60, 72, 84]

NknotsDict = {
18: 2,#3,
24: 2,#3,
36: 2,#3,
48: 2,#None, #3,
60: 2,#None, #3,
72: 2,#None, #3,
84: 2,#None, #3,
}
#Nknots = None 

xminDict = {
18: 2,
24: 2,
36: 2,
48: 2,
60: 2,#2,
72: 2,#2,
84: 2,#2,
}

xmaxDict = {
18: 1,
24: 1,
36: 1,
48: 1,
60: 1,
72: 1,
84: 1,
}

# this is where we're going to save our results for Tc and its error
Tcs, Tces, Bcs, Bces = [], [], [], []

Ntcolors=getColorGradient(len(Nss),map='winter')

icolor = 0
for Ns in Nss:

    Nknots = NknotsDict[Ns]
    xmin   = xminDict[Ns]
    xmax   = xmaxDict[Ns]


    # INPUT: Here we need these tables that have beta, polyakov loop and error,
    #        susceptibility and error, Nconf, etc. We only use the polyakov loop
    #        and its error in this script.
    data=readTable(f'Nt{Nt}/Nt{Nt}_Ns{Ns}.txt')

    # translate from the "data" array to beta, polyakov loop mean, and its error
    if xmax is None:
        beta, P, Pe = data[0][xmin:]     , data[1][xmin:]     , data[2][xmin:]
    else:
        beta, P, Pe = data[0][xmin:-xmax], data[1][xmin:-xmax], data[2][xmin:-xmax]


    # Convert beta to T using the same parameterization as in the literature.
    T = []
    for b in beta:
        lp = latticeParams(Ns,Nt,b,scaleType='r0',scaleYear=2014,paramYear=2015)
        T.append( lp.getT() )
    tt = np.linspace(T[0]   ,T[-1]   ,1000)
    bb = np.linspace(beta[0],beta[-1],1000)

    def getTc(pm,pe=None) -> float:
        """ Use a natural spline to find inflection point (Tc)

        Args:
            pm (array): Polyakov loop means 
            pe (array): Polyakov loop errors 

        Returns:
            float: Tc 
        """
        if Nknots is None:
            spl  = getSpline(xdata=T,ydata=pm,natural=True)
        else:
            spl  = getSpline(xdata=T,ydata=pm,edata=pe,num_knots=Nknots,natural=True)
        dPdT = diff_deriv(tt,spl)
        maxIndex = np.argmax(dPdT)
        return tt[maxIndex]

    def getBetac(pm,pe=None) -> float:
        """ Use a natural spline to find inflection point (beta_c)

        Args:
            pm (array): Polyakov loop means 

        Returns:
            float: beta_c 
        """
        if Nknots is None:
            spl  = getSpline(xdata=beta,ydata=pm,natural=True)
        else:
            spl  = getSpline(xdata=beta,ydata=pm,edata=pe,num_knots=Nknots,natural=True)
        dPdT = diff_deriv(bb,spl)
        maxIndex = np.argmax(dPdT)
        return bb[maxIndex]

    # We tried varying the number of bootstrap samples between 300 and 1000 and found
    # no significant changes. i set err_by_dist=False here, which then computes the
    # error using a typical standard deviation rather than using the inner 68-percentiles.
    # i had to do this because otherwise the bootstrap was finding zero error bar
    Tc, Tce = bootstr_from_gauss(getTc   ,P,Pe,1000,args=(Pe,),err_by_dist=False)
    Bc, Bce = bootstr_from_gauss(getBetac,P,Pe,1000,args=(Pe,),err_by_dist=False)

    # For each Nt, plot the spline, data, and estimated Tc band.
    if Nknots is None:
        Pspl = getSpline(xdata=T,ydata=P,natural=True)
    else:
        Pspl = getSpline(xdata=T,ydata=P,edata=Pe,num_knots=Nknots,natural=True)
    Tspl = np.linspace(T[0],T[-1],1001)
    plot_lines(Tspl,Pspl(Tspl),marker=None,color=Ntcolors[icolor],label='$N_s='+str(Ns)+'$')
    plot_dots(T,P,Pe,color=Ntcolors[icolor])
    plt.axvspan(Tc-Tce,Tc+Tce,color=Ntcolors[icolor],alpha=0.3)
    Tcs.append(Tc)
    Tces.append(Tce)
    Bcs.append(Bc)
    Bces.append(Bce)

    logger.info(Ns,f": beta_c ={get_err_str(Bc,Bce)}  T_c = {get_err_str(Tc,Tce)}")
    icolor += 1


# Finish the plot    
set_params(xlabel='$T$ [MeV]',legendpos=6,ylabel='$\\langle|P|\\rangle$')
plt.savefig('otherFigures/splineTcs_varyNs.pdf')
plt.show()  