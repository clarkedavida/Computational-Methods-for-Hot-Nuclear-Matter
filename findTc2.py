# 
# findTc.py                                                               
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

# make a list of the Ntau i want to look at
#Nts = [6,8,10,12,14,16,18,20]
Nts = [6,8,10,12,16,18,20]

# this is where we're going to save our results for Tc and its error
Tcs, Tces, Bcs, Bces = [], [], [], []

Ntcolors=getColorGradient(len(Nts))

icolor = 0
for Nt in Nts:

    # latticeParams is an object that expects an Ns when it gets instantiated.
    # but actually it doesn't matter.
    Ns = Nt*3
#    if Nt==12:
#        Ns = Nt*2

    # INPUT: Here we need these tables that have beta, polyakov loop and error,
    #        susceptibility and error, Nconf, etc. We only use the polyakov loop
    #        and its error in this script.
    data=readTable('Nt'+str(Nt)+'/Nt'+str(Nt)+'_Ns'+str(Ns)+'.txt')

    # translate from the "data" array to beta, polyakov loop mean, and its error
    beta  = data[0]
    P, Pe = data[1], data[2]

    T = []
    for b in beta:
        lp = latticeParams(Ns,Nt,b,scaleType='r0',scaleYear=2014,paramYear=2015)
        T.append( lp.getT() )
    tt = np.linspace(T[0]   ,T[-1]   ,1001)
    bb = np.linspace(beta[0],beta[-1],1001)

    def getTc(pm) -> float:
        """ Use a natural spline to find inflection point (Tc)

        Args:
            pm (array): Polyakov loop means 

        Returns:
            float: Tc 
        """
        spl  = getSpline(T,pm,natural=True)
        dPdT = diff_deriv(tt,spl)
        maxIndex = np.argmax(dPdT)
        return tt[maxIndex]

    def getBetac(pm) -> float:
        """ Use a natural spline to find inflection point (beta_c)

        Args:
            pm (array): Polyakov loop means 

        Returns:
            float: beta_c 
        """
        spl  = getSpline(beta,pm,natural=True)
        dPdT = diff_deriv(bb,spl)
        maxIndex = np.argmax(dPdT)
        return bb[maxIndex]

    # We tried varying the number of bootstrap samples between 300 and 1000 and found
    # no significant changes
    Tc, Tce = bootstr_from_gauss(getTc,data[1],data[2],300)
    Bc, Bce = bootstr_from_gauss(getBetac,data[1],data[2],300)

    # For each Nt, plot the spline, data, and estimated Tc band.
    Pspl = getSpline(T,P,natural=True)
    Tspl = np.linspace(T[0],T[-1],1001)
    plot_lines(Tspl,Pspl(Tspl),marker=None,color=Ntcolors[icolor],label='$N_\\tau='+str(Nt)+'$')
    plot_dots(T,P,Pe,color=Ntcolors[icolor])
    plt.axvspan(Tc-Tce,Tc+Tce,color=Ntcolors[icolor],alpha=0.3)
    Tcs.append(Tc)
    Tces.append(Tce)
    Bcs.append(Bc)
    Bces.append(Bce)

#    logger.info(Nt,":",get_err_str(Tc,Tce))
    logger.info(Nt,":",get_err_str(Bc,Bce))
    icolor += 1

# Write out to table so we can do continuum-limit extrapolations
writeTable('cont_extrap.d',Nts,Tcs,Tces,header=['Nt','Tc [MeV]','err_Tc [Mev]'])
writeTable('cont_extrap_beta.d',Nts,Bcs,Bces,header=['Nt','betac','err_beta'])

# Finish the plot    
set_params(xlabel='$T$ [MeV]',legendpos=6,ylabel='$\\langle|P|\\rangle$')
plt.show()  