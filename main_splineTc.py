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
Nts = [6,8,10,12,14,16,18,20]

# there are two possible spline strategies. if you set Nknots=None, it will do a standard
# natural cubic spline. this is of course a solve to the data points, and will overfit to
# the extent that the spline will pass through every datum by definition. if you set
# Nknots to anything else, it will do what i call a 'natural-like' spline. my idea here was
# the following: i would like to do a spline fit that cares about error bars, but i would also
# like to impose zero curvature at the endpoints. this is achieved by creating a fake datum
# before the left-most point that is colinear with the left-most point and the second
# left-most point. the same is done on the right side. this strategy will favor zero
# curvature while still taking the error bars into account when fitting, which prevents
# any overfitting.
Nknots = 3

# this is where we're going to save our results for Tc and its error
Tcs, Tces, Bcs, Bces = [], [], [], []

Ntcolors=getColorGradient(len(Nts))

icolor = 0
for Nt in Nts:

    # latticeParams is an object that expects an Ns when it gets instantiated.
    # but actually it doesn't matter.
    Ns = Nt*3

    # INPUT: Here we need these tables that have beta, polyakov loop and error,
    #        susceptibility and error, Nconf, etc. We only use the polyakov loop
    #        and its error in this script.
    data=readTable(f'Nt{Nt}/Nt{Nt}_Ns{Ns}.txt')

    # translate from the "data" array to beta, polyakov loop mean, and its error
    beta, P, Pe = data[0][1:-1], data[1][1:-1], data[2][1:-1]


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
    plot_lines(Tspl,Pspl(Tspl),marker=None,color=Ntcolors[icolor],label='$N_\\tau='+str(Nt)+'$')
    plot_dots(T,P,Pe,color=Ntcolors[icolor])
    plt.axvspan(Tc-Tce,Tc+Tce,color=Ntcolors[icolor],alpha=0.3)
    Tcs.append(Tc)
    Tces.append(Tce)
    Bcs.append(Bc)
    Bces.append(Bce)

    logger.info(Nt,f": beta_c ={get_err_str(Bc,Bce)}  T_c = {get_err_str(Tc,Tce)}")
    icolor += 1

# Write out to table so we can do continuum-limit extrapolations
writeTable('cont_extrap.d',Nts,Tcs,Tces,header=['Nt','Tc [MeV]','err_Tc [Mev]'])
writeTable('cont_extrap_beta.d',Nts,Bcs,Bces,header=['Nt','betac','err_beta'])

# Finish the plot    
set_params(xlabel='$T$ [MeV]',legendpos=6,ylabel='$\\langle|P|\\rangle$')
plt.savefig('otherFigures/splineTcs.pdf')
plt.show()  