# 
# continuumExtrapolate.py                                                               
# 
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.readWrite import readTable
from latqcdtools.statistics.statistics import gaudif, error_prop, BAIC, modelAverage
from latqcdtools.physics.continuumExtrap import continuumExtrapolate, _powerSeries
from latqcdtools.physics.constants import r0_phys
from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.physics.referenceScales import r0_div_a, ignoreBetaRange
from latqcdtools.base.printErrorBars import get_err_str
from latqcdtools.base.plotting import clearPlot, set_params
from latqcdtools.base.utilities import toNumpy
from latqcdtools.base.check import ignoreUnderflow

ignoreUnderflow()
ignoreBetaRange()


def extrapolate(filename,order=2,Ncut=0,show_results=False,plot_results=False):

    data      = readTable(filename)
    Nts       = data[0][Ncut:] 
    betacs    = data[1][Ncut:]
    betacerrs = data[2][Ncut:]
    
    Tds, Tderrs, spacings = [], [], []
    r0 = r0_phys(year=2014, units="MeVinv")
    
    for i, Nt in enumerate(Nts):
    
        def getTd(beta):
            # We are only using latticeParams to get a, so Ns doesn't matter
            lp = latticeParams(Nsigma=1, Ntau=Nt, coupling=beta, scaleType='r0', scaleYear=2014, paramYear=2015)
            return lp.getT()*r0
    
        spacings.append(1/r0_div_a(betacs[i],year=2015))

        Td, Tderr = error_prop( getTd, [betacs[i]], [betacerrs[i]] ) 
        Tds.append(Td)
        Tderrs.append(Tderr)

    Tds, Tderrs, spacings = toNumpy(Tds, Tderrs, spacings)    
    
    # Perform continuum-limit extrapolation
    result, result_err, chidof = continuumExtrapolate( spacings, Tds, Tderrs, order=order, xtype="a",
                                                       show_results=show_results, plot_results=plot_results, )
    clearPlot()

    # The Bayesian Akaike information criterion tries to estimate P(M|D), the probability of your "model" 
    # (in this case the set of all analysis choices you made, including the continuum fit function form, 
    # whether you threw out any points, etc) given the data. The relationship is like P(M|D) ~ exp(-BAIC/2). 
    # Many ingredients go into the BAIC including the chi^2. If the fit is poor, P(M|D) goes down. If you 
    # cut data, P(M|D) goes down. If your fit has lots of parameters, P(M|D) goes down.
    IC = BAIC(xdata=spacings, ydata=Tds, cov=np.diag(Tderrs), func=_powerSeries, params=result, Ncut=Ncut)

    return result[0], result_err[0], IC


# Now what we're going to do is try a bunch of different continuum fits, since the situation is
# a bit ambiguous. We build up a list of results and corresponding errors and information criteria.
# The BAIC will give a systematic way of dealing with systematic uncertainty due to model choice.
def getContinuumEstimate(filename,Ncut_mid,Ncut_max):
    BAICs, Tdr0s, Tdr0es = [], [], []

    # We try some quadratic (in a^2) models
    for Ncut in range(Ncut_mid+1):
        logger.info('quadratic, Ncut =',Ncut)
        result, result_err, IC = extrapolate(filename,order=2,Ncut=Ncut)
        Tdr0s.append(result)
        Tdr0es.append(result_err)
        BAICs.append(IC)

    # and some linear (in a^2) ones
    for Ncut in range(Ncut_mid,Ncut_max+1):
        logger.info('linear, Ncut =',Ncut)
        result, result_err, IC = extrapolate(filename,order=1,Ncut=Ncut)
        Tdr0s.append(result)
        Tdr0es.append(result_err)
        BAICs.append(IC)

    BAICs, Tdr0s, Tdr0es = toNumpy(BAICs, Tdr0s, Tdr0es)

    # This is essentially a weighted average, where each continuum-limit result gets weighted by
    # the model probability discussed above, P(M|D). The error is computed as the weighted average
    # of the errors plus an additional term due to model spread. Loosely we think of that first
    # part as the "statistical error" and the additional term as the "systematic error".
    Tdr0, Tdr0e = modelAverage(Tdr0s, Tdr0es, BAICs)
    return Tdr0, Tdr0e



Tdr0_lit, Tdr0_lite = 0.7457, 0.0045

# First let's see what happens when we apply our Bayesian model averaging to the literature.
# How will what we find compare with what they find? 
Tdr0, Tdr0e = getContinuumEstimate('literature.d',Ncut_mid=3,Ncut_max=6)
logger.info('       lit-BMA =',get_err_str(Tdr0    ,Tdr0e))
logger.info('           lit =',get_err_str(Tdr0_lit,Tdr0_lite))
logger.info('q(lit vs. lit) =',round( gaudif(Tdr0,Tdr0e,Tdr0_lit,Tdr0_lite),4 ))

Tdr0_lit2, Tdr0_lite2 = Tdr0, Tdr0e 

# What about our RW data?
Tdr0, Tdr0e = getContinuumEstimate('RWbetas.txt',Ncut_mid=2,Ncut_max=4)
logger.info('               SRI =',get_err_str(Tdr0    ,Tdr0e))
logger.info('               lit =',get_err_str(Tdr0_lit,Tdr0_lite))
logger.info('    q(SRI vs. lit) =',round( gaudif(Tdr0,Tdr0e,Tdr0_lit,Tdr0_lite),4 ))
logger.info('q(SRI vs. lit-BMA) =',round( gaudif(Tdr0,Tdr0e,Tdr0_lit2,Tdr0_lite2),4 ))

# What about our spline data?
Tdr0, Tdr0e = getContinuumEstimate('cont_extrap_beta.d',Ncut_mid=2,Ncut_max=4)
logger.info('               SRI =',get_err_str(Tdr0    ,Tdr0e))
logger.info('               lit =',get_err_str(Tdr0_lit,Tdr0_lite))
logger.info('    q(SRI vs. lit) =',round( gaudif(Tdr0,Tdr0e,Tdr0_lit,Tdr0_lite),4 ))
logger.info('q(SRI vs. lit-BMA) =',round( gaudif(Tdr0,Tdr0e,Tdr0_lit2,Tdr0_lite2),4 ))

