# 
# continuumExtrapolate.py                                                               
# 
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.readWrite import readTable
from latqcdtools.math.num_deriv import diff_deriv
from latqcdtools.math.spline import getSpline
from latqcdtools.statistics.statistics import gaudif, error_prop
from latqcdtools.statistics.bootstr import bootstr_from_gauss
from latqcdtools.physics.continuumExtrap import continuumExtrapolate
from latqcdtools.physics.constants import r0_phys
from latqcdtools.physics.lattice_params import latticeParams
from latqcdtools.physics.referenceScales import r0_div_a
from latqcdtools.base.printErrorBars import get_err_str
from latqcdtools.base.plotting import clearPlot, latexify, set_params

latexify()

def extrapolate(filename):

    Nts, betacs, betacerrs = readTable(filename)
    
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
        logger.info(Nt,get_err_str(Td,Tderr))
    
    
    # Perform O(a^4) continuum-limit extrapolation
    result, result_err, chidof = continuumExtrapolate( spacings, Tds, Tderrs, order=2, xtype="a",
                                                       show_results=True, plot_results=True )
    
    # Do a Z-test against literature result,  
    Tdr0, Tdr0e = result[0], result_err[0]
    Tdr0_lit, Tdr0_lite = 0.7457, 0.0045
    logger.info('lit =',get_err_str(Tdr0_lit,Tdr0_lite))
    logger.info('q(ours vs. lit) =',gaudif(Tdr0,Tdr0e,Tdr0_lit,Tdr0_lite))
    set_params(title=filename) 
    clearPlot()

#extrapolate('RWbetas.txt')
#extrapolate('cont_extrap_beta.d')
extrapolate('literature.d')