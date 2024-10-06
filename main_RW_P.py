
import sys
import numpy as np
from latqcdtools.base.readWrite import readTable
from latqcdtools.physics.statisticalPhysics import reweight
import latqcdtools.base.logger as logger
from latqcdtools.statistics.jackknife import jackknife
from latqcdtools.base.plotting import set_params,plt,plot_dots,set_default_param,\
    BACKGROUND,plot_vspan
from latqcdtools.base.printErrorBars import get_err_str
from latqcdtools.base.utilities import find_nearest_idx

Nt=int(sys.argv[1])

if Nt==12:
    Ns=2*Nt
else:
    Ns=3*Nt

logger.info('Nt=',Nt)


NBINS = 20
AXIS  = 1
NRW   = 301

# Here are our critical beta values from literature, expressed as characters and
# a floating point number.
oldBeta = {
  6  : 5.89425,
  8  : 6.06239,  # Need a point between the given two.
#  8  : 6.05,
  10 : 6.20873,
  12 : 6.33514,  # Clean for Ns=24
  14 : 6.4473,   
  16 : 6.5457,
  18 : 6.6331,
  20 : 6.7132,
}
oldBetac = {
  6  : '589425',
  8  : '606239',
#  8  : '605',
  10 : '620873',
  12 : '633514',
  14 : '64473',
  16 : '65457',
  18 : '66331',
  20 : '67132',
}
deltaRW = {
  6  : 0.002,
  8  : 0.001,
  10 : 0.0003,
  12 : 0.001,
  14 : 0.0003,
  16 : 0.0003,
  18 : 0.0002,
  20 : 0.0003
}

# This is the starting point around which we reweight.
beta0 = oldBeta[Nt]  

# INPUT: We need a table of Polyakov loop measurements and corresponding average plaquette
#        on every configuration for a given beta value. You should be able to replace this
#        path with the corresponding table for whatever configuration you want to analyze.
PL0, act0 = readTable('Nt'+str(Nt)+'/Nt'+str(Nt)+'_Ns'+str(Ns)+'b'+oldBetac[Nt]+'.txt')

V = 6*Ns**3*Nt # number of plaquettes (b/c we printed average plaquette) 
S = act0*V     # the original action. Remember the idea of reweighting is that we rescale our 
               #   measurements according to exp(- change in S as you change beta to the RW point).
               #   this is equivalent to rescaling (reweighting) the measurement by the likelihood
               #   of having a configuration with action SRW (as opposed to the original action) 


def RWP(data,xRW,x0) -> float:
    """ Reweight the Polyakov loop

    Args:
        data (list): a list [polyakov loop, action] 
        xRW (float): the point we are RWing to 
        x0 (float): the starting point

    Returns:
        float: reweighted Polyakov loop 
    """
    X = data[0]
    S = data[1]
    return reweight(X,xRW,x0,S) 


def RWSUSC(data,xRW,x0) -> float:
    """ Reweight the susceptibility. The susceptibility is an observable that is
    defined in terms of expectation values. At the same time, we think of the
    reweight() method as a redefined expectation value.

    Args:
        data (list): an list [polyakov loop, action] 
        xRW (float): the point we are RWing to 
        x0 (float): the starting point

    Returns:
        float: reweighted susceptibility 
    """
    X = data[0]
    S = data[1]
    return Ns**3*( reweight(X**2,xRW,x0,S) - reweight(X,xRW,x0,S)**2 )


SUSCRW     = []
SUSCRW_err = []
RWBetas = np.linspace((1-deltaRW[Nt])*beta0,(1+deltaRW[Nt])*beta0,NRW)


# findBetaMax has trouble resolving the mean under the jackknife if the number of RW points is
# too small. But increasing the number of points massively slows things down. I decided to use
# this loop to get a more accurate picture of the mean, then jackknife(findBetaMax) only to
# estimate the error bar
chimax=-1
bmax=-1
for pRW in RWBetas:
    SUSCRWm, SUSCRWe = jackknife( RWSUSC, [PL0, S], args=(pRW,beta0), numb_blocks=NBINS, conf_axis=AXIS )
    if SUSCRWm > chimax:
        chimax = SUSCRWm
        bmax = pRW 
    SUSCRW.append(SUSCRWm)
    SUSCRW_err.append(SUSCRWe)


def findBetaMax(PL0S) -> float:
    """ In the infinite volume limit, at the critical beta value, the system undergoes
    a phase transition. The order parameter (Polyakov loop) as a function of the temperature
    experiences a jump discontinuity, and correspondingly the Polyakov loop susceptibility, 
    which can be thought of as a temperature derivative of the Polyakov loop, spikes to
    infinity. (This is the statement that the slope of the Polyakov loop as a function of
    the temperature becomes infinite.)

    At finite volume, there are no phase transitions. (One way to think about this is the 
    Lee-Yang theorem.) The jump discontinuity softens to a rapid change in slope whose
    inflection point we use to estimate a "pseudo-critical" temperature. (As opposed to a
    critical temperature, which is what we would have in the infinite volue limit. We don't
    usually call something pseudo-critical unless it reduces to true, critical behavior in
    some appropriate limit.) The inflection point corresponds to maximum slope, and hence
    the susceptibility will exhibit a peak there. Therefore we estimate the location of
    the "psuedo-critical" point by the location of the peak in susceptibility.

    This routine finds the beta value that maximizes the susceptibilty, which we call
    beta_c (although we may more aptly call it beta_pc; I along with many other physicists
    are often sloppy in this way). This corresponds to the temperature that maximizes the
    susceptibility since
        T = 1 / a(beta) Nt 
    The way we find the maximum is by reweighting to several beta values within a tenth of
    a percent the original. (The range in beta is small because the lattice spacing a is
    a complicated function of beta that turns out to be extremely sensitive to it.) This is
    why we needed the reweighting in the first place.

    Args:
        PL0S (array): original Polyakov loop

    Returns:
        float: critical beta 
    """
    SUSCmax = -1
    betamax = -1
    for pRW in RWBetas:
        SUSCRWm = RWSUSC(PL0S,pRW,beta0)
        if SUSCRWm > SUSCmax:
            SUSCmax = SUSCRWm
            betamax = pRW 
    return betamax


# Find critical beta and corresponding error in susceptibility
_, bmaxerr = jackknife( findBetaMax, [PL0,S], numb_blocks=NBINS, conf_axis=AXIS) 
logger.info('beta_c(JACK) =',get_err_str(bmax,bmaxerr))
maxsusc,maxsuscerr = jackknife( RWSUSC, [PL0, S], args=(bmax,beta0), numb_blocks=NBINS,conf_axis=AXIS )
logger.info('chi_max =',get_err_str(maxsusc,maxsuscerr))

# Read in literature results
litNt, litbetac, litbetaerr = readTable('literature.d')
iNt = find_nearest_idx(litNt,Nt)

# Plotting...
set_default_param(font_size=8)
plot_dots(RWBetas,SUSCRW,SUSCRW_err,color='green',marker=None,alpha=0.2,ZOD=BACKGROUND+1) # RW points
plot_dots([bmax],[maxsusc],[maxsuscerr],xedata=[bmaxerr],color='red',marker=None,
    label='$\\beta_c(N_\\tau)$') # Estimate for max
plot_vspan(minVal=litbetac[iNt]-litbetaerr[iNt], maxVal=litbetac[iNt]+litbetaerr[iNt],color='yellow',
    alpha=0.3, ZOD=BACKGROUND, label='literature') # literature
set_params(xlabel='$\\beta$',ylabel='$\\chi_{|P|}$',title='RW Pure SU(3), $N_\\tau='+str(Nt)+'$')
plt.tight_layout() # this repositions things to make the plot look nicer
plt.savefig('RW_figures/Nt'+str(Nt)+'.pdf')
plt.show() 
