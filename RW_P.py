
import numpy as np
from latqcdtools.base.readWrite import readTable
from latqcdtools.physics.statisticalPhysics import reweight
from latqcdtools.physics.lattice_params import latticeParams
import latqcdtools.base.logger as logger
from latqcdtools.statistics.jackknife import jackknife
from latqcdtools.statistics.bootstr import bootstr
from latqcdtools.base.plotting import set_params,plt,plot_dots,set_default_param
from latqcdtools.base.printErrorBars import get_err_str

# Input the N_tau you wanna look at
Nt=16
Ns=3*Nt

# Here are our critical beta values from literature, expressed as characters and
# a floating point number.
oldBeta = {
  6  : 5.89425,
  8  : 6.06239,
#  8  : 6.026,
  10 : 6.20873,
  12 : 6.33514,
  14 : 6.4473,
  16 : 6.5457,
  18 : 6.6331,
  20 : 6.7132,
}
oldBetac = {
  6  : '589425',
  8  : '606239',
#  8  : '6026',
  10 : '620873',
  12 : '633514',
  14 : '64473',
  16 : '65457',
  18 : '66331',
  20 : '67132',
}

# This is the starting point around which we reweight.
beta0 = oldBeta[Nt]  

# INPUT: We need a table of Polyakov loop measurements and corresponding average plaquette
#        on every configuration for a given beta value.
PL0, act0 = readTable('Nt'+str(Nt)+'/Nt'+str(Nt)+'_Ns'+str(Ns)+'b'+oldBetac[Nt]+'.txt')

V = 6*Ns**3*Nt # number of plaquettes (b/c we printed average plaquette) 
S = act0*V     # the original action. Remember the idea of reweighting is that we rescale our 
               #   measurements according to exp(- change in S as you change beta to the RW point).
               #   this is equivalent to rescaling (reweighting) the measurement by the likelihood
               #   of having a configuration with action SRW (as opposed to the original action) 

def RWP(data,xRW,x0) -> float:
    """ Reweight the Polyakov loop

    Args:
        data (list): an list [polyakov loop, action] 
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
    return reweight(X**2,xRW,x0,S) - reweight(X,xRW,x0,S)**2


betas   = []
MRW     = []
MRW_err = []
SUSCRW     = []
SUSCRW_err = []
newBetas = np.linspace(beta0-0.010,beta0+0.010,31)
for beta in newBetas:
    pRW        = beta 
    MRWm, MRWe = jackknife( RWP, [PL0, S], args=(pRW,beta0) )
    SUSCRWm, SUSCRWe = jackknife( RWSUSC, [PL0, S], args=(pRW,beta0) )
    MRW.append(MRWm)
    MRW_err.append(MRWe)
    SUSCRW.append(SUSCRWm)
    SUSCRW_err.append(SUSCRWe)
    betas.append(beta)


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
    RWBetas = np.linspace((1-0.001)*oldBeta[Nt],(1+0.001)*oldBeta[Nt],1001)
    SUSCmax = 0
    betamax = -1
    for beta in RWBetas:
        pRW     = beta 
        SUSCRWm = RWSUSC(PL0S,pRW,beta0)
        if SUSCRWm > SUSCmax:
            SUSCmax = SUSCRWm
            betamax = beta
    return betamax

def findTc(PL0S) -> float:
    """ Same thing as above but for the temperature.

    Args:
        PL0S (array): original Polyakov loop

    Returns:
        float: critical temperature 
    """
    RWBetas = np.linspace((1-0.001)*oldBeta[Nt],(1+0.001)*oldBeta[Nt],1001)
    SUSCmax = 0
    betamax = -1
    for beta in RWBetas:
        pRW     = beta 
        SUSCRWm = RWSUSC(PL0S,pRW,beta0)
        if SUSCRWm > SUSCmax:
            SUSCmax = SUSCRWm
            betamax = beta
    lp = latticeParams(Ns,Nt,betamax,None,None,None,scaleType='r0')
    return lp.getT() 

bmax, bmaxerr = bootstr( findBetaMax, [PL0,S], numb_samples=300 ) 
Tc, Tcerr     = bootstr( findTc, [PL0,S], numb_samples=300 ) 
logger.info('beta_c =',get_err_str(bmax,bmaxerr))
logger.info('T_c =',get_err_str(Tc,Tcerr))

set_default_param(font_size=8)
plot_dots(betas,SUSCRW,SUSCRW_err,color='blue',marker=None)
set_params(xlabel='$\\beta$',ylabel='$\\chi_{|P|}$',title='RW Pure SU(3), $N_\\tau='+str(Nt)+'$')
plt.tight_layout() # this repositions things to make the plot look nicer
plt.show() 
