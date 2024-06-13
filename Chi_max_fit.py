import numpy as np
import matplotlib.pyplot as plt
from latqcdtools.statistics.fitting import Fitter
from latqcdtools.base.initialize import initialize, finalize
import latqcdtools.base.logger as logger

# Initialize logging
initialize('fitting_linear.log')

Nss = np.array([36, 48, 60])  # Spatial extents
max_X = np.array([0.602900546, 1.19202461, 2.15226992])
max_X_err = np.array([0.00596052962, 0.015628063, 0.060702249])  # Example max_X values corresponding to each Ns
Ns_cubed = (Nss**3)  # Compute Ns^3 for the fitting

# Define the linear fit function
def fit_func(x, params):
    A = params[0]
    B = params[1]
    return A + B * x


fitter = Fitter(fit_func, Ns_cubed, max_X, max_X_err)

# Perform the fit
res, res_err, chi_dof, stats = fitter.try_fit(start_params=[1, 0], algorithms=['curve_fit'], detailedInfo=True)

# Log the results
logger.info(" A, B : ", res)
logger.info(" Ae, Be: ", res_err)
logger.info("chi2/d.o.f.: ", chi_dof)
logger.info("     logGBF: ", stats['logGBF'])
logger.info("       pcov: \n\n", stats['pcov'], "\n")

# Plot the fit and data
fitter.plot_fit(domain=(min(Ns_cubed), max(Ns_cubed)))
fitter.plot_data()
plt.xlabel('$N_s^3$')
plt.ylabel('$\\chi_{max}$')
plt.title('Linear Fit of $\\chi_{max}$ vs $N_s^3$')
plt.xticks(rotation=30)
plt.savefig('Chi_max_fit.png')

# Finalize and close any resources
finalize()
