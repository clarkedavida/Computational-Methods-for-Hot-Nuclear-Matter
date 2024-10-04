# 
# litCompareByBeta.py                                                               
#
# Compare RW with spline and literature. We go Nt by Nt and do a Z-test.
#  
# D. Clarke
# 

import numpy as np
from latqcdtools.base.readWrite import readTable
from latqcdtools.base.cleanData import intersectAtCol
from latqcdtools.base.printErrorBars import get_err_str
from latqcdtools.statistics.statistics import gaudif
import latqcdtools.base.logger as logger
from latqcdtools.base.utilities import printClean

lit    = readTable('literature.d')
sriRW  = readTable('RWbetas.txt')
sriSPL = readTable('cont_extrap_beta.d')

def compare(A,B):

    _A, _B = intersectAtCol(A,B,col=0)
    
    for i in range(len(_A[0])):
        Nt     = int(_A[0][i])
        b_lit  = _B[1][i] 
        be_lit = _B[2][i] 
        b_sri  = _A[1][i] 
        be_sri = _A[2][i]
        dBeta  = np.abs(b_lit-b_sri)/b_lit
    
        printClean(Nt, get_err_str(b_sri,be_sri), get_err_str(b_lit,be_lit), gaudif(b_lit,be_lit,b_sri,be_sri), dBeta )


logger.info('lit vs SRI-RW:')
compare(lit,sriRW)
logger.info()
logger.info('lit vs SRI-SPL')
compare(lit,sriSPL)
logger.info()
logger.info('SRI-RW vs SRI-SPL')
compare(sriRW,sriSPL)