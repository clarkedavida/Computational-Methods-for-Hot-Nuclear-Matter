# 
# litCompareByBeta.py                                                               
# 
# D. Clarke
# 

from latqcdtools.base.readWrite import readTable
from latqcdtools.base.cleanData import intersectAtCol
from latqcdtools.base.printErrorBars import get_err_str
from latqcdtools.statistics.statistics import gaudif
import latqcdtools.base.logger as logger


lit = readTable('literature.d')
sri = readTable('RWbetas.txt')
#sri = readTable('cont_extrap_beta.d')

lit, sri = intersectAtCol(lit,sri,col=0)

logger.info('Nt  SRI    LITERATURE    q')
for i in range(len(lit[0])):
    Nt     = int(lit[0][i])
    b_lit  = lit[1][i] 
    be_lit = lit[2][i] 
    b_sri  = sri[1][i] 
    be_sri = sri[2][i]

    logger.info(Nt, get_err_str(b_sri,be_sri), get_err_str(b_lit,be_lit), round(gaudif(b_lit,be_lit,b_sri,be_sri),4) )