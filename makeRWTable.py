# 
# makeRWTable.py                                                               
# 
# C. Rohde 
# 

import numpy as np
import glob,sys
import latqcdtools.base.logger as logger
from latqcdtools.base.readWrite import writeTable

betac=sys.argv[1]

outFiles = list(glob.iglob(f'*{betac}*'))

betals=[]
Pmls=[]
Pels=[]
Ntls=[]
Nsls=[]
Nconfls=[]
Xmls=[]
Xels=[]
Bmls=[]
Bels=[]
Nt=None
Ns=None

RePls=[]
ImPls=[]
Bls=[]

for filename in outFiles:
    try:
        logger.info('opening',filename)
        inFile=open(filename,'r')
    except:
        logger.TBError('error opening file',filename)
   
    #isolate polyakov loop absolute value
    for line in inFile:
        if not line.startswith('#'):
            continue
    
        col=line.split()
    
        if col[-3]=='loop':
            #absPls.append((eval(col[-1])[0]**2+eval(col[-1])[1]**2)**0.5)
            RePls.append(eval(col[-1])[0])
            ImPls.append(eval(col[-1])[1])
        
        if col[-3]=='Plaquette':
            Bls.append(float(eval(col[-1])))

        if col[-3]=='beta':
            beta=float(col[-1])
        
        if len(col)>6:
            if col[-6]=='Lattice':
                Nt=col[-1]
                Ntls.append(Nt)
        
                Ns=col[-2]
                Nsls.append(Ns)

                if len(set(Ntls))!=1:
                    logger.info('ERROR','multiple Nt values found')
                    sys.exit(-1)

                if len(set(Nsls))!=1:
                    logger.info('ERROR','multiple Ns values found')
                    sys.exit(-1)

RePls=np.array(RePls)
ImPls=np.array(ImPls)
Bls=np.array(Bls)
    
writeTable('Nt'+Nt+'_Ns'+Ns+'b'+str(beta).replace('.','')+'.txt',\
        np.sqrt(RePls**2+ImPls**2),Bls,header=['|P|','[]'])
logger.info('wrote table:','Nt'+Nt+'_Ns'+Ns+'b'+str(beta).replace('.','')+'.txt')