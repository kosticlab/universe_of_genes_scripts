###cat together summary data
##20190221
#Tierney

import os
import pandas as pd

covDirs=[x for x in os.listdir('.') if 'coverage'==x[:8] and 'gz' not in x]

dfs=[]
for c in covDirs:
    iterationDirs=os.listdir(c)
    iterationDirs=[x for x in iterationDirs if 'png' not in x and 'csv' not in x and 'gz' not in x and '_' not in x]
    for d in iterationDirs:
        df=pd.read_csv(c+'/'+d+'/summary_data.csv')
        df.loc[:,'iteration']=[d]
        df.loc[:,'coverage_range']=[c]
        dfs.append(df)

output=pd.concat(dfs)
output.to_csv('summary_data_coverage_analysis.csv')
