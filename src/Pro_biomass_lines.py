'''

name:   Pro_biomass_lines.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_3_Morris_2011_ROS/src

author: DKM


goal: import and visualise Pro data from HOOH and no HOOH trials      ---- lines 

'''



import pandas as pd
import numpy as np
from matplotlib import *
import matplotlib.pyplot as plt
from scipy.integrate import *
from scipy import *
from pylab import *



############################

#  Data Import from csv   

############################

pro_hooh =  pd.read_csv('../data/UH18301_HEPES.txt', sep=",", header=None)
pro_hooh.columns = ["times", "biomass"]
#pro_hooh['exp_value'] = np.exp(pro_hooh['biomass'])
pro_hooh['biomass'] = 10**pro_hooh['biomass']   #10** bc of way data theif stored values 

pro_detoxed = pd.read_csv('../data/UH18301_HEPES_EZ55.txt', sep=",", header=None)
pro_detoxed.columns = ["times", "biomass"]
#pro_detoxed['exp_value'] = np.exp(pro_detoxed['biomass'])
pro_detoxed['biomass'] = 10**pro_detoxed['biomass']

#print(pro_hooh)
#print(pro_detoxed)





######################################

#  Graphing Data                 

#####################################

plt.figure()           #graphed all treatments in same rep column each time.

plt.scatter(x = pro_hooh['times'], y = pro_hooh['biomass'], marker = 'D',  s = 50,label = 'HOOH trial')
plt.plot(pro_hooh['times'],pro_hooh['biomass'], color = 'k', linestyle='solid', linewidth=0.25)

plt.scatter(x = pro_detoxed['times'], y = pro_detoxed['biomass'], marker = 's', s = 50, label = 'HOOH detoxed (EZ55)')
plt.plot(pro_detoxed['times'], pro_detoxed['biomass'], linestyle='solid', linewidth=0.25)


plt.semilogy()

plt.legend(loc='center right',prop={'size': 10}, fontsize=22)
plt.title('UH18301 Pro with HOOH', fontsize = '22')
plt.xlabel('Time (days)',fontsize = '18')
#plt.legend(prop={"size":14})
plt.ylabel('Cell Abundance (ml$^{-1}$)',fontsize = '18')
plt.xticks(fontsize = 14) 
plt.yticks(fontsize = 14)



plt.show()




print("Done")

