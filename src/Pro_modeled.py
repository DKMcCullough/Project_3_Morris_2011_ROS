'''

name:   Pro_modeled.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_3_Morris_2011_ROS/src

author: DKM


goal: model  Pro data from HOOH and no HOOH trials

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

#plt.scatter(x = pro_hooh['times'], y = pro_hooh['biomass'], label = 'HOOH trial')

#plt.scatter(x = pro_detoxed['times'], y = pro_detoxed['biomass'], label = 'HOOH detoxed (EZ55)')

plt.scatter(x = pro_hooh['times'], y = pro_hooh['biomass'], marker = 'D',  s = 50,label = 'HOOH trial')
plt.plot(pro_hooh['times'],pro_hooh['biomass'], color = 'k', linestyle='solid', linewidth=0.25)

plt.scatter(x = pro_detoxed['times'], y = pro_detoxed['biomass'], marker = 's', s = 50, label = 'HOOH detoxed (EZ55)')
plt.plot(pro_detoxed['times'], pro_detoxed['biomass'], linestyle='solid', linewidth=0.25)


plt.semilogy()

plt.legend(loc='lower left',prop={'size': 10}, fontsize=22)
plt.title('UH18301 Pro with HOOH', fontsize = '22')
plt.xlabel('Time (days)',fontsize = '18')
#plt.legend(prop={"size":14})
plt.ylabel('Cell Abundance (ml$^{-1}$)',fontsize = '18')
plt.xticks(fontsize = 14) 
plt.yticks(fontsize = 14)



#plt.show()



#######################################

#   Model detoxed 

#######################################

P = 1e5
N = 1e5
k2 = mumax = 0.5
k1 = alpha = 0.09
step = 0.1 #delta 
ndays = 7
mtimes = np.linspace(0,ndays,int(ndays/step))
Ps = np.array([]) 
Ns = np.array([])


for t in mtimes:
    Ps = np.append(Ps,P)
    Ns = np.append(Ns,N)
    dPdt = k2 * P * N /( (k2/k1) + N) 
    dNdt =-P*( k2*N)/((k2/k1)+N)
    if N+dNdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
        N = 0.00000000000000000004
    else:
        N = N + dNdt*step 
    P = P + dPdt*step

plt.plot(mtimes,  Ps, marker = '_', color = 'k',  label = 'Model')

plt.legend()



#######################################

#   Model NOT detoxed

#######################################
P = 1e5
N = 1e5
k2 = mumax = 0.88
k1 = alpha = 0.2
step = 0.1 #delta 
ndays = 7
mtimes = np.linspace(0,ndays,int(ndays/step))
Ps = np.array([]) 
Ns = np.array([])

#HOOH_df = pd.read_csv('../data/hooh_blank.txt', delimiter =',', header= None, names =( 'Time (1/days)','HOOH concentration'))
#HOOH = 10**np.array(HOOH_df.iloc[:,1]
HOOH = 2e6
kdam = (HOOH)*0.089


for t in mtimes:
    Ps = np.append(Ps,P)
    Ns = np.append(Ns,N)
    dPdt = k2 * P * N /( (k2/k1) + N) - kdam+P
    dNdt =-P*( k2*N)/((k2/k1)+N)
    if N+dNdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
        N = 0.00000000000000000004
    else:
        N = N + dNdt*step 
    P = P + dPdt*step

plt.plot(mtimes,  Ps, marker = '_', color = 'k'  )#label = 'Pro modeled')

plt.legend(loc = 'right middle')





plt.show()




print("Done")






