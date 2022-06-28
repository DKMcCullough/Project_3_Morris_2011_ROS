'''

name:   Pro_modeled.py 

location: /Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_3_Morris_2011_ROS/src

author: DKM


goal: model  Pro data from HOOH and no HOOH trials

'''



import pandas as pd
import numpy as np
#from matplotlib import *
import matplotlib.pyplot as plt
from scipy import *
from scipy.integrate import odeint
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

plt.scatter(x = pro_hooh['times'], y = pro_hooh['biomass'],marker = 's', s = 50, c = 'b', label = 'HOOH trial')

plt.scatter(x = pro_detoxed['times'], y = pro_detoxed['biomass'], marker = 's', s = 50,color = 'orange', label = 'HOOH detoxed')

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

#   Model detoxed   (via EZ55)

#######################################
#initial values 

#P0 = pro_detoxed['biomass'][0]  #P0 set to initial data point 
P0 = pro_detoxed['biomass'][1]  #P0 set to second data point b/c of dip in abundance btwn t0 and t30
N0 = 1.22e5
inits = (P0,N0)

#parameters
k2 =  0.88    #mumax    #this is high.........0.7 is kind upper limit shown in literature? 
k1 =  0.00002     #alpha
d= 0.0000000000001
params = [k1,k2,d]

#time window 
step = 0.01 #delta 
ndays = 7
mtimes = np.linspace(0,ndays,int(ndays/step))

#empty P and N arrays 
P = np.array([])
N = np.array([])
y = [P,N]



#function set up for ode int


def Pdetox(y,t,params):
    k1,k2,d = params[0], params[1], params[2]
    P,N = y[0],y[1]
    dPdt = k2 * P * N /( (k2/k1) + N) -d*P
    dNdt =-P*( k2*N)/((k2/k1)+N)
    return [dPdt,dNdt]

#solve ODEs via odeint
detox = odeint(Pdetox, inits, mtimes,args = (params,))

#redefine where P and N are in returned matrix from ode int
Ps = detox[:,0]
Ns = detox[:,1]

#plot P
plt.plot(mtimes, Ps, linestyle = ':', linewidth = 3, color = 'orange', label = 'd Model')


#######################################

#   Model NOT detoxed

#######################################
P0 = pro_hooh['biomass'][0]
inits = (P0,N0)

HOOH = 3.75e6    #Hepes buffer of 3.75 micromolar put in   #HOOH in nanomolar? 
kdam = (HOOH)*0.051
params = [k1,k2,kdam]


#function set up for ode int


def Ptrial(y,t,params):
    k1,k2,kdam = params[0], params[1], params[2]
    P,N = y[0],y[1]
    dPdt = k2 * P * N /( (k2/k1) + N) - kdam+P
    dNdt =-P*( k2*N)/((k2/k1)+N)
    return [dPdt,dNdt]


#solve ODEs via odeint
trial = odeint(Ptrial, inits, mtimes,args = (params,))

#redefine where P and N are in returned matrix from ode int
Ps = trial[:,0]
Ns = trial[:,1]

plt.plot(mtimes, Ps, linestyle = ':', linewidth = 3, color = 'b', label = 'kdam model')

#plt.legend(loc = 'center right')
plt.semilogy()
plt.show()

print("Done")

