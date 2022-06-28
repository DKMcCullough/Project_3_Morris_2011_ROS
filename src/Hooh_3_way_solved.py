'''
Hooh_3_way_solved.py

Trying to match HOOH production data from Morris et al 2011 using an analytical solution, euler's aproximation, and ODEint 

Using Morris_et_al_2011 data about HOOH production from Hepes buffer to get a 'productuion rate of HOOH vi Hepes buffer'  

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_3_Morris_2011_ROS/scripts/
'''

from scipy import *
from scipy.integrate import *
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   




###############################

#  Importing the data frame

###############################

 
HOOH_df = pd.read_csv('../data/hooh_blank.txt', delimiter =',', header= None, names =( 'Time (1/days)','HOOH concentration'))




##############################

#   data arrays

###############################

HOOH_data = 10**np.array(HOOH_df.iloc[:,1])   #raised to 10 because ot the way data thief recorded the data and its scale
#HOOHs_nm = 1000**np.array(HOOH_data)   #TRIED TO GET IN nM to match better with other things 
#print(HOOH_data)





##############################

#   time arrays

##############################


HOOH_times =np.array(HOOH_df.iloc[:,0])
#print(HOOH_times)





##############################

#    graphing the data 

##############################


plt.scatter(HOOH_times, HOOH_data, marker = 's', s=50, c = 'r', label = 'Measured HOOH') 
plt.xlabel('Time (days)', fontsize=16) 
plt.ylabel('HOOH concentration (\u03BCM)',fontsize=16)

plt.yscale('log')
plt.tick_params(labelsize=12)
#plt.show()


Hs = (HOOH_data)
Ts = HOOH_times

Hsd = np.std(Hs)     #getting numpy to find the stdv of the data (but really will probs be 0.1 bc the data set doesn't have triplicates in it. 




####################################

#analytical solution

####################################

#initial values and creating time array

delta = 0.5
S_HOOH = 2.3
step = 0.05 #delta t
ndays = 7
times = np.linspace(0,ndays,int(ndays/step))

def f(t, S_HOOH, delta):
	H = (S_HOOH/delta)*(1-e**(-delta*t))
	return H

Hs = f(times,S_HOOH,delta)
#print(times,Hs) 


plt.plot(times,Hs,c='g',linestyle = '-.',label='Analytical Solution')
#, marker='*'


#############################

#Euler's Integration

############################

HsEuler = np.array([]) 
H = 0.0
t0 = 0

for t in times: 
	HsEuler = np.append(HsEuler,H)
	dHdt = S_HOOH - delta*H
	H = H + dHdt*step
	
#plt.plot(times,HsEuler,c='red',label = 'Model')#,label = "Euler's Aproximation")

plt.legend()

#plt.show()



####################################

#ODE int

####################################


from scipy.integrate import odeint

def HsODEint(H,t):
	dHdt = S_HOOH-delta*H
	#print(dHdt,t,H)
	return dHdt


ode_solutions = odeint(HsODEint,0,times)


plt.plot(times,ode_solutions,c='purple', linestyle = ':', label = 'Odeint Approximation')
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)


plt.legend(loc = 'lower right')
plt.show()

	




	
	
