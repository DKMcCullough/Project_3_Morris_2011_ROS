'''
Hooh_3_way_solved.py

Trying to match HOOH production data from Morris et al 2011 using an analytical solution, euler's aproximation, and ODEint 

*******	then using error sum of squares to find the best line to fit the data ****** 

Using Morris_et_al_2011 data about HOOH production from Hepes buffer to get a 'productuion rate of HOOH vi Hepes buffer'  

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_3_Morris_2011_ROS/scripts/
'''

from scipy import *
from scipy.integrate import odeint
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   




##############################

#   data 

###############################

#  Importing the data frame

trial_df = pd.read_csv('../data/hooh_blank.txt', delimiter =',', header= None, names =( 'Time (1/days)','HOOH concentration'))
detox_df = pd.read_csv('../data/hooh_ez55.txt', delimiter =',', header= None, names =( 'Time (1/days)','HOOH concentration'))

#   slicing df to get data arrays

trial_data = 10**np.array(trial_df.iloc[:,1])   #raised to 10 because ot the way data thief recorded the data and its scale
detox_data = 10**np.array(detox_df.iloc[:,1])   #raised to 10 because ot the way data thief recorded the data and its scale

#   slicing df to get time arrays

trial_times = np.array(trial_df.iloc[:,0])
detox_times = np.array(detox_df.iloc[:,0])

#    graphing the data 


plt.scatter(trial_times, trial_data, marker = 's',s=35,c = 'r', label = 'abiotic HOOH') 
plt.scatter(detox_times, detox_data, marker = 's',s=35,c = 'g', label = 'detoxed HOOH')

plt.xlabel('Time (days)', fontsize = 20) 
plt.ylabel('HOOH concentration ($\u03BC$M)', fontsize = 20)
plt.suptitle('HOOH dynamics', fontsize = 25)
#plt.legend()


####################################

#   Model
####################################

#global model values set up 
h_convert = 0.65                      #term used to convert input Hepes concentration to HOOH felt by cells
S_HOOH = (3.75)*h_convert       #3.75 micromolar Hepes used for this actually; K_h ac calculated in morris and Zinser 2013 is between 10 and 100 nM
#HOOH in micromolar on graph and HEPES input in micromolar as well...but maybe per L? vs what HOOH is? Not sure.
step = 0.3 #delta t
ndays = 7.2
times = np.linspace(0,ndays,int(ndays/step))

#abiotic model

delta = 0.47
H = np.array([])

def HsODEint(H,t):
    dHdt = S_HOOH-delta*H
    #print(dHdt,t,H)
    return dHdt

ambient_solutions = odeint(HsODEint,0,times)


plt.plot(times,ambient_solutions,c='red', linewidth = 2, label = 'abiotic model')

#detoxed model  (additive d) 
ddetox = 8.2
H = np.array([])

def HsODEint(H,t):
    dHdt = S_HOOH-delta*H-ddetox*H
    #print(dHdt,t,H)
    return dHdt

detox_solutions = odeint(HsODEint,0,times)

plt.plot(times,detox_solutions, c='green', linewidth = 2, label = 'detoxed model')


#multiplicative d model
ddetox = 21
H = np.array([])

def HsODEint(H,t):
    dHdt = S_HOOH-delta*ddetox*H
    #print(dHdt,t,H)
    return dHdt

detox_solutions = odeint(HsODEint,0,times)

plt.plot(times,detox_solutions, c='purple', linewidth = 2, label = 'm_d detoxed model')



plt.legend(loc = 'center right')


plt.show()

print('Done')

