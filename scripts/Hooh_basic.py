'''
Hooh_basic

retrying to make this HOOH script from scrap

Using Morris_et_al_2011 data about HOOH production from Hepes buffer to get a 'productuion rate of HOOH vi Hepes buffer'  

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_3_Morris_2011_ROS/scripts/
'''

from scipy.integrate import *
from scipy import *
from pylab import *
from numpy import * 
from matplotlib import *
import sys    #what does this do? 
import pandas as pd
import numpy as np      # done again to give it a nickname 
import matplotlib.pyplot as plt   #done again to specifically import objects from pyplot within Matplotlib





###############################

#  Importing the data frame

###############################

 
HOOH_df = pd.read_csv('../data/hooh_blank.txt', delimiter =',', header= None, names =( 'Time (1/days)','HOOH concentration'))




##############################

#   data arrays

###############################

HOOH_data = 10**array(HOOH_df.iloc[:,1])   #raised to 10 because ot the way data thief recorded the data and its scale

print(HOOH_data)





##############################

#   time arrays

##############################

#HOOH_times = 10**array(HOOH_df.iloc[:,0])

HOOH_times = array(HOOH_df.iloc[:,0])
print(HOOH_times)





##############################

#    graphing the data 

##############################


plt.scatter(HOOH_times, HOOH_data, marker = 'x', c = 'r', label = 'HOOH produced') 
plt.xlabel('Time (day $^{-1}$)') 
plt.ylabel('HOOH concentration')

#show()


Hs = HOOH_data
Ts = HOOH_times

Hsd = np.std(Hs)     #getting numpy to find the stdv of the data (but really will probs be 0.1 bc the data set doesn't have triplicates in it. 


#############################

#Euler's Integration

############################

HsModel = np.array([]) 
tsModel = np.array([])
step = 0.3 #delta t
ndays = 4
times = np.linspace(0,ndays,ndays/step)
Hepes = 3.75
delta = 0.8
kHOOH = 1.5
H = 0.0

for t in times: 
	dHdt = kHOOH*Hepes - delta*H
	H = H + dHdt*step
	HsModel = np.append(HsModel,H)
	
plt.plot(times,HsModel,c='r')

plt.show()

