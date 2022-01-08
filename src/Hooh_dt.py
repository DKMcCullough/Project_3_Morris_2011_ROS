'''
Hooh_dt

comparison of how the HOOH 

Using Morris_et_al_2011 data 

created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_2_MOrris_2011_ROS/scripts/
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

#  Importing the data

###############################

 
HOOH = pd.read_csv('../data/hooh_blank.txt', delimiter =',', header= None, names =( 'Time (1/days)','HOOH concentration'))

##############################

#   data arrays

###############################

HOOH_data = 10**array(blank.iloc[:,1])   #raised to 10 because ot the way data thief recorded the data and its scale

#print(HOOH_data)




##############################

#   time arrays 

##############################

HOOH_times = array(blank.iloc[:,0])

#print(HOOH_times)





################################

# graphing the data

################################




f1, (ax1,ax2) = plt.subplots(1,2)
ax1.scatter(HOOH_times,HOOH_data,marker='x',c='r', label = 'HooH blank')
ax1.set_xlabel('Time (day $^{-1}$)')

#need to set y axes to be the same thing....maybe I should just graph them on the same axis? 

show()
 
Hs = HOOH_data
Ts = HOOH_times    #or should I use an average of the blank and ez55 times? 



############################

# STDV

############################

#what to do for STDV???

Hsd = np.std(Hs)



################################

#  Integration

###############################


# function to be integrated

def f(u,y,kdetox,ez55):
        H = u[0]
        dHst = Shooh*(Hepes) - delta*(H)
        return concatenate([r_[[dHst]]])

'''

# calling ODEint 
   
def integrate (params,inits,htimes,ptimes,forshow=False,delt=1.0 / 24.0):
        ndays = 10
        mtimes = linespace(0,ndays,into(ndays/delt)) 
        kdetox = exp(params[0])
        u = odeint(f,inits,mtimes,args=kdetox)
        if forshow==False: 
                Hinds = r_[[where(abs(a-mtimes)==min(abs(a-mtimes)))[0][0] for a in blank_times]]
                Pinds = r_[[where(abs(a-mtimes)==min(abs(a-mtimes)))[0][0] for a in pro_times]]
                Hnt = u.T[0][Hinds]
                Pnt = u.T[1][Pinds]
        else: 
                Hnt = u.T[0]
                Pnt = u.T[1]
        
        return Hnt, Pnt



###################################

# generic arrays and optimization parameters

###################################



stds = zeros(3) + .05 # this controls how large the random steps are in the parameter search (see below)

opt = r_[[1]] # this allows you to control whether each parameter is imposed or fittedi

names = ['kdetox'] # name each parameter array - useful for printing the output


nits = 1000 # number of iterations

pits = 100  # frequency with which to print results to the user about the progress of the algorithm

burnin = 100 # burnin period - don't 'store' the results for the first number of iterations as these are likely to be very far from the optimal ssolutions


   '''




  
