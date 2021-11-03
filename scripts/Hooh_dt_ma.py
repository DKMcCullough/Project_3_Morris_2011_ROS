'''
Hooh_dt_ma

modeling HOOH curve from Morris data using Metropolis Algorithum

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

HOOH_data = 10**array(HOOH.iloc[:,1])   #raised to 10 because ot the way data thief recorded the data and its scale

#print(HOOH_data)




##############################

#   time arrays 

##############################

HOOH_times = array(HOOH.iloc[:,0])

#print(HOOH_times)





################################

# graphing the data

################################




f1, (ax1) = plt.subplots()
ax1.scatter(HOOH_times,HOOH_data,marker='x',c='r', label = 'HooH blank')
ax1.set_xlabel('Time (day $^{-1}$)')
ax1.set_ylabel('HooH concentration (\u03BCM)')

#need to set y axes to be the same thing....maybe I should just graph them on the same axis? 

show()

 
Hs = HOOH_data
Ts = HOOH_times    #or should I use an average of the blank and ez55 times? 



############################

# STDV

############################

#what to do for STDV???

Hsd = 0.1     #don't have stdv of this data so just using 0.1 for it. 



################################

#  Integration

###############################


# function to be integrated

def f(u,t,Shooh,delta):
        H = u[0]
        dHdt = Shooh*(Hepes) - delta*(H)      #Do I need an H value in the production portion of the equation? 
        print (H, dHdt)
        return concatenate([r_[[dHdt]]])


# calling ODEint 
   
def integrate (params,inits,Ts,forshow=False,delt=1.0 / 24.0):
        ndays = 10.0
        mtimes = linspace(0,ndays,int(ndays/delt)) 
        Shooh,delta = exp(params[0]),exp(params[1])
        u = odeint(f,inits,mtimes,args=(Shooh,delta))
        if forshow==False: 
                Hinds = r_[[where(abs(a-mtimes)==min(abs(a-mtimes)))[0][0] for a in Ts]]
                Hnt = u.T[0][Hinds]
        else: 
                Hnt = u.T[0]
        
        return Hnt



###################################

# generic arrays and optimization parameters

###################################



stds = zeros(3) + .05 # this controls how large the random steps are in the parameter search (see below)

opt = r_[[1]] # this allows you to control whether each parameter is imposed or fittedi

names = ['Shooh,delta'] # name each parameter array - useful for printing the output


nits = 1000 # number of iterations

pits = 100  # frequency with which to print results to the user about the progress of the algorithm

burnin = 100 # burnin period - don't 'store' the results for the first number of iterations as these are likely to be very far from the optimal ssolutions



################################
 
# Model Fitting 

################################




#first guesses for parameters

Shooh = 0.200  #just a guess
Hepes = 3.75*1000    #3.75mM Hepes uesed in Morris 2011 #*0.001 to get to microM
delta = 0.100   #just a guess 


#puttting in arrays for easier manipulation

params = r_[[Shooh, delta]]

params = log(params) 

npars = params.shape[0]     #number of parameters being searched through


#initial conditions

#inits = r_[[(Hs,Hepes)]]   #this made the entire data array of Hs the initial value
inits = r_[[0]]      #making H = 0 the initial value for the model.


m_times = linspace(0,10,int(10*24))

#do I need mtimes down here too?


# first run of model

Hnt = integrate(params,inits,m_times)

#chi = sum((Hnt - Hs)**2/(Hsd**2))      #Hnt and Hs are vastly different sizes so this won't run...but how to make the sizes match or have it work without the sizes matching? 

print('Hnt = ' + str(Hnt)) #has lots of values
print('Hs = ' + str(Hs))
print('Hsd = ' + str(Hsd)) #only one value (0.1)  






##################################

#graphing the model

#################################


#print(Hnt)
#print(m_times)



f2, (ax2) = plt.subplots() 
ax2.scatter(HOOH_times,HOOH_data,marker='x',c='r', label = 'HOOH data')
#ax3 = ax2.twin()
ax2.plot(m_times,Hnt,c='y', label = 'HOOH model')
ax2.set_xlabel('Time (day $^{-1}$)')
ax2.set_ylabel('HOOH concentration (\u03BCM)')
#ax3.set_xlabel('Time (day $^{-1}$')


show() 



'''
# distribution arrarys and acceptance ratios 

ar = 0.0     #acceptance ratio

ars = r_[[]]

Shoohs,deltas = r_[[]], r_[[]]

pall = [Shoohs,deltas] 





#run the fitting

for it in arange(1,nits,1):    
        parsnew = params + opt*normal(0,stds,npars)
        Shooh = 0.3
        delta = 0.1
        inits = r_[[Shooh,delta]]
        Hnt = integrate(parsnew,inits,mtimes)
        chinew = sum((

'''


