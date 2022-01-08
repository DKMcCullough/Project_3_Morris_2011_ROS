'''
Practice of finding the solution for and ODE and handcoding Eulers method for estimation

/Users/dkm/Documents/Talmy_research/Basics/scripts/mu_P.py

by DKM

'''

from matplotlib import *
from scipy import *
from math import *
import numpy as np
import matplotlib.pyplot as plt


########################

# initial conditions

########################

step = 0.5
ndays = 30
times = np.linspace(0,ndays,int(ndays/step))    #use int() function to make the ndays/step and integer output not a float
mu = 0.2
Pi = 1.0*10**4    #p initial


##########################

# function

##########################

def f(t, Pi):
        P =  Pi*e**(mu*t)
        return P



#########################

#finding solutions for P

########################

Ps = f(times,Pi)
#print(Ps)


#######################

#graphing calculated solutions

#######################

plt.plot(times,Ps,c='r',marker='*',label="Analytical solution") 
plt.xlabel('Times')
plt.ylabel('Phytoplankton Biomass')

#plt.show()



###########################

#Euler's integration 

###########################

P = 1.0*10**4  
t0 = 0
PsModel = r_[[]]

for t in times:
	dPdt = mu*P
	P = P + dPdt*step
	PsModel = np.append(PsModel,P)


#############################

#Graphing model solutions

#############################


plt.plot(times,PsModel,c='b',marker='.',linestyle='--',label="Euler's Approximation")
#plt.xlabel('Model Times')       #commented out to get a single set of axis labels on teh singular graph
#plt.ylabel('Modeled Phytoplankton Biomass')

#plt.show()



#############################

# using OdeInt

############################

from scipy.integrate import odeint

P = 1.0*10**4

def  diffD(P, t): 	#differential fuction
	dPdt = mu*P
	print(dPdt,P)
	return dPdt

solutions = odeint(diffD,[P],times)
#print(solutions)
#print(times)

##############################

#graphing OdeInt solutiions

####################################

plt.plot(times,solutions,c = 'g',marker = '.',linestyle = ':', label = 'ODEint Approximation')

plt.legend()

plt.show()




