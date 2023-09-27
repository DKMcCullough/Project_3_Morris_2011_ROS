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

 #UH18301 Pro in 3.75 mM HEPES or Taps buffer. High light (24uC in a Sunbox -  noon maximum of about 250 quanta m) 


df_all = pd.read_csv('../data/Buffers_Morris_2011_f1.csv')

df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time (days)':'time'}, axis=1)    #'renaming column to make it callable by 'times'

#creating log and stats for data 
df_all['log1'] = np.log(df_all['rep1'])
df_all['log2'] = np.log(df_all['rep2'])
df_all['log3'] = np.log(df_all['rep3'])

df_all['abundance'] =  np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['log_abundance'] = np.nanmean(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)
df_all['log_sigma'] = np.std(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)


##############################

#   data arrays

###############################


orgs = df_all['organism'].unique()

assays = df_all['ID'].unique()
nassays = assays.shape[0]


#HOOHs_nm = 1000**np.array(HOOH_data)   #TRIED TO GET IN nM to match better with other things 
#print(HOOH_data)
colors = ('r','g')
markers = ('o','*')

##############################

#    graphing the data 

##############################
for a,n in zip(assays,range(nassays)):
    df = df_all[(df_all['ID'] == a)]
    treats  = df['Treatment'].unique()
    ntreats = treats.shape[0]
    fig1,(ax1)= plt.subplots(2,1, figsize = (12,8))
    for t,nt in zip(treats,range(ntreats)): 
        count = nt
        print(t)
        df = df[(df['ID']==a) & (df['Treatment']==t)]
        pdf = df[df['organism'] == 'P']
        hdf = df[df['organism'] == 'H']
        print(pdf)
        ax1[0].plot(pdf['time'], pdf['abundance'], marker= markers[count], markersize= 10, label =(str(t)+' produced HOOH'), color = colors[count] ) 
        ax1[1].plot(hdf['time'], hdf['abundance'], marker= markers[count], markersize= 10, label =(str(t)+' produced HOOH'), color = colors[count] ) 
    fig1.suptitle('Abiotic Dynamics '+ str(a))
    ax1[0].set_ylabel('HOOH concentration (\u03BCM)')
    ax1[1].set_ylabel('Pro abundance (cell/ml)')
    fig1.supxlabel('Time (days)')
    l1 = ax1[0].legend(loc = 'center right', prop={"size":13}) 
    l1.draw_frame(False)#print(df)
    ax1[0].semilogy()
    ax1[1].semilogy()
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.show()
    fig1.savefig('../figures/Hproduction_'+str(a)+'_data.png')



    
'''
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

	

'''
print('\n ~~~****~~~****~~~ \n')
print('done with singular hepes')
print('\n ~~~****~~~****~~~ \n')
