#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 23:50:19 2022

name:Pro_mono_morris_2011.py

location: /Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/ROS_focused/Project_6_cocultures/src

author: DKM

to do: Get model working for the treatments of both species

"""

#Pro UH18301

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import odeint



#SynWH7803 (Vol 28 and 52) is KatG possitive
#SynWH8102 (Vol 54) is KatG negative 


#data importmulti_H_Pros.csv"

df_all = pd.read_csv("../data/multi_H_Pros.csv")
df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all.dropna()
#df_all['avg_exp'] = df_all['avg_exp'].fillna(value = 0.0) #filling Nans with 0.0 in 'avg' column 
df_all = df_all.rename({'Time (days)':'times'}, axis=1)    #'renaming column to make it callable by 'times'

#making log of data to look at error
df_all['log1'] = np.log(df_all['rep1'])
df_all['log2'] = np.log(df_all['rep2'])
df_all['log3'] = np.log(df_all['rep3'])

df_all['abundance'] =  np.nanmean(np.r_[[df_all[i] for i in ['rep1','rep2','rep3']]],axis=0)
df_all['log_abundance'] = np.nanmean(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)
df_all['log_sigma'] = np.std(np.r_[[df_all[i] for i in ['log1','log2','log3']]],axis=0)

####################################

# Slicing Data

####################################

strains = df_all['Strain'].unique()

nstrains = strains.shape[0]

colors = ('green', 'blue','c', 'orange', 'r', 'k') 
markers = ('s','v','o','*','d','v')

for s,ns in zip(strains,range(nstrains)): 
    df_s = df_all[(df_all['Strain'] == s)]
    treats  = df_s['Treatment(HOOH (uM))'].unique()
    ntreats = treats.shape[0]
    fig1,(ax1)= plt.subplots( figsize = (10,6))
    for t,nt in zip(treats,range(ntreats)): 
        count = nt
        df = df_s[(df_s['Treatment(HOOH (uM))'] == t)]
        times = (df['times'])
        pdata = (df['abundance'])
        ax1.plot(times, pdata, marker= markers[count], markersize= 10, label =(str(t)+' nM HOOH'), color = colors[count] ) 
    fig1.suptitle('Prochlorococcus Monoculture Dynamics '+ str(s))
    ax1.set_ylabel('Pro cells (per ml)')
    fig1.supxlabel('Time (days)')
    l1 = ax1.legend(loc = 'center right', prop={"size":13}) 
    l1.draw_frame(False)#print(df)
    ax1.semilogy()
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.show()
    fig1.savefig('../figures/Pro_'+str(s)+'_data.png')
    
    #inits = pd.read_csv(("../data/inits/"+s+"_inits.csv")) 
    
    
    
    
inits = pd.read_csv(("../data/inits/hepes"+s+".csv")) 
#inits = pd.read_csv(("../data/inits/hepes"+s+".csv")) 
'''




#
##################################################3
# parameter and variable Set UP 
#############################################


step = 0.001
ndays = 2
mtimes = np.linspace(0,ndays,int(ndays/step))


Qnp = 1#(9.4e-15*(1/(14.0))*1e+9)  #Nitrogen Quota for Pro from Bertilison 


k1p =  0.00002     #Pro alpha
k2 = 0.88    #************
ksp = k2/k1p
dp = 0.2   #pro delta
kdam = 0.0035  #hooh mediated damage rate of Pro  
deltah = 0.002       #decay rate of HOOH via Syn 
phip = 2.5e-7    #0007  #detoxification-based decay of HOOH via Syn in this case
rho =  0.002
SN = 0
Sh = 8


#empty arrays to be populated by odeint when calling leak function
P = np.array([])
N = np.array([])
H = np.array([])


#initial values to be used for odeint start 
P0 = df_4['avg_exp'].iloc[0]
N0 = 1.0e8        #nM 
H0 = 400    #400 actual.....need to get data h0    #nM
#inits = (P0,D0,N0,H0)

#H0 is 400 or 0, SH is 0 both times???
params = (ksp,k2,dp,deltah,phip,SN,Sh,kdam)
y = [P,N,H]
inits = (P0, N0, H0)

def pro_mono(y,t,params):
    ksp,k2,dp,deltah,phip,SN,Sh,kdam = params[0], params[1], params[2], params[3],params[4], params[5], params[6], params[7]
    P,N,H = y[0],y[1],y[2]
    dPdt = (k2 * N /( (ksp) + N) )*P - (dp *P) - kdam*H*P     
    dNdt =  SN - ((k2 * N /( (ksp) + N) )*P* Qnp) - rho*N    
    dHdt = Sh - deltah*H -phip*H*P #phi being S cell-specific detox rate
    return [dPdt,dNdt,dHdt]

Pro = odeint(pro_mono,inits, mtimes, args = (params,))
Ps = Pro[:,0]
Ns = Pro[:,1]
Hs = Pro[:,2]

#Model graphing 


ax1.plot(mtimes, Ps , linewidth = 3, color = 'g', label = 'model P ') #'k1 =' + str(k1p))
ax1.text(0.02,1.5e4,s= ('kdam = '+str(kdam)+ ' and phi = '+str(phip)))

ax1.semilogy()

ax1.legend(loc = 'best')

plt.show()


print('*** Done ***')


treatments = [50,200,400,800,10000]  #HOOH nM []   #1000 or 10000 nM HOOH for last treatment....kdam no longer saturation if later case
treatments_smol = treatments[1:3]
dc = dict()

####slicing df into treatments and saving into dictionary (dc)######
for i in treatments:
    df_i = df_all[df_all["treatment"].isin([i])]
    name = ('df_' + str(i))
    dc.update({name : df_i})  #update dictionary of dfs with each loop itteration. 
    #print(dc)

########################

#graphing

#########################
plt.rcParams['font.size'] = 14

fig1,ax1 = plt.subplots()

colors = ('green', 'c', 'orange', 'r', 'k') #make into a set in loop? for c in colors, color = count(c)?????
markers = ('s','v','o','*','d')
#for df_i in dc:
for t in treatments: 
    count = treatments.index(t)
    #print(count)
    df = dc['df_'+str(t)]
    times = df['times']
    data = df['avg_exp'] #data was loggeed in original graph; transformed in excel before read in
    ax1.plot(times, data, linestyle = 'None', marker= markers[count], markersize= 10, label = (str(t) +' nM HOOH'), color = colors[count])  #color = colors(i))
    #ax1.plot(times,data,linestyle='-', linewidth=0.25, color='black', marker = 'None')
#plt.show()
ax1.set(xlabel= 'Time (days)', ylabel='Biomass ( cells  ml$^{-1}$)')
ax1.set_title = ('Prochlorococcus U18301 w/ HOOH') #not showing up for some reason? 
ax1.legend(loc = 'lower left', prop={"size":10})
ax1.semilogy()
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)


############################

#model

#############################


#dPdt = (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration) (Population size) 
#dSdt = (Supply of nutriet) - (max growth rate)(nutient concentration)/michelis-mention constant) + (nutrinet concentration)
#Qn = (9.4e-15*(1/14.0)*1e+9)  #Nitrogen Quota for Pro 

step = 0.01
ndays = 2.1
times = np.linspace(0,ndays,int(ndays/step))
Qn = (9.4e-15*(1/(14.0))*1e+9)   #Nitrogen Quota for Pro from Bertillison? 

#S_base =   30.0             #3.0 ~160nM from BCC paper calculations 
     
P = (5e5) # units are cells per mL
#S = (S_base + 0)    #nM N per ml for units      #    0.164 micromolar rediual N from Calfee_et_al 2022

#k1= 0.002   #ALPHA  
#k2 = 0.22  #VMAX
#nrdelta = 0.06      #nutrient replete delta  #Kdam for this 
#nddelta = 0.19       #nutrient deplete delta kddam for this
#mumax = k2    #set max/min for this known from lit? 

#parameter values (k1 = alpha, k2 = Vmax, kdam = initial HOOH damage, kddam = HOOH damage when N runs out)
k1s = np.r_[[0.02, 0.02, 0.02, 0.02, 0.02]]*1e+100     
k2s = [0.32, 0.32, 0.32, 0.32, 0.32]
kdams = [0.07,0.2, 0.6, 1.3, 3.1]
kddams = [0.03, 0.5, 2, 2, 1.4]

params = list(zip(k1s,k2s,kdams,kddams))

fig2,ax2 = plt.subplots()
#df = dc['df_50']
#data = df['avg_exp']
#S_base =  (data.max()- data.iloc[0])*Qn 

for t in treatments: 
    count = treatments.index(t)
    k1 = params[count][0]
    k2 = params[count][1]
    kdam = params[count][2]
    kddam = params[count][3]
    ks = (k2/k1)   #set max value for this that is know from lit? (between 0.01 and 0.015 for N metabolism in )
    SsEuler = np.array([])
    PsEuler = np.array([])
    P = 1.3e5
    #print(kdam,kddam)
    S_base =   162.0             #3.0 ~160nM from BCC paper calculations 
    df = dc['df_50']
    data = df['avg_exp']
    S =  (data.max()- data.iloc[0])*Qn    #using QN and 50 treatment as base N for test.     
    # units are cells per mL
    #S = S_base    #nM N per ml for units      #    0.164 micromolar rediual N from Calfee_et_al 2022
    for t in times:
            PsEuler = np.append(PsEuler,P)
            SsEuler = np.append(SsEuler,S)
            #if (S>2e-3):
            #    delta = kdam
            #else:
            #    delta = kddam
            delta = kdam
            #print(S)
            #print(P)
            dPdt = k2 * P * S /((ks) + S) - delta*P
            dSdt = -P*(k2*(S)/((ks)+S))*Qn
            if S+dSdt*step <0:                    #making sure S isnt taken negative and therefore causing issues when we log transform the data
                    S = S + dSdt*step  #S = 4e-47
            else:
                    S = S + dSdt*step
            P = P + dPdt*step
    ax1.plot(times,(PsEuler), linestyle = 'dashed', color = colors[count]) 
    ax2.plot(times,(SsEuler), linestyle = 'dashed', color = colors[count])


    
#ax2_legend = [zip(treatments,colors)]  #ipz doesn't give a printable object for legend to show
ax2.set(xlabel= 'Time (days)', ylabel='Nitrogen (nM)')
ax2.set_title = ('Nutrient dynamics') #not showing up for some reason? 
#ax2.legend([(list(a2_legend))] ,loc = 'lower left')
ax2.semilogy()

fig3,ax3 = plt.subplots()

colors = ('green', 'c', 'orange', 'r', 'k') #make into a set in loop? for c in colors, color = count(c)?????
markers = ('s','v','o','*','d')
for x, y, c, m, in zip(treatments,kdams, colors, markers):
    ax3.plot(x, y,linestyle = 'none', color = c, marker = m, markersize = 10,label = ('kdam = '+str(y)))
    #ax3.plot(x, y,linestyle='-', linewidth=0.25, color='black', marker = 'None'
#ax3.plot(kdams,treatments,linestyle = 'none', marker=markers[kdam[]],color=colors[count])
#ax3.set(xlabel= 'HOOH treatment (nM)', ylabel='kdam (day$^{-1}$)')
#ax3.plot(treatments,kdams,linewidth = 0.25, c = 'k')
#ax3.semilogy()
ax3.legend(loc = 'lower right')
plt.show()


'''

print('*** Done ***')