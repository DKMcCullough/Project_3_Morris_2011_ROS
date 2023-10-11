'''
viz_Pro_help_and_light.py

visualise different buffer proiductions of H from Morris 2013 (miliQ_bufferH_2013MorrisSIfig1b) 
  HOOH production in 10 mM Buffers plus seawater media incubated  in lightexposed milli-Q water (i.e., without seawater solutes). rep1-3 are bio1 rep4-6 are bio2

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


df_all = pd.read_csv('../data/Pro_help_and_light_Morris2011_5.csv')


df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df_all = df_all.rename({'Time':'time'}, axis=1)    #'renaming column to make it callable by 'times'
df_all.fillna(0)



df_all['rep1'] = df_all['rep1'].astype('float')
df_all['rep2']  = df_all['rep2'].astype('float')
df_all['rep3']  = df_all['rep3'].astype('float')
df = df_all 

#logging data for latter graphing 
df['log1'] = np.log(df['rep1'])
df['log2'] = np.log(df['rep2'])
df['log3'] = np.log(df['rep3'])

#####

#raw avgs reps

df['abundance'] =  np.nanmean(np.r_[[df[i] for i in ['rep1','rep2','rep3']]],axis=0)
df['sigma'] = np.nanstd(np.r_[[df[i] for i in ['rep1','rep2','rep3']]],axis=0)


#log avgs and stdvs

df['log_abundance'] = np.nanmean(np.r_[[df[i] for i in ['log1','log2','log3']]],axis=0)
df['log_sigma'] =  np.nanstd(np.r_[[df[i] for i in ['log1','log2','log3']]],axis=0)


##############################

#   data arrays

###############################


orgs = df['Strain'].unique()
norgs = orgs.shape[0]


treats = df['Treatment'].unique()
ntreats = treats.shape[0]

lights = df['Light'].unique()
nlightss = lights.shape[0]
##############################

#    graphing the data 

##############################


fig1,(ax1)= plt.subplots(2,2, figsize = (11,8))
fig1.suptitle('Pro with different Help and Light', size = 22)
fig1.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
ax1[1,0].set_ylabel('Pro (cells/ml)')
ax1[1,0].set_xlabel('Time' )
ax1[0,0].semilogy()
ax1[1,1].set_ylabel('STDV')
ax1[1,1].set_xlabel('Mean' )

fig2,(ax2) = plt.subplots(2,2,figsize = (11,8))
fig2.suptitle(' Log Pro with different Help and Light',size = 22)
fig2.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
ax2[1,0].set_ylabel('Pro (cells/ml)')
ax2[1,0].set_xlabel('Time' )
ax2[0,0].semilogy()
ax2[1,1].set_ylabel('Log STDV')
ax2[1,1].set_xlabel('Log Mean' )

for s,ns in zip(orgs,range(norgs)):
    dfw = df_all[(df_all['Strain'] == s)]
    df_low =  dfw[(dfw['Light']== 'Low')]
    df_high = dfw[(dfw['Light']== 'High')]
    df_ax = dfw[(dfw['Treatment']== 'EZ55')]
    df_help = dfw[(dfw['Treatment']== 'EZ55')]
    print(dfw)
    ax1[0,0].errorbar(df_low.time,df_low.abundance, yerr=df_low.sigma, marker = 'd', c='yellow',label =  'low light raw mean')
    ax1[0,0].scatter(df_low.time,df_low.rep1, c='g', label = 'rep1')
    ax1[0,0].scatter(df_low.time,df_low.rep2, c='b', label = 'rep2')
    ax1[0,0].scatter(df_low.time,df_low.rep3, c='r',label = ' rep3')
    ax1[0,0].text(1.2,0.5,'Low Light',horizontalalignment='center', verticalalignment='center', transform=ax1[0,1].transAxes)
    ax1[0,0].semilogy()
    l1 = ax1[0].legend(loc = 'upper left')
    l1.draw_frame(False)
    ax1[0,1].scatter(df_low.abundance, df_low.sigma, c = 'purple')
    
    '''
    ax1[0,0].errorbar(df.time,df.abundance, yerr=df.sigma, marker= '.',markersize= 10, label =('Mean'), color = 'c' )
    ax1[0,0].errorbar(df.time,df.avg1, yerr=df.stdv1, marker= '.',markersize= 10, label =('avg1'), color = 'pink' )
    ax1[0,0].plot(df.time,df.avg2, yerr=df.stdv2, marker= '.',markersize= 10, label =('avg2'), color = 'purple' )
    #annotate side of large fig
    ax1[nt,0].text(1.2,0.5,'Buffer: '+ str(t),horizontalalignment='center', verticalalignment='center', transform=ax1[nt,1].transAxes)
    #log y axis and add legend to dynamics graph 
    ax1[nt,0].semilogy()
    l1 = ax1[nt,0].legend(loc = 'lower center')
    l1.draw_frame(False)
    #graph stats
    ax1[nt,1].scatter(df.abundance,df.sigma, c='c')
    ax1[nt,1].scatter(df.avg1,df.stdv1, c='pink')
    ax1[nt,1].scatter(df.avg2,df.stdv2, c='purple')
    #graph logged dynamics 
    ax2[nt,0].errorbar(df.time,df.log_abundance, yerr=df.sigma, marker= '.',markersize= 10, label =('Log Mean'), color = 'c' )
    ax2[nt,0].errorbar(df.time,df.lavg1, yerr=df.stdlog1, marker= '.',markersize= 10, label =('Log avg1'), color = 'pink' )
    ax2[nt,0].errorbar(df.time,df.lavg2, yerr=df.stdlog2, marker= '.',markersize= 10, label =('Log avg2'), color = 'purple' )
    #log y axis and add legend to dynamics graph 
    ax2[nt,0].semilogy()
    l2 = ax2[nt,0].legend(loc = 'upper right')
    l2.draw_frame(False)
    #graph loggedstats
    ax2[nt,1].scatter(df.log_abundance,df.log_sigma, c='c')
    ax2[nt,1].scatter(df.lavg1,df.stdlog1, c='pink')
    ax2[nt,1].scatter(df.lavg2,df.stdlog2, c='purple')
    #annotate sidf large figs
    ax2[nt,0].text(1.2,0.5,'Buffer: '+ str(t),horizontalalignment='center', verticalalignment='center', transform=ax1[nt,1].transAxes)

    #stdv vs mean
    

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.show()
fig1.savefig('../figures/Pro_help_with_light_raw.png')
fig2.savefig('../figures/logDynamics_Pro_help_with_light.png')
'''


print('\n ~~~****~~~****~~~ \n')
print('done with Pro light and help assays')
print('\n ~~~****~~~****~~~ \n')
