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
df_all = df_all.fillna(0)



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

for s,ns in zip(orgs,range(norgs)):
    for t,nt in zip(treats,range(ntreats)):
        dfw = df_all[(df_all['Strain'] == s) & (df_all['Treatment']==t)]
        df_low =  dfw[(dfw['Light']== 'Low')]
        df_high = dfw[(dfw['Light']== 'High')]
    #df_ax = dfw[~(dfw['Treatment']== 'EZ55')]
    #df_help = dfw[(dfw['Treatment']== 'EZ55')]
        #print(dfw)

        fig1,(ax1)= plt.subplots(2,2, figsize = (11,8))
        fig1.suptitle(str(s) + ' ' + str(t), size = 22)
        fig1.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
        ax1[0,0].set_ylabel('Pro (cells/ml)')
        ax1[0,0].set_xlabel('Time' )
        ax1[0,0].semilogy()
        ax1[1,0].set_ylabel('Pro (cells/ml)')
        ax1[1,0].set_xlabel('Time' )
        ax1[1,0].semilogy()
        ax1[1,1].set_ylabel('STDV')
        ax1[1,1].set_xlabel('Mean' )
        ax1[0,1].set_ylabel('STDV')
        ax1[0,1].set_xlabel('Mean' )

        fig2,(ax2) = plt.subplots(2,2,figsize = (11,8))
        fig2.suptitle(' Log ' +  str(s) + ' ' + str(t),size = 22)
        fig2.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
        ax2[0,0].set_ylabel('Pro (cells/ml)')
        ax2[0,0].set_xlabel('Time' )
        ax2[0,0].semilogy()
        ax2[1,0].set_ylabel('Pro (cells/ml)')
        ax2[1,0].set_xlabel('Time' )
        ax2[1,0].semilogy()
        ax2[1,1].set_ylabel('Log STDV')
        ax2[1,1].set_xlabel('Log Mean' )
        ax2[0,1].set_ylabel('Log STDV')
        ax2[0,1].set_xlabel('Log Mean' )

        ax1[0,0].errorbar(df_low.time,df_low.abundance, yerr=df_low.sigma, marker = 'd', c='brown',label =  'raw mean')
        ax1[0,0].scatter(df_low.time,df_low.rep1, c='c', label = 'rep1')
        ax1[0,0].scatter(df_low.time,df_low.rep2, c='b', label = 'rep2')
        ax1[0,0].scatter(df_low.time,df_low.rep3, c='purple',label = ' rep3')
        ax1[0,0].text(1.2,0.5,'Low Light',horizontalalignment='center', verticalalignment='center', transform=ax1[0,1].transAxes)
        l1 = ax1[0,0].legend(loc = 'lower right')
        l1.draw_frame(False)
        ax1[0,1].scatter(df_low.abundance, df_low.sigma, c = 'purple')
    
        ax1[1,0].errorbar(df_high.time,df_high.abundance, yerr=df_high.sigma, marker = 'd', c='brown',label =  'raw mean')
        ax1[1,0].scatter(df_high.time,df_high.rep1, c='c', label = 'rep1')
        ax1[1,0].scatter(df_high.time,df_high.rep2, c='b', label = 'rep2')
        ax1[1,0].scatter(df_high.time,df_high.rep3, c='purple',label = ' rep3')
        ax1[1,1].text(1.2,0.5,'High Light',horizontalalignment='center', verticalalignment='center', transform=ax1[1,1].transAxes)
        #l1 = ax1[1,0].legend(loc = 'upper left')
        #l1.draw_frame(False)
        ax1[1,1].scatter(df_high.abundance,df_high.sigma, c = 'purple')
        
        #logged graph
        ax2[0,0].errorbar(df_low.time,df_low.log_abundance, yerr=df_low.log_sigma, marker = 'd', c='brown',label =  'log mean')
        ax2[0,0].scatter(df_low.time,df_low.log1, c='c', label = 'log1')
        ax2[0,0].scatter(df_low.time,df_low.log2, c='b', label = 'log2')
        ax2[0,0].scatter(df_low.time,df_low.log3, c='purple',label = ' log3')
        ax2[0,0].text(1.2,0.5,'Low Light',horizontalalignment='center', verticalalignment='center', transform=ax2[0,1].transAxes)
        l2 = ax2[0,0].legend(loc = 'lower right')
        l2.draw_frame(False)
        ax2[0,1].scatter(df_low.log_abundance, df_low.log_sigma, c = 'purple')
    
        ax2[1,0].errorbar(df_high.time,df_high.log_abundance, yerr=df_high.log_sigma, marker = 'd', c='brown',label =  'log mean')
        ax2[1,0].scatter(df_high.time,df_high.log1, c='c', label = 'log1')
        ax2[1,0].scatter(df_high.time,df_high.log2, c='b', label = 'log2')
        ax2[1,0].scatter(df_high.time,df_high.log3, c='purple',label = ' log3')
        ax2[1,1].text(1.2,0.5,'High Light',horizontalalignment='center', verticalalignment='center', transform=ax2[1,1].transAxes)
        #l3 = ax2[1,0].legend(loc = 'upper left')
        #l3.draw_frame(False)
        ax2[1,1].scatter(df_high.log_abundance,df_high.log_sigma, c = 'purple')
        
        fig1.savefig('../figures/'+str(s)+'raw'+str(t)+'.png')
        fig2.savefig('../figures/'+str(s)+'log'+str(t)+'.png')


    

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.show()
#fig1.savefig('../figures/Pro_help_with_light_raw.png')
#fig2.savefig('../figures/logDynamics_Pro_help_with_light.png')



print('\n ~~~****~~~****~~~ \n')
print('done with Pro light and help assays')
print('\n ~~~****~~~****~~~ \n')
