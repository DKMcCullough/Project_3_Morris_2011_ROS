'''
viz_ocean_H_insitu.py


created by DKM

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_3_Morris_2011_ROS/scripts/
'''

from scipy import *
from scipy.integrate import *
import pandas as pd
import numpy as np      
import matplotlib.pyplot as plt   
#from Basemap import * #asemap not working yet 



###############################

#  Importing the data frame

###############################

 #UH18301 Pro in 3.75 mM HEPES or Taps buffer. High light (24uC in a Sunbox -  noon maximum of about 250 quanta m) 


df_all = pd.read_csv('../data/ocean_H_morris2011f2a.csv')


df_all.drop(df_all.columns[df_all.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)



df = df_all 


##############################

#   data arrays

###############################

lats = df['Latitude'].unique()
deps = df['Depth'].unique()

df.plot(kind="scatter", x = 'Depth', y='HOOH', alpha = 1, title = 'in situ HOOH')

df.plot(kind="scatter", x='Latitude', y='Depth', alpha=(abs((df.HOOH)/np.max(df.HOOH))))
#df.plot(kind="scatter", x = (abs((df.HOOH)/np.max(df.HOOH))), y='Depth', alpha = 1)

##############################

#    graphing the data 

##############################
'''

fig1,(ax1)= plt.subplots(ntreats,2, figsize = (11,8))
fig1.suptitle('Raw H dynamics in the Ocean', size = 22)
fig1.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
ax1[1,0].set_ylabel('HOOH concentration (\u03BCM)')
ax1[3,0].set_xlabel('Time' )
ax1[0,0].semilogy()
ax1[1,1].set_ylabel('STDV')
ax1[3,1].set_xlabel('Mean' )

fig2,(ax2) = plt.subplots(ntreats,2,figsize = (11,8))
fig2.suptitle(' Log  H dynamcis in the ocean',size = 22)
fig2.subplots_adjust(right=0.85, left=0.10,wspace = 0.25, hspace = 0.30)
ax2[1,0].set_ylabel('HOOH concentration (\u03BCM)')
ax2[2,0].set_xlabel('Time' )
ax2[0,0].semilogy()
ax2[1,1].set_ylabel('Log STDV')
ax2[3,1].set_xlabel('Log Mean' )

for t,nt in zip(treats,range(ntreats)):
    df = df_all[(df_all['Buffer'] == t)]
    ax1[nt,0].errorbar(df.time,df.abundance, yerr=df.sigma, marker= '.',markersize= 10, label =('Mean'), color = 'c' )
    ax1[nt,0].errorbar(df.time,df.avg1, yerr=df.stdv1, marker= '.',markersize= 10, label =('avg1'), color = 'pink' )
    ax1[nt,0].errorbar(df.time,df.avg2, yerr=df.stdv2, marker= '.',markersize= 10, label =('avg2'), color = 'purple' )
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
fig1.savefig('../figures/ocean_H_insitu_dynamics.png')
fig2.savefig('../figures/ocean_H_insitu_stats.png')


'''
print('\n ~~~****~~~****~~~ \n')
print('done with singular hepes')
print('\n ~~~****~~~****~~~ \n')
