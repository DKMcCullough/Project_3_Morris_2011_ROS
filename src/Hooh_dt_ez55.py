'''
Hooh_vs_hepes

comparison of how the HOOH and proclorococcus curves look based on if the HOOH is added in all at once in the begining of the expiriment or if Hepes buffer is used. 

Using Morris_et_al_2011 data 

created by DKm

location:/Users/dkm/Documents/Talmy_research/Zinser_and_Ben/Project_2_MOrris_2011_ROS/scripts/
'''

from scipy.integrate import *
from scipy import *
from pylab import *
from numpy import * 
from matplotlib import *
import sys    #what does this do? 
import pandas as pd
import numpy as np #do I need to import this twice or can I do it as above and just add np part? 
import matplotlib.pyplot as plt

###############################

#  Importing the data

###############################

 
pro_hooh = pd.read_csv('../data/hooh_blank.txt', delimiter =',', header= None, names =( 'Time','UH1301 abundance'))
pro_ez55_hooh = pd.read_csv('../data/hooh_ez55.txt', delimiter =',', header= None, names =( 'Time','UH1301 abundance'))
pro_hepes = pd.read_csv('../data/UH18301_HEPES.txt', delimiter =',', header= None, names =( 'Time','UH1301 abundance'))
pro_ez55_hepes = pd.read_csv('../data/UH18301_HEPES_EZ55.txt', delimiter =',', header= None, names =( 'Time','UH1301 abundance'))


##############################

#   data arrays

###############################

hooh_data = 10**array(pro_hooh.iloc[:,1])
ez_hooh_data = 10**array(pro_ez55_hooh.iloc[:,1])
hepes_data = 10**array(pro_hepes.iloc[:,1])
ez_hepes_data = 10**array(pro_ez55_hepes.iloc[:,1])

print(hooh_data,'\n',pro_hooh,'\n',ez_hooh_data,'\n',pro_ez55_hooh, '\n', hepes_data,'\n', pro_hepes, '\n',  ez_hepes_data, '\n', pro_ez55_hepes)




##############################

#   time arrays 

##############################

hooh_times = array(pro_hooh.iloc[:,0])
ez_hooh_times = array(pro_ez55_hooh.iloc[:,0])
hepes_times = array(pro_hepes.iloc[:,0])
ez_hepes_times = array(pro_ez55_hepes.iloc[:,0])

print(hooh_times,'\n',ez_hooh_times,'\n', hepes_times, '\n',  ez_hepes_times,)


#these time arrays aare in log space...so must exponentiate each item in the arrays to have time in (1/day)
'''

hooh_t = numpy.exp(hooh_times)
ez_hooh_t = numpy.exp(ez_hooh_times)
hepes_t = numpy.exp(hepes_times)
ez_hepes_t = numpy.exp(ez_hepes_times)

print(hooh_t,ez_hooh_t, hepes_t, ez_hooh_t)

'''


################################

# graphing the data

################################




f1, (ax1,ax2) = plt.subplots(1,2)
ax1.scatter(hooh_times,hooh_data,marker='x',c='r', label = '0.8 \u03BCM Hooh and Prochlorococcus')
ax1.set_xlabel('Time (day $^{-1}$)')
ax1.set_ylabel('Prochlorococcus cell density (?)')
ax2.scatter(ez_hooh_times,ez_hooh_data,marker='x',c='b',label = '0.8 \u03BCM Hooh and Prochlorococcus and EZ55')
ax2.set_xlabel('Time (day $^{-1}$)')
ax2.set_ylabel('Prochlorococcus cell density (?)')



f2, (ax3,ax4) = plt.subplots(1,2)
ax3.scatter(hepes_t,hepes_data,marker='x',c='y', label = 'Hepes and Prochlorococcus')
ax3.set_xlabel('Time (day $^{-1}$)')
ax3.set_ylabel('Prochlorococcus cell density (?)')
ax4.scatter(ez_hepes_t,ez_hepes_data,marker='x',c='g',label = 'Hepes and Prochlorococcus and EZ55')
ax4.set_xlabel('Time (day $^{-1}$)')
ax4.set_ylabel('Prochlorococcus cell density (?)')



show ()



#something is wrong with the times...idk what though...









