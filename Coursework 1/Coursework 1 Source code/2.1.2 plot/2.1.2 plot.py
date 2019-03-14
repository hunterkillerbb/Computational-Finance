import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt

df = pd.read_csv('Data_value.csv')
a = df.values

plt.figure(figsize = (16,9))
plt.title('Figure 1: Portfolio value at time = 0 with S0 = X1',fontsize = 20 )
plt.xlabel('Path number N',fontsize =12)
plt.ylabel('Portvolio values',fontsize =12)
plt.plot(a[:,0],a[:,1],label = 'S0 = 60000',lw = 0.8, color = 'red',marker = 'x')
plt.hlines(-660.801,xmin = -1000, xmax = 101000,lw = 0.6, linestyles = '--',label ="L1 = -660.801",color ='black')
plt.legend( prop ={'size' : 15 })


plt.figure(figsize = (16,9))
plt.title('Figure 2: Portfolio value at time = 0 with S0 = X2',fontsize = 20)
plt.xlabel('Path number N',fontsize = 12)
plt.ylabel('Portvolio values',fontsize =12)
plt.plot(a[:,0],a[:,2],label = 'S0 = 105000',lw = 0.8, color = 'blue',marker = 'x')
plt.hlines(-98961.7,xmin = -1000, xmax = 101000, lw = 0.5, linestyles = '--',label ="L2 = -98961.7",color ='black')
plt.legend(prop = {'size' : 15 })

plt.show()