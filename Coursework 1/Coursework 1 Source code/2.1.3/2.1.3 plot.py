import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt

a = pd.read_csv('Data_CI.csv')
b = pd.read_csv('Data_CI2.csv')
a = a.values
b= b.values


plt.figure(figsize = (16,9))
plt.title('Figure 3: Confidence interval against N for two methods',fontsize = 20 )
plt.xlabel('Path number N',fontsize = 12)
plt.ylabel('Confidence interval',fontsize = 12)
plt.vlines(a[:,0],a[:,1],a[:,2],color ='black',lw = 1.0 ,linestyles = '--' ,label = 'Intervals(Basic)')
plt.vlines(b[:,0],b[:,1],b[:,2],color ='#FF8080',lw = 1.0 ,linestyles ='--',label = 'Intervals(Antithetic)')
plt.hlines(-660.801,-1000,101000,lw = 0.5,color ='blue',label ='L3 = -660.801')
#plt.plot(a[:,0],a[:,1],color = '#FF8080',lw = 0.8,label = 'Bounds(Basic)' )
#plt.plot(a[:,0],a[:,2],color='#FF8080',lw = 0.8 )
plt.plot(b[:,0],b[:,1],color = '#FF8080',lw = 0.8,label = 'Bounds(Antithetic) ' )
plt.plot(b[:,0],b[:,2],color='#FF8080',lw = 0.8 )
plt.legend(prop ={'size' : 15 })

plt.show()