import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt

df = pd.read_csv('Data_K_N.csv')
a =df.values


plt.figure(figsize = (16,9))
plt.title('Figure 5: Option values with different K & N ',fontsize = 20 )
plt.xlabel('Path number N',fontsize =12)
plt.ylabel('Option values',fontsize = 12)
plt.plot(a[0:20,1],a[0:20,2],label = 'K = 35 ',lw = 0.5,marker = 'x')
#plt.plot(a[20:40,1],a[20:40,2],label = 'K = 70 ',lw = 0.5,marker = 'x')
#plt.plot(a[40:60,1],a[40:60,2],label = 'K = 105 ',lw = 0.5,marker = 'x')
plt.plot(a[60:80,1],a[60:80,2],label = 'K = 140 ',lw = 0.5,marker = 'x')
#plt.plot(a[80:100,1],a[80:100,2],label = 'K = 175 ',lw = 0.5,marker = 'x')
#plt.plot(a[100:120,1],a[100:120,2],label = 'K = 210 ',lw = 0.5,marker = 'x')
plt.plot(a[120:140,1],a[120:140,2],label = 'K = 245 ',lw = 0.5,marker = 'x')
#plt.plot(a[140:160,1],a[140:160,2],label = 'K = 280 ',lw = 0.5,marker = 'x')
#plt.plot(a[160:180,1],a[160:180,2],label = 'K = 315 ',lw = 0.5,marker = 'x')
plt.plot(a[180:200,1],a[180:200,2],label = 'K = 350 ',lw = 0.5,marker = 'x')
plt.hlines(4900, 1000, 20000, lw =0.6, color = 'black', label = 'L4 = 4900', linestyles = '--')
plt.legend(loc = 'upper left', bbox_to_anchor = (1.02,1.02),prop ={'size' : 15 })
plt.show()