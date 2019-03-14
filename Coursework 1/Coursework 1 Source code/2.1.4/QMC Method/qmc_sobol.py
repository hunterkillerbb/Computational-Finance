import matplotlib.pyplot as plt
import numpy             as np
np.random.seed(12345678)
from numpy.random import randn
from numpy import ones, zeros, sqrt, exp, log2
from scipy.stats import norm
from Option_selection import binary_call,european_put,european_call 
from bb import bb
from sobol import sobol, scramble 

def binary(S,K):
    temp = []
    for i in range(len(S)):
        if S[i] > K:
            temp.append(1)
        else:
            temp.append(0)
     
    return np.array(temp)
        
       
    
#
# Test weak convergence of Euler method for 
# European put option using Sobol points
#
# Test problem:   dS   = (r - D0)*S dt + sigma*S dW
#

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

plt.close("all")
plt.ion()

#
# problem parameters and exact solution
#

r   = 0.03
D   = 0.03 
sigma = 0.27 
T   = 1.75 
S0  = 60000
K_1 = 60000
K_2 = 105000

#Ms   = zeros((11,))
#err1 = zeros((11,))
#err2 = zeros((11,))
M_2 = []
V_2 = []

Ve  = european_put (r,D,sigma,T,S0,K_1) - european_call(r,D,sigma,T,S0,K_2) - 2 * K_2*binary_call(r,D,sigma,T,S0,K_2)
print(Ve)

for it in range(1,4):

    #
    # Quasi-Monte Carlo simulation comparing to exact solution
    #

    N  = 128   # number of timesteps
    M2 = 64    # number of scrambles

    h  = T/N

    for p in range(2,18):
        print(p)
        M = 2**p  # number of paths in each "family"
        M_2.append(M)
        unscrambled = sobol(m=p, s=N, scramble=False)

        sum1 = 0.
        sum2 = 0.

        for m in range(1,M2+1):
            #
            # Sobol quasi-random number generation with scrambling,
            # and Brownian Bridge construction of Brownian increments
            #
            if it==1:
                U = scramble(unscrambled).T
                Z  = norm.ppf(U)
                dW = bb(Z,T)

            #
            # Sobol quasi-random number generation with scrambling,
            # but no Brownian Bridge construction of Brownian increments
            #
            elif it==2:
                U = scramble(unscrambled).T
                Z  = norm.ppf(U)
                dW = sqrt(h)*Z

            #
            # alternative standard random number generation
            #
            else:
                dW = sqrt(h)*randn(N,M)
                
            K_1 = 60000
            K_2 = 105000
            S = S0*ones((M,))
            K_21 = K_2
            K_1 = K_1*ones((M,))
            K_2 = K_2*ones((M,))   # M = path number，N = intervals

            for n in range(N):
                S  = S*(1+r*h+sigma*dW[n,:])                

            P = exp(-r*T)*(np.maximum(K_1 - S,0.) - np.maximum(S - K_2 , 0.) - 2 * K_21 * binary(S,K_21))  #payoff
            P = np.sum(P)/M

            sum1 = sum1 + np.sum(P)
            sum2 = sum2 + np.sum(P**2)

        V  = sum1/M2
        V_2.append(V)
        #sd = sqrt((sum2/M2 - (sum1/M2)**2)/(M2-1))

    plt.figure(figsize = (11,6))
    plt.xlabel('Path number N');plt.ylabel('portfolio value')
    plt.plot( M_2,V_2, label = 'Value',lw = 0, marker = '.',markersize = 2)
    plt.legend(loc = 'upper left',fontsize = 14)
    if it == 1:
        plt.title('Figure xxx: Case--Sobol sequence(with BB)')
    elif it == 2:
        plt.title('Figure xxx: Case--Sobol sequence(with no BB)')
    else:
        plt.title('Figure xxx:Case--Peusdo random numbers')
    



        #Ms[p-2]   = M
       # err1[p-2] = V-Ve
        #err2[p-2] = 3*sd

    #plt.figure(figsize = (11,6))
    #plt.loglog(Ms,abs(err1),'b-x',Ms,err2,'r-x')
    #plt.xlabel('N'); plt.ylabel('Error')
    #plt.legend((' Error',' MC error bound'), loc='upper right', fontsize=14)
    #plt.axis([1, 4096, 0.005, 10])
    #if   it==1:
    
       # plt.title('comparison to exact solution')
        
   # elif it==2:
       # plt.title('comparison to exact solution')
        
    #else:
        #plt.title('comparison to exact solution')


       
        
# python2.7 and python3 compatibility
if hasattr(__builtins__, 'raw_input'):
    raw_input("Press Enter to continue...")
else:
    input("Press Enter to continue...")