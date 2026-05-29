import numpy as np
import matplotlib.pyplot as plt
import random
import time
start = time.time()

Nlist = [10,20,30,40]
K = 1.0
T = 1.0
J = 1.0
burnin = 300



def onesweep(N,T,n):
    
    
    
    deltaE_b1 = 2*J*n[::2,1::2] * (np.roll(n,1,axis=0)[::2,1::2] + np.roll(n,-1,axis=0)[::2,1::2] + np.roll(n,1,axis=1)[::2,1::2] + np.roll(n,-1,axis=1)[::2,1::2])
    deltaE_b2 = 2*J*n[1::2,::2] * (np.roll(n,1,axis=0)[1::2,::2] + np.roll(n,-1,axis=0)[1::2,::2] + np.roll(n,1,axis=1)[1::2,::2] + np.roll(n,-1,axis=1)[1::2,::2])
    
    n[::2,1::2] *= -1
    n[1::2,::2] *= -1
    
    b1 = n[::2,1::2].copy()
    b2 = n[1::2,::2].copy()
    
    b1[np.random.random(deltaE_b1.shape) > np.exp(-deltaE_b1/(K*T))] *= -1
    b2[np.random.random(deltaE_b2.shape) > np.exp(-deltaE_b2/(K*T))] *= -1
    
    n[::2,1::2] = b1 
    n[1::2,::2] = b2
    
    deltaE_w1 = 2*J*n[::2,::2] * (np.roll(n,1,axis=0)[::2,::2] + np.roll(n,-1,axis=0)[::2,::2] + np.roll(n,1,axis=1)[::2,::2] + np.roll(n,-1,axis=1)[::2,::2])
    deltaE_w2 = 2*J*n[1::2,1::2] * (np.roll(n,1,axis=0)[1::2,1::2] + np.roll(n,-1,axis=0)[1::2,1::2] + np.roll(n,1,axis=1)[1::2,1::2] + np.roll(n,-1,axis=1)[1::2,1::2])
    
    n[::2,::2] *= -1
    n[1::2,1::2] *= -1
    
    w1 = n[::2,::2].copy()
    w2 = n[1::2,1::2].copy()
    
    w1[np.random.random(deltaE_w1.shape) > np.exp(-deltaE_w1/(K*T))] *= -1
    w2[np.random.random(deltaE_w2.shape) > np.exp(-deltaE_w2/(K*T))] *= -1
    
    n[::2,::2] = w1
    n[1::2,1::2] = w2
    

for N in Nlist:

    n = np.ones([N,N], dtype=int)

    temp = np.arange(1,4,0.01)

    mavg = []
    eavg = []
    C = []
    X = []


    for t in temp:
        
        absmvalue = []
        mvalue = []
        evalue = []
        
        steps = 0

        while steps < burnin:
            
            steps += 1
            
            onesweep(N,t,n)
            
        ki = 0

        for k in range(8000):
            
            ki += 1
            
            onesweep(N,t,n)
            
            #E = -J*n*(np.roll(n,1,axis=0) + np.roll(n,1,axis=1))
            #evalue.append(np.sum(E))
            
            if ki%10 == 0:
                m = np.sum(n)
                absmvalue.append(abs(m))
                mvalue.append(m)
        
        #C.append(np.var(evalue)/((N**2)*(t**2)))
        #X.append(np.var(mvalue)/((N**2)*t))
        #eavg.append(np.mean(evalue))
        mavg.append(np.mean(absmvalue))
        
    plt.plot(temp,np.array(mavg)/(N**2), label=N)

plt.axvline(2.2691853)
plt.xlabel('Temperature')
plt.ylabel('Magnetization')
plt.title('Temperature vs Magnetization')
plt.legend()    
plt.show()
            
print(time.time() - start)

    