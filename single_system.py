import numpy as np
import matplotlib.pyplot as plt

def onesweep(n,N,T,J=1,K=1):
    
    """   
    Parameters:
        n: An array of random intial state of system.
        N: The size of the Square Grid(N x N).
        T: Temperature.
        J: J is the coupling constant —
           it sets the strength of interaction between neighboring spins.
           When J>0 the system is ferromagnetic (spins want to align).
        K: Boltzman Constant.
        
    I used the checkerboard technique for flipping the spins, i.e, first only black squares are flipped,
    and then white squares. The condition for flipping is either the energy change should be negative or
    if the energy is positive then randomly it is selected with the probability of exp(-deltaE/TK) if it
    is going to flip.
    
    When every point is allowed to flip once, it is known as One Sweep of Simulation.
    
    Returns:
        The modified grid n with flipped spins.
        
    """
    
    # Computing the Energy change for every black square on the checkerboard using the formula: E(i,j) = 2*J*s(i,j)*(sum of adjacent squares(whites for this case))
    deltaE_b1 = 2*J*n[::2,1::2] * (np.roll(n,1,axis=0)[::2,1::2] + np.roll(n,-1,axis=0)[::2,1::2] + np.roll(n,1,axis=1)[::2,1::2] + np.roll(n,-1,axis=1)[::2,1::2])
    deltaE_b2 = 2*J*n[1::2,::2] * (np.roll(n,1,axis=0)[1::2,::2] + np.roll(n,-1,axis=0)[1::2,::2] + np.roll(n,1,axis=1)[1::2,::2] + np.roll(n,-1,axis=1)[1::2,::2])
    
    n[::2,1::2] *= -1  #Flipping the spins before checking:
    n[1::2,::2] *= -1  # Because we can check the ones that are rejected and unflip them.
    
    b1 = n[::2,1::2].copy()
    b2 = n[1::2,::2].copy()
    
    # Flipping back the rejected spins.
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
    
def simulation(n,N,t,burnin,steps, J=1.0, K=1.0):
    
    """
    Parameters:
        n: An array of random intial state of system.
        N: The size of the Square Grid(N x N).
        t: Temperature
        burnin: The initial steps required for the system to stabalize are insignificant
                for our measurements of observables so we skip them.
        steps: The number of steps after the burnin steps.
        
    Three observable's list is initiated at start 'Absolute Magnetization values', 'Magnetization values',
    and 'Energy values' for different states of the system.
    
    After skippin the burnin steps we let the system undergo the main steps and calculate the Energy
    and total Magnetization.
    

    Returns:
        np.mean(absmvalue): The average of the absolute magnetization value over all the steps.
        np.var(evalue)/((N**2)*(t**2)): The Specific Heat of the system over all the steps.
                                        From the following formula: C = Var(E) / (N^2 * T^2).
        np.var(mvalue)/((N**2)*t): The Magnetic Susceptibility of the system over all the steps.
                                   From the following formula: X = Var(M) / (N^2 * T).
    """
    
    absmvalue = []
    mvalue = []
    evalue = []
    
    s = 0

    while s < burnin:
        
        s += 1
        
        onesweep(n,N,t,J,K)

    for k in range(steps):
        
        onesweep(n,N,t,J,K)
    
        E = -J*n*(np.roll(n,1,axis=0) + np.roll(n,1,axis=1)) 
        m = np.sum(n)
        
        evalue.append(np.sum(E))
        absmvalue.append(abs(m))
        mvalue.append(m)
    
    return np.mean(absmvalue), np.var(evalue)/((N**2)*(t**2)), np.var(mvalue)/((N**2)*t)

def simulate_temperatures(n,N,temp,burnin,steps):
    
    """   
    Parameters:
        n: An array of random intial state of system.
        temp: An array of temperatures for the system to simulate the process.
        N: The size of the Square Grid(N x N).
        burnin: The initial steps required for the system to stabalize are insignificant
                for our measurements of observables so we skip them.
        steps: The number of steps after the burnin steps.
        
    Simulates the process over a range of temperatures provided and record the average Magnetization,
    Specific Heat and Magnetic Susceptibility of the System at every temperatures in a list.

    Returns:
        mavg: The average Magnetization of the system over a range of temperatures.
        C: The Specific Heat of the system over a range of temperatures.
        X: The Magnetic Susceptibility of the system over a range of temperatures.
    """
    
    mavg = []
    C = []
    X = []
    
    for t in temp:
        
        mag, specifich, suscep = simulation(n,N,t,burnin,steps)
        
        mavg.append(mag)
        C.append(specifich)
        X.append(suscep)
    
    return mavg, C, X
    
def plot(temp,N,mavg,C,X):
    
    """Plots the Magnetization, Specific Heat and Magnetic Susceptibility over the 
       range of the temperatures.
    """
    
    fig, (ax1, ax2, ax3) = plt.subplots(1,3)

    ax1.plot(temp,np.array(mavg)/(N**2))
    ax1.axvline(2.2691853, color='red', label="Onsager's Tc value")
    ax1.set_xlabel('Temperature')
    ax1.set_ylabel('Magnetization')
    ax1.set_title('Temperature vs Magnetization')

    ax2.plot(temp, C)
    ax2.axvline(2.2691853, color='red', label="Onsager's Tc value")
    ax2.set_xlabel('Temperature')
    ax2.set_ylabel('Specific Heat')
    ax2.set_title('Temperature vs Specific Heat')

    ax3.plot(temp, X)
    ax3.axvline(2.2691853, color='red', label="Onsager's Tc value")
    ax3.set_xlabel('Temperature')
    ax3.set_ylabel('Susceptibility')
    ax3.set_title('Temperature vs Susceptibility')
    
    plt.style.use('seaborn-v0_8')
    plt.legend()
    plt.tight_layout()
    plt.show()
            
def main():
    
    N = 40
    K = 1.0
    J = 1.0
    burnin = 300
    steps = 800
    temp = np.arange(1,4,0.01)
    n = np.ones([N,N], dtype=int)
    
    mavg, C, X =simulate_temperatures(n,N,temp,burnin,steps)
    
    plot(temp,N,mavg,C,X)

if __name__ == "__main__":
    main()