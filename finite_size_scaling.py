import numpy as np
import matplotlib.pyplot as plt
from single_system import onesweep

def simulation_ki(n,N,t,burnin,steps,ki=10,J=1,K=1):
    
    """
    Same structure as the 'simulation' function in 'single_system.py' but only measures
    the magnetization and only every ki steps to reduce autocorrelation.
    
    Parameters:
        ki: Number of steps after which the magnetization of the system is measured.
    

    Returns:
        np.mean(absmvalue): The average of the absolute magnetization value over all the steps.
    """
    
    absmvalue = []
        
    s = 0

    while s < burnin:
            
        s += 1
            
        onesweep(n,N,t,J,K)
            
    _ = 0

    for k in range(steps):
            
        _ += 1
            
        onesweep(n,N,t,J,K)
            
        if _%ki == 0:
            absmvalue.append(abs(np.sum(n)))
    
    return np.mean(absmvalue)

def run_and_plot(Nlist, temp, burnin, steps):
    
    """
    Parameters:
        Nlist: A list of different grid sizes for different systems.
        temp: An array of temperatures for the system to simulate the process.
        burnin: The initial steps required for the system to stabalize are insignificant
                for our measurements of observables so we skip them.
        steps: The number of steps after the burnin steps.
        
    Simulates the process for different grid sizes over a range of temperatures, and
    plots the Magnetization for every temperature for single grid and moves to next
    grid size.   
    
    """
    
    for N in Nlist:

        n = np.ones([N,N], dtype=int)

        mavg = []

        for t in temp:
            
            mag = simulation_ki(n,N,t,burnin,steps)
            
            mavg.append(mag)
            
        plt.plot(temp,np.array(mavg)/(N**2), label=N)
    
def show_plot():
    
    """
    Adds labels and title to the plot and shows the plot.
    """
    
    plt.axvline(2.2691853, label="Onsager's Tc value")
    plt.xlabel('Temperature')
    plt.ylabel('Magnetization')
    plt.title('Temperature vs Magnetization')
    plt.legend()    
    plt.show()          

def main():
    Nlist = [10,20,30,40]
    temp = np.arange(1,4,0.01)
    K = 1.0
    J = 1.0
    burnin = 300
    steps = 8000
    
    run_and_plot(Nlist,temp,burnin,steps)
    
    show_plot()
    
if __name__ == "__main__":
    main()