import numpy as np
import datetime
from fractions import Fraction 
import matplotlib.pyplot as plt
import warnings

#The program below is an updated version of the code I originally developed for my senior research project.  
#It is more efficient, utilizes classes, and allows users to input parameters for the calculation.

#See ElectronicTransmissionPeriodicSolid.pdf in https://mjbiamskiorg.wordpress.com/code/ for a detailed description my research as well as a derivation of the physics used below.

#This program calculates the probability of an electron passing through a 1D crystal with a disordered 
# lattice structure (disorder determined using random number). 
#The calculation is performed for a single dimensionless energy value (which is a multiple of pi defined here as evalue). 
#The program allows for the user to input parameters where potentials represent the number of 
# lattices, disorder the non-uniformity in the lattice spacing, and realizations the number of calculations for 
# a given disorder (realizations are averaged for final caluclation).
#Once the calculations have been performed a single plot is created displaying the results.   


warnings.filterwarnings("ignore", category=RuntimeWarning)

#To see a clear depiction of how disorder effects transmission probability enter potentials > 10^3 
#Average runtime for 10^3 potentials is 0.02 s per realization


#collect input data
disorder_array = list(map(float, input("Enter up to 5 different amounts of disorder (as decimal percentages, separated by a space): ").split(" ")))
evalue = Fraction(input("Enter the desired energy value (as int multiple or fraction of pi): "))
pots = int(input("Enter number of potentials: "))
reals = int(input("Enter number of realizations for each set of data: "))


class TransmissionProb:
    def __init__(self, pevalue, energy, potentials, realizations):
        self.pevalue = pevalue
        self.energy = energy
        self.potentials = potentials
        self.realizations = realizations

        energyvalue = (np.pi*self.energy/2)**2 #define energy to be used for kd value

        #generate empty arrays
        T_tot = []
        T_sum = []
        self.T_avg = []
        
        n = 0 #define counter
        
        print("-----------------------")
        print("At {p} percent error ".format(p=self.pevalue*100))
        print("-----------------------")

        while n < self.realizations:
            start = datetime.datetime.now()
            N = 1
            g_array = [] #define empty random number array
            self.N_array = [] #define empty array for potentials
            MDtot_array = [] #define empty total transfer matrix array
            T_array = [] #define transmission coefficient array 
            g = 0  
            for i in range(self.potentials):
                g = np.random.uniform(-1*self.pevalue,self.pevalue) #generate random number
                g_array.append(g) #populate g_array
            while N < len(g_array):
                self.N_array.append(N) #populate N_array
                N += 1
            MDtot = np.matrix([[1,0],[0,1]])
            for i in range(len(self.N_array)):
                D = 2
                d = 1
                k = D*np.sqrt(energyvalue)/d
                B = 2*D/d
                #define transfer matrix for one potential
                MD0 = np.matrix([[(1-((1j*B)/(2*k)))*np.exp(-1j*d*(1+g_array[i])*k), 
                                (-1j*B)/(2*k)], 
                            [(1j*B)/(2*k), 
                                (1+((1j*B)/(2*k)))*np.exp(1j*d*(1+g_array[i])*k)]]) 
                MDtot = np.matmul(MDtot,MD0) #compute total transfer matrix
                MDtot_array.append(MDtot) #populate total transfer matrix array
            for i in range(len(self.N_array)):
                MD = MDtot_array[i]
                mdtot11 = MD[1,1]
                mdtot11star = np.conj(mdtot11)
                t1 = 1/(mdtot11*mdtot11star) #compute transmission coefficient
                T_array.append(t1) #populate transmission coefficient array
            T_tot.append(T_array)
            n += 1 #increment counter
            print(str(n) + " realization(s): " + str(datetime.datetime.now()-start) + 
                " seconds") #show time per realization

        T_sum = [sum(i) for i in zip(*T_tot)] #sum transmission coefficients

        #average transmission coefficients for realizations
        for i in range(len(T_sum)):
            self.T_avg.append(T_sum[i]/self.realizations) 



transprob_array = [[] for i in range(len(disorder_array))]

for i in range(len(disorder_array)):
    disordervalues = TransmissionProb(disorder_array[i], evalue, pots, reals)
    transprob_array[i].append(disordervalues.N_array)
    transprob_array[i].append(disordervalues.T_avg)



plt.figure(figsize=(8, 5),dpi=200)

for i in range(len(disorder_array)):
    plt.plot(transprob_array[i][0],transprob_array[i][1],label=str(str(disorder_array[i]*100)+"%"),linewidth=0.7)

plt.xscale("log")
plt.xlabel("N")
plt.ylabel("$T_{N}$")
plt.gca().add_artist(plt.legend(loc='center left'))
plt.gca().add_artist(plt.legend(["$kd$ =" +str(evalue) + "$\pi$"],handlelength=0,loc='lower left'))
plt.show()