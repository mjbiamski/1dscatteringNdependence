#import libraries
import numpy as np
import matplotlib.pyplot as plt

#generate empty arrays
T_tot = []
T_sum = []
T_avg = []
n = 0 #define counter
realizations = 10
potentials = 1000000

while n < realizations:
    E = (27*np.pi/(5*2))**2 #define energy to be used for kd value
    N = 1
    g_array = [] #define empty random number array
    N_array = [] #define empty array for potentials
    kd_array = [] #define empty array for kd values
    MDtot_array = [] #define empty total transfer matrix array
    T_array = [] #define transmission coefficient array 
    g = 0  
    for i in range(potentials):
        g = np.random.uniform(-0.01,0.01) #generate random number
        g_array.append(g) #populate g_array
    while N < len(g_array):
        N_array.append(N) #populate N_array
        N += 1
    MDtot = np.matrix([[1,0],[0,1]])
    for i in range(len(N_array)):
        D = 2
        d = 1
        k = D*np.sqrt(E)/d
        B = 2*D/d
        #define transfer matrix for one potential
        MD0 = np.matrix([[(1-((1j*B)/(2*k)))*np.exp(-1j*d*(1+g_array[i])*k), 
                          (-1j*B)/(2*k)], 
                       [(1j*B)/(2*k), 
                        (1+((1j*B)/(2*k)))*np.exp(1j*d*(1+g_array[i])*k)]]) 
        MDtot = np.matmul(MDtot,MD0) #compute total transfer matrix
        MDtot_array.append(MDtot) #populate total transfer matrix array
    for i in range(len(N_array)):
        MD = MDtot_array[i]
        mdtot11 = MD[1,1]
        mdtot11star = np.conj(mdtot11)
        t1 = 1/(mdtot11*mdtot11star) #compute transmission coefficient
        T_array.append(t1) #populate transmission coefficient array
    T_tot.append(T_array)
    n += 1 #increment counter

T_sum = [sum(i) for i in zip(*T_tot)] #sum transmission coefficients

#average transmission coefficients for realizations
for i in range(len(T_sum)):
    T_avg.append(T_sum[i]/realizations) 

T_totB = []
T_sumB = []
T_avgB = []
n = 0
realizations = 10
potentials = 1000000

while n < realizations:
    E = (27*np.pi/(5*2))**2
    N = 1
    g_array = []
    N_array = []
    kd_array = []
    MDtot_array = []
    T_array = []
    g = 0  
    for i in range(potentials):
        g = np.random.uniform(-0.05,0.05)
        g_array.append(g)
    while N < len(g_array):
        N_array.append(N)
        N += 1
    MDtot = np.matrix([[1,0],[0,1]])
    for i in range(len(N_array)):
        D = 2
        d = 1
        k = D*np.sqrt(E)/d
        B = 2*D/d
        MD0 = np.matrix([[(1-((1j*B)/(2*k)))*np.exp(-1j*d*(1+g_array[i])*k), 
                          (-1j*B)/(2*k)], 
                       [(1j*B)/(2*k), 
                        (1+((1j*B)/(2*k)))*np.exp(1j*d*(1+g_array[i])*k)]])
        MDtot = np.matmul(MDtot,MD0)
        MDtot_array.append(MDtot)
    for i in range(len(N_array)):
        MD = MDtot_array[i]
        mdtot11 = MD[1,1]
        mdtot11star = np.conj(mdtot11)
        t1 = 1/(mdtot11*mdtot11star)
        T_array.append(t1)
    T_totB.append(T_array)
    n += 1

T_sumB = [sum(i) for i in zip(*T_totB)]

for i in range(len(T_sumB)):
    T_avgB.append(T_sumB[i]/realizations)


T_totC = []
T_sumC = []
T_avgC = []
n = 0
realizations = 10
potentials = 1000000

while n < realizations:
    E = (27*np.pi/(5*2))**2
    N = 1
    g_array = []
    N_array = []
    kd_array = []
    MDtot_array = []
    T_array = []
    g = 0  
    for i in range(potentials):
        g = np.random.uniform(-0.1,0.1)
        g_array.append(g)
    while N < len(g_array):
        N_array.append(N)
        N += 1
    MDtot = np.matrix([[1,0],[0,1]])
    for i in range(len(N_array)):
        D = 2
        d = 1
        k = D*np.sqrt(E)/d
        B = 2*D/d
        MD0 = np.matrix([[(1-((1j*B)/(2*k)))*np.exp(-1j*d*(1+g_array[i])*k), 
                          (-1j*B)/(2*k)], 
                       [(1j*B)/(2*k), 
                        (1+((1j*B)/(2*k)))*np.exp(1j*d*(1+g_array[i])*k)]])
        MDtot = np.matmul(MDtot,MD0)
        MDtot_array.append(MDtot)
    for i in range(len(N_array)):
        MD = MDtot_array[i]
        mdtot11 = MD[1,1]
        mdtot11star = np.conj(mdtot11)
        t1 = 1/(mdtot11*mdtot11star)
        T_array.append(t1)
    T_totC.append(T_array)
    n += 1

T_sumC = [sum(i) for i in zip(*T_totC)]

for i in range(len(T_sumC)):
    T_avgC.append(T_sumC[i]/realizations)

plt.figure(figsize=(8, 5),dpi=200)
plt.plot(N_array,T_avgC,label='10%',linewidth=0.7)
plt.plot(N_array,T_avgB,label='5%',linewidth=0.7)
plt.plot(N_array,T_avg,label='1%',linewidth=0.7)
plt.xscale("log")
plt.xlabel("N")
plt.ylabel("$T_{N}$")
plt.gca().add_artist(plt.legend(loc='lower right'))
plt.gca().add_artist(plt.legend(["$kd$ = $27 \pi/5$"],handlelength=0,loc='upper right'))
plt.ylim([0,1])
plt.show()