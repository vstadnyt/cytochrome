import cytochrome_lib #This is a cytochrome library
import matplotlib.pyplot as plt
import numpy as np



version = "Last update: Aug 8, 2017"
desription = "This code calculates population distribution in the cytochrome b6f protein and plots kinetic profiles for two different models: \n'nn' and 'np' models \n The outputs are: \n Figure 1: \n Figure 2: The ppulation distributions for different oxydations states of the cytochrome proteins. \n Figure 3: the resulting absorbance and circular dichroism kinetics for two different models"
print desription
print version


#the eclusions_lst is a list of hemes that are taken into account during calculations (1 - include; 0 - exclude); 
#There are 8 values for 4 hemes and 2 dipoles per heme: [Qx_p1, Qy_p1, Qx_n1, Qy_n1, Qx_p2, Qy_p2, Qx_n2, Qy_n2]


##This is a main part of a code
#This part creates two lists of several instances of a cyt class (see cytochrome library) with different input files
exclusions_lst = []
exclusions_lst.append([0,0,0,0,0,0,0,0])
exclusions_lst.append([0,0,1,1,0,0,0,0])
exclusions_lst.append([0,0,1,1,0,0,1,1])
exclusions_lst.append([1,1,1,1,0,0,1,1])
exclusions_lst.append([1,1,1,1,1,1,1,1])
cyt_b6f_np = []
for excl in exclusions_lst:
    cyt_b6f_np.append(cytochrome_lib.cyt('cytochrome_b6f.txt',excl))
for i in range(len(exclusions_lst)):
    cyt_b6f_np[i].read_structure_file()
    cyt_b6f_np[i].Hamiltonian()
    cyt_b6f_np[i].D_and_R_strength()
    cyt_b6f_np[i].spectra_plot()

exclusions_lst = []
exclusions_lst.append([0,0,0,0,0,0,0,0])
exclusions_lst.append([0,0,1,1,0,0,0,0])
exclusions_lst.append([1,1,1,1,0,0,0,0])
exclusions_lst.append([1,1,1,1,0,0,1,1])
exclusions_lst.append([1,1,1,1,1,1,1,1])
cyt_b6f_nn = []
for excl in exclusions_lst:
    cyt_b6f_nn.append(cytochrome_lib.cyt('cytochrome_b6f.txt',excl))
for i in range(len(exclusions_lst)):
    cyt_b6f_nn[i].read_structure_file()
    cyt_b6f_nn[i].Hamiltonian()
    cyt_b6f_nn[i].D_and_R_strength()
    cyt_b6f_nn[i].spectra_plot()


x_range_nm = cyt_b6f_nn[0].x_range_nm

plt.figure(1)
plt.ion()
plt.subplot(2,2,1)

for i in range(len(exclusions_lst)):
    plt.plot(x_range_nm,np.sum(cyt_b6f_nn[i].specR,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)
#plt.legend(['n1p1','n1n2','n1p2','p1n2','p1p2','n2p2']);
plt.title('cytochrome b6f np model')

plt.subplot(2,2,2)
for i in range(len(exclusions_lst)):
    plt.plot(x_range_nm,np.sum(cyt_b6f_np[i].specR,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)


plt.title('cytochrome b6f nn model')

plt.subplot(2,2,3)
for i in range(len(exclusions_lst)):
    plt.plot(x_range_nm,np.sum(cyt_b6f_nn[i].specD,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)



plt.subplot(2,2,4)
for i in range(len(exclusions_lst)):
    plt.plot(x_range_nm,np.sum(cyt_b6f_np[i].specD,axis = 0),linewidth=2)


plt.show()


length = 10000
population = cytochrome_lib.kinetics_solve(np.array([1,1,1,1,0,0,0]),length)

plt.figure(2)
plt.ion()
for i in range(5):
    plt.plot(range(length),population[i,:])
plt.title("Population distribution of proteins in different oxydation states")
plt.legend(['0e- state (fully oxydized)','1e- state','2e- state','3e- state','4e- state(fully reduced)'])
plt.show()



Absorbance_lst_b6f_nn = []
Circular_Dichroism_lst_b6f_nn = []
for i in range(5):
    Absorbance_lst_b6f_nn.append(population[i,:]*np.sum(np.sum(cyt_b6f_nn[i].specD,axis = 0)))
    Circular_Dichroism_lst_b6f_nn.append(population[i,:]*np.sum(np.abs(np.sum(cyt_b6f_nn[i].specR,axis = 0))))
Absorbance_b6f_nn = np.asarray(Absorbance_lst_b6f_nn)
Circular_Dichroism_b6f_nn = np.asarray(Circular_Dichroism_lst_b6f_nn)

Absorbance_lst_b6f_np = []
Circular_Dichroism_lst_b6f_np  = []
for i in range(5):
    Absorbance_lst_b6f_np.append(population[i,:]*np.sum(np.sum(cyt_b6f_np[i].specD,axis = 0)))
    Circular_Dichroism_lst_b6f_np.append(population[i,:]*np.sum(np.abs(np.sum(cyt_b6f_np[i].specR,axis = 0))))
Absorbance_b6f_np = np.asarray(Absorbance_lst_b6f_np)
Circular_Dichroism_b6f_np = np.asarray(Circular_Dichroism_lst_b6f_np)


plt.figure(3)

plt.ion()
plt.title('cytochrome b6f nn and np models')
plt.plot(range(length),np.sum(Absorbance_b6f_nn, axis = 0)/np.max(np.sum(Absorbance_b6f_nn, axis = 0)))
plt.plot(range(length),np.sum(Absorbance_b6f_np, axis = 0)/np.max(np.sum(Absorbance_b6f_np, axis = 0)))
plt.plot(range(length),np.sum(Circular_Dichroism_b6f_nn, axis = 0)/np.max(np.sum(Circular_Dichroism_b6f_nn, axis = 0)))
plt.plot(range(length),np.sum(Circular_Dichroism_b6f_np, axis = 0)/np.max(np.sum(Circular_Dichroism_b6f_np, axis = 0)))
plt.legend(['OD_nn','OD_np','CD_nn','CD_np'])

plt.show()

print "\nCalculations are finished. Please, see figures 1-3"

