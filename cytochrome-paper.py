import sys
import os
sys.path.insert(1,os.getcwd()+'\\libraries')

import cytochrome_lib #This is a cytochrome library
import matplotlib.pyplot as plt
plt.close("all") #close all preexisting figures
import numpy as np

from mpl_toolkits.axes_grid1 import make_axes_locatable

version = "Last update: Aug 22, 2017"
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
exclusions_lst.append([1,1,1,1,0,0,0,0])
exclusions_lst.append([1,1,1,1,0,0,1,1])
exclusions_lst.append([1,1,1,1,1,1,1,1])
cyt_dict = {}
cyt_dict['b6f np'] = []
for excl in exclusions_lst:
    cyt_dict['b6f np'].append(cytochrome_lib.cyt('inputs\cytochrome_b6f.txt',excl))
    
exclusions_lst = []
exclusions_lst.append([0,0,0,0,0,0,0,0])
exclusions_lst.append([0,0,1,1,0,0,0,0])
exclusions_lst.append([0,0,1,1,0,0,1,1])
exclusions_lst.append([1,1,1,1,0,0,1,1])
exclusions_lst.append([1,1,1,1,1,1,1,1])

cyt_dict['b6f nn'] = []
for excl in exclusions_lst:
    cyt_dict['b6f nn'].append(cytochrome_lib.cyt('inputs\cytochrome_b6f.txt',excl))

exclusions_lst_pairs = []
#exclusions_lst_pairs.append([0,0,0,0,0,0,0,0])
exclusions_lst_pairs.append([0,0,1,1,0,0,1,1])
exclusions_lst_pairs.append([1,1,1,1,0,0,0,0])
exclusions_lst_pairs.append([0,0,1,1,1,1,0,0])
exclusions_lst_pairs.append([0,0,0,0,1,1,1,1])
exclusions_lst_pairs.append([1,1,0,0,0,0,1,1])
exclusions_lst_pairs.append([1,1,0,0,1,1,0,0])
exclusions_lst_pairs.append([1,1,1,1,1,1,1,1])

cyt_dict['b6f pairs'] = []
for excl in exclusions_lst_pairs:
    cyt_dict['b6f pairs'].append(cytochrome_lib.cyt('inputs\cytochrome_b6f.txt',excl))
    
x_range_nm = cyt_dict['b6f nn'][0].x_range_nm

plt.figure(1)
plt.ion()
plt.subplot(2,2,1)

for i in range(len(exclusions_lst)):
    plt.plot(x_range_nm,np.sum(cyt_dict['b6f nn'][i].specR,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)
#plt.legend(['n1p1','n1n2','n1p2','p1n2','p1p2','n2p2']);
plt.title('cytochrome b6f np model')

plt.subplot(2,2,2)
for i in range(len(exclusions_lst)):
    plt.plot(x_range_nm,np.sum(cyt_dict['b6f np'][i].specR,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)


plt.title('cytochrome b6f nn model')

plt.subplot(2,2,3)
for i in range(len(exclusions_lst)):
    plt.plot(x_range_nm,np.sum(cyt_dict['b6f np'][i].specD,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)



plt.subplot(2,2,4)
for i in range(len(exclusions_lst)):
    plt.plot(x_range_nm,np.sum(cyt_dict['b6f np'][i].specD,axis = 0),linewidth=2)


plt.show()


length = 10000
population = cytochrome_lib.kinetics_solve(np.array([1,1,1,1,0,0,0]),length)
x_time = np.linspace(0, 10, 10000)
plt.figure(2, figsize=(3, 5))
plt.ion()
sum_population = []
for i in range(5):
    plt.plot(x_time,population[i,:])
    sum_population.append(population[i,:]*i)
plt.plot(x_time,np.sum(sum_population,axis = 0)/4)
plt.title("Population distribution of proteins in different oxydation states")
plt.legend(['$n_0$ (Oxydized)','$n_1$','$n_2$','$n_3$','$n_4$(Reduced)','$n_{\sum i}$'])
plt.xlabel('1/k, seconds')
plt.ylabel('Populations')
#plt.subplots_adjust(bottom=0.14)   # <--
#plt.tight_layout() 
plt.show()


Absorbance_lst_b6f_nn = []
Circular_Dichroism_lst_b6f_nn = []
for i in range(5):
    Absorbance_lst_b6f_nn.append(population[i,:]*np.sum(np.sum(cyt_dict['b6f nn'][i].specD,axis = 0)))
    Circular_Dichroism_lst_b6f_nn.append(population[i,:]*np.sum(np.abs(np.sum(cyt_dict['b6f nn'][i].specR,axis = 0))))
Absorbance_b6f_nn = np.asarray(Absorbance_lst_b6f_nn)
Circular_Dichroism_b6f_nn = np.asarray(Circular_Dichroism_lst_b6f_nn)

Absorbance_lst_b6f_np = []
Circular_Dichroism_lst_b6f_np  = []
for i in range(5):
    Absorbance_lst_b6f_np.append(population[i,:]*np.sum(np.sum(cyt_dict['b6f np'][i].specD,axis = 0)))
    Circular_Dichroism_lst_b6f_np.append(population[i,:]*np.sum(np.abs(np.sum(cyt_dict['b6f np'][i].specR,axis = 0))))
Absorbance_b6f_np = np.asarray(Absorbance_lst_b6f_np)
Circular_Dichroism_b6f_np = np.asarray(Circular_Dichroism_lst_b6f_np)


plt.figure(3, figsize=(3, 3))

plt.ion()
plt.title('cytochrome b6f nn and np models')
plt.plot(x_time,np.sum(Absorbance_b6f_nn, axis = 0)/np.max(np.sum(Absorbance_b6f_nn, axis = 0)))
plt.plot(x_time,np.sum(Absorbance_b6f_np, axis = 0)/np.max(np.sum(Absorbance_b6f_np, axis = 0)))
plt.plot(x_time,np.sum(Circular_Dichroism_b6f_nn, axis = 0)/np.max(np.sum(Circular_Dichroism_b6f_nn, axis = 0)))
plt.plot(x_time,np.sum(Circular_Dichroism_b6f_np, axis = 0)/np.max(np.sum(Circular_Dichroism_b6f_np, axis = 0)))
plt.legend(['OD_nn','OD_np','CD_nn','CD_np'])
plt.xlabel('1/k, seconds')
plt.ylabel('Abs and CD, a.u.')
plt.subplots_adjust(bottom=0.14)   # <--
plt.tight_layout() 
plt.show()

plt.figure(4, figsize=(3, 4))

plt.figure(4)

plt.ion()
plt.subplot(2,1,1)
plt.title('Stick spectrum')
markerline, stemlines, baseline  = plt.stem(np.round(10**7/cyt_dict['b6f np'][4].Eigenvalues,2),cyt_dict['b6f np'][4].DipoleS, markerfmt=' ')
plt.setp(stemlines, 'color', 'r', 'linewidth', 5)
plt.setp(baseline, 'color', 'k', 'linewidth', 1)
plt.legend(['b6f'])
plt.ylabel('dipole strength, a.u.')


plt.subplot(2,1,2)

markerline, stemlines, baseline  = plt.stem(np.round(10**7/cyt_dict['b6f np'][4].Eigenvalues,2),cyt_dict['b6f np'][4].RotationalS, markerfmt=' ')
plt.setp(stemlines, 'color', 'r', 'linewidth', 5)
plt.setp(baseline, 'color', 'k', 'linewidth', 1)
plt.legend(['b6f'])
plt.xlabel('excitonic state position, nm')
plt.ylabel('rotational strength, a.u.')
plt.subplots_adjust(bottom=0.14)   # <--
plt.tight_layout() 
plt.show()


plt.figure(5)

plt.ion()

plt.title('stick spectrum state positions')



plt.plot(np.round(10**7/cyt_dict['b6f np'][4].Eigenvalues,2),'-o')
plt.legend(['b6f'])
plt.xlabel('excitonic state')
plt.ylabel('state position in wavelength, nm')
plt.axhline(y=433, xmin=0, xmax=1, hold=None)
plt.show()

plt.figure(6)

plt.subplot(3,2,(1,2))
plt.title('rotational strength vs excitonic state')
plt.plot(cyt_dict['b6f np'][4].RotationalS,'-o')

plt.legend(['b6f'])

plt.axhline(y=0, xmin=0, xmax=1, hold=None)


plt.subplot(3,2,(5,6))
plt.plot(10**7/cyt_dict['b6f np'][4].Eigenvalues,'-o')
plt.axhline(y=433, xmin=0, xmax=1, hold=None)
plt.show()





plt.figure(7, figsize=(3, 4))
plt.title('Eigenvectors')
lst = ['b6f nn']
ax = plt.gca()

im = ax.imshow(100*np.multiply(cyt_dict[lst[0]][4].Eigenvectors,cyt_dict[lst[0]][4].Eigenvectors), cmap=plt.cm.Blues, alpha=1.0)
plt.title('')
plt.yticks(range(8),  ['Bx_p1', 'By_p1', 'Bx_n1', 'By_n1', 'Bx_p2', 'By_p2', 'Bx_n2', 'By_n2'] ) 
plt.xticks(range(8),np.round(10**7/cyt_dict['b6f np'][4].Eigenvalues,2),rotation = 90)
plt.xlabel('excitonic state wavelength, nm')
plt.ylabel('contribution from a site state')


# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)

plt.show()
plt.subplots_adjust(bottom=0.14)   # <--
plt.tight_layout() 
plt.savefig('test.png', dpi=600) 




plt.figure(8, figsize=(3, 4))

plt.ion()
plt.subplot(2,1,1)
plt.title('Dipole and Rotional strength')

plt.plot(np.round(10**7/cyt_dict['b6f np'][4].Eigenvalues,2),cyt_dict['b6f np'][4].DipoleS,'-o')
plt.legend(['b6f'])
plt.ylabel('dipole strength, a.u.')


plt.subplot(2,1,2)

plt.plot(np.round(10**7/cyt_dict['b6f np'][4].Eigenvalues,2),cyt_dict['b6f np'][4].RotationalS,'-o')
plt.legend(['b6f'])
plt.xlabel('excitonic state position, nm')
plt.ylabel('rotational strength, a.u.')

plt.subplots_adjust(bottom=0.14)   # <--
plt.tight_layout() 
plt.show()




plt.figure(9,figsize=(3,4))
plt.ion()
plt.subplot(2,1,1)

for i in range(len(exclusions_lst_pairs)):
    plt.plot(x_range_nm,np.sum(cyt_dict['b6f pairs'][i].specD,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)
plt.legend(['n1n2','n1p1','n1p2','n2p2','n2p1','p1p2'], loc = 'right')

plt.ylabel('Absorbance, a.u.')

plt.subplot(2,1,2)
for i in range(len(exclusions_lst_pairs)):
    plt.plot(x_range_nm,33000*np.sum(cyt_dict['b6f pairs'][i].specR,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)

plt.ylabel('Circular Dihroism, mdeg')
plt.xlabel('wavelenth,nm')

plt.subplots_adjust(bottom=0.14)   # <--
plt.tight_layout() 
plt.show()


plt.figure(10)
plt.ion()
plt.subplot(2,1,1)

for i in range(len(exclusions_lst_pairs)):
    plt.plot(x_range_nm,np.sum(cyt_dict['b6f pairs'][i].specD,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)
plt.legend(['n1n2','n1p1','n1p2','n2p2','n2p1','p1p2','n1p1n2p2'])

plt.ylabel('Absorbance, a.u.')

plt.subplot(2,1,2)
for i in range(len(exclusions_lst_pairs)):
    plt.plot(x_range_nm,33000*np.sum(cyt_dict['b6f pairs'][i].specR,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)

plt.ylabel('Circular Dihroism, mdeg')
plt.xlabel('wavelenth,nm')

plt.show()

print "\nCalculations are finished. Please, see figures 1-3"

