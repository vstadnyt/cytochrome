import sys
import os
sys.path.insert(1,os.getcwd()+'\\libraries')

import cytochrome_lib #This is a cytochrome library
import matplotlib.pyplot as plt
import numpy as np



#the eclusions_lst is a list of hemes that are taken into account during calculations (1 - include; 0 - exclude); 
#There are 8 values for 4 hemes and 2 dipoles per heme: [Qx_p1, Qy_p1, Qx_n1, Qy_n1, Qx_p2, Qy_p2, Qx_n2, Qy_n2]

#here we are creating a list of desired combinations of different interaction heme pairs
exclusions_lst = []
exclusions_lst.append([1,1,1,1,0,0,0,0])
exclusions_lst.append([1,1,0,0,1,1,0,0])
exclusions_lst.append([1,1,0,0,0,0,1,1])
exclusions_lst.append([0,0,1,1,1,1,0,0])
exclusions_lst.append([0,0,1,1,0,0,1,1])
exclusions_lst.append([0,0,0,0,1,1,1,1])



##This is a main part of a code
#This part creates two lists of several instances of a cyt class (see cytochrome library) with different input files

cyt_dict = {}
cyt_dict['b6f'] = []
for excl in exclusions_lst:
    cyt_dict['b6f'] .append(cytochrome_lib.cyt('inputs\cytochrome_b6f.txt',excl))
for i in range(len(exclusions_lst)):
    cyt_dict['b6f'] [i].read_structure_file()
    cyt_dict['b6f'] [i].Hamiltonian()
    cyt_dict['b6f'] [i].D_and_R_strength()
    cyt_dict['b6f'] [i].spectra_plot()
    
cyt_dict['bc1'] = []
for excl in exclusions_lst:
    cyt_dict['bc1'].append(cytochrome_lib.cyt('inputs\cytochrome_bc1.txt',excl))
for i in range(len(exclusions_lst)):
    cyt_dict['bc1'][i].read_structure_file()
    cyt_dict['bc1'][i].Hamiltonian()
    cyt_dict['bc1'][i].D_and_R_strength()
    cyt_dict['bc1'][i].spectra_plot()




x_range_nm = cyt_dict['b6f'][0].x_range_nm

plt.figure(1)
plt.ion()
plt.subplot(2,2,1)

for i in range(len(exclusions_lst)):
    plt.plot(x_range_nm,np.sum(cyt_dict['b6f'][i].specR,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)
#plt.legend(['n1p1','n1n2','n1p2','p1n2','p1p2','n2p2']);
plt.title('cytochrome b6f')
plt.legend(['n1p1','n1n2','n1p2','p1n2','p1p2','n2p2'])
plt.subplot(2,2,2)
for i in range(len(exclusions_lst)):
    plt.plot(x_range_nm,np.sum(cyt_dict['bc1'][i].specR,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)


plt.title('cytochrome bc1')

plt.subplot(2,2,3)
for i in range(len(exclusions_lst)):
    plt.plot(x_range_nm,np.sum(cyt_dict['b6f'][i].specD,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)



plt.subplot(2,2,4)
for i in range(len(exclusions_lst)):
    plt.plot(x_range_nm,np.sum(cyt_dict['bc1'][i].specD,axis = 0),linewidth=2)
#plt.plot(x_range_nm,np.sum(specR_full,axis = 0),linewidth=5)
#handles = []

#plt.subplots_adjust(left=0.07, right=0.93, wspace=0.25, hspace=0.35)
#fig.suptitle('Comparison of CD and OD for two different cytochromes',size=16)

plt.show()
plt.savefig('pairs-cytochrome-b6f-bc1.jpg')
print 'execution complete'



