import numpy as np
import matplotlib.pyplot as plt
import os
import re

##This code is a part of researcher published in two papers:
##Paper1: A Map of Dielectric Heterogeneity in a Membrane Protein: the Hetero-Oligomeric Cytochrome b6f Complex S. Saif Hasan, Stanislav D. Zakharov, Adrien Chauvet, Valentyn Stadnytskyi, Sergei Savikhin, and William A. Cramer The Journal of Physical Chemistry B 2014 118 (24), 6614-6625 DOI: 10.1021/jp501165k
##Paper2: Pathways of Transmembrane Electron Transfer in Cytochrome bc Complexes: Dielectric Heterogeneity and Interheme Coulombic Interactions S. Bhaduri, V. Stadnytskyi, S. D. Zakharov, S. Saif Hasan, L. Bujnowicz, M. Sarewicz, S. Savikhin, A. Osyczka, and W. A. Cramer The Journal of Physical Chemistry B 2017 121 (5), 975-983 DOI: 10.1021/acs.jpcb.6b11709
## The author of this library is Valentyn Stadnytskyi, v.stadnytskyi@gmail.com. This library was produced as a part of Ph.D. thesis with Prof. Savikhin sergei@purdue.edu and Prof. Cramer waclab@purdue.edu at Purdue University

####################################
##IMPORTANT#########################
####################################
#This library was written and tested on:
#Python 2.7.12 |Anaconda 4.1.1 (64-bit)| (default, Jun 29 2016, 11:07:13) [MSC v.1500 64 bit (AMD64)]
####################################
str_success = 'The Cytochrome Library has been succesfully imported'

str_version = '0.1alpha'
str_date = 'August 17, 2017'
def read_config_file(filename,method = None):
    #filename = 'config.cyt'
    return_dict = {}
    with open (filename, "r") as myfile:
        data = myfile.read()
    if method == None:
        config_dict = {}
        config_dict['author name'] = re.search('author:(.*)\n', data).group(1)
        config_dict['author email'] = re.search('email:(.*)\n', data).group(1).replace(' ', '')
        config_dict['version'] = re.search('version:(.*)\n', data).group(1).replace(' ', '')
        config_dict['date'] = re.search('date:(.*)\n', data).group(1).replace(' ', '')
        config_dict['input folder'] = re.search('input folder:(.*)\n', data).group(1).replace(' ', '').replace('\n', '')
        return_dict = config_dict
    if method == 'commands':
        command_dict = {}
        command_dict['help'] = (data.split('<help>\n'))[1].split('\n<help>')[0]
        command_dict['exit'] =(data.split('<exit>\n'))[1].split('\n<exit>')[0]
        command_dict['other'] =(data.split('<other>\n'))[1].split('\n<other>')[0]
        return_dict = command_dict
    return return_dict
        
class cyt(object):
     
    def __init__(self, filename, exceptions):
        for i in exceptions:
            if type(i) !=  int:
                raise ValueError("values in the exception list have to be an integer but ",str(type(i))," is found")
            if i != 1 and i != 0:
                raise ValueError("all values in the exception list have to be 1 or 0." + "Values is " + str(i))
        if os.path.exists(filename):
            pass
        else:
            raise ValueError("the file with name " + filename + " doesn't exist or the path to the file is wrong.")
            
        """Return a ochrome object whose structural data comes fomr file with filename and has 
        populated hemes defined by exceptions"""
        self.filename = filename
        self.exceptions = exceptions

    def read_structure_file(self):
        
        with open(self.filename, "r") as myfile:
            self.data=myfile.readlines()
        self.wavelength = float(self.data[0].partition('wavelength(nm):')[2])

        self.width = float(self.data[1].partition('width(nm):')[2])

        bp1 = np.zeros((5,3))
        bn1 = np.zeros((5,3))
        bp2 = np.zeros((5,3))
        bn2 = np.zeros((5,3))
        for i in range(5):
            coord = self.data[5+i].partition(':')[2]
            bp1[i,:] = np.asarray([float(x) for x in coord.split(',')])
        for i in range(5):
            coord = self.data[11+i].partition(':')[2]
            bn1[i,:] = np.asarray([float(x) for x in coord.split(',')])
        for i in range(5):
            coord = self.data[17+i].partition(':')[2]
            bp2[i,:] = np.asarray([float(x) for x in coord.split(',')])
        for i in range(5):
            coord = self.data[23+i].partition(':')[2]
            bn2[i,:] = np.asarray([float(x) for x in coord.split(',')])
        self.coord = np.vstack([bp1,bn1,bp2,bn2])
        return self.wavelength , self.width, self.coord


    def Hamiltonian(self):
        
        self.Ham = np.zeros((8,8))

        self.transition_dipole_moments = np.zeros((8,3))
        self.Fe_coordinates = np.zeros((8,3))

        for i in range(4):
            self.transition_dipole_moments[2*i,:] = (self.coord[3+5*i,:] - self.coord[1+5*i,:]) / np.linalg.norm(self.coord[3+5*i,:] - self.coord[1+5*i,:])
            self.transition_dipole_moments[2*i+1,:] = (self.coord[2+5*i,:] - self.coord[0+5*i,:]) / np.linalg.norm(self.coord[2+5*i,:] - self.coord[0+5*i,:])
            self.Fe_coordinates[2*i,:] = self.coord[4+5*i,:]
            self.Fe_coordinates[2*i+1,:] = self.coord[4+5*i,:]
        for k in range(8):
            for kk in range(8):
                mu_i = self.transition_dipole_moments[k,:];
                mu_j = self.transition_dipole_moments[kk,:];
                r = self.Fe_coordinates[k,:] - self.Fe_coordinates[kk,:]
                r_mag = np.linalg.norm(r);
                if r_mag == 0:
                     self.Ham[k,kk] = 0;
                else:
                    r_hat = r/r_mag
                    self.Ham[k,kk] = 4.6*4.6*5040*(np.dot(mu_i,mu_j) - 3*np.dot(mu_i,r_hat)*np.dot(mu_j,r_hat)) / (r_mag**3)
        for i in range(8):
            self.Ham[i,i] = 10**7/self.wavelength
        exceptions_matrix = np.outer(self.exceptions,self.exceptions)
        self.Ham = np.multiply(self.Ham,exceptions_matrix)
        self.Eigenvalues, self.Eigenvectors = np.linalg.eig(self.Ham)
        return self.Ham, self.Eigenvalues, self.Eigenvectors
    
    def D_and_R_strength(self):
        #here we will create empty numpy arrays of defined size
        self.DipoleS  = np.zeros(8)
        self.RotationalS  = np.zeros(8)
        self.Density_matrix = np.zeros((8,8,8))
        self.Dipole_strength_matrix = np.zeros((8,8))
        self.Rotational_strength_matrix = np.zeros((8,8))
    
        for e in range(8):
            for i in range(8):
                for j in range(8):
                    self.Density_matrix[i,j,e] = self.Eigenvectors[i,e]*self.Eigenvectors[j,e]
        for i in range(8):
            for j in range(8):
                self.Dipole_strength_matrix[i,j] = np.dot(self.transition_dipole_moments[i],self.transition_dipole_moments[j])       
        for i in range(8):
            for j in range(8):
                self.Rotational_strength_matrix[i,j] = np.dot((self.Fe_coordinates[i,:]-self.Fe_coordinates[j,:]),np.cross(self.transition_dipole_moments[i],self.transition_dipole_moments[j]))  
        #vec[:,i] corresponds to value[i]
        for e in range(8):
            for i in range(8):
                for j in range(8):
                    self.DipoleS[e] = self.DipoleS[e] + self.Density_matrix[i,j,e] * self.Dipole_strength_matrix[i,j] 
                    self.RotationalS[e] = self.RotationalS[e] + self.Density_matrix[i,j,e] * self.Rotational_strength_matrix[i,j]
        return self.DipoleS, self.RotationalS
    
    def spectra_plot(self):
        
        self.x_range_cm = np.linspace(22500, 23700, 400)
        self.x_range_nm = range(400)
        for i in range(len(self.x_range_cm)):
            self.x_range_nm[i] = 10**7/self.x_range_cm[i]
    
        self.specD = np.zeros((8,len(self.x_range_cm)))
        self.specR = np.zeros((8,len(self.x_range_cm)))
        for i in range(8):
            self.specD[i,:] = self.DipoleS[i]*np.exp(-np.power(self.x_range_cm - self.Eigenvalues[i], 2.) / (2 * np.power(221.38, 2.)))
            self.specR[i,:] = 10**-7*1.57*np.multiply(self.RotationalS[i]*np.exp(-np.power(self.x_range_cm - self.Eigenvalues[i], 2.) / (2 * np.power(221.38, 2.))),self.x_range_cm)
        return self.specD, self.specR, self.x_range_nm

    
def kinetics_solve(k, length):
    import numpy as np
    import matplotlib.pyplot as plt
    """Return the balance remaining after withdrawing *amount*
    dollars."""
    #k = np.array([1,1,1,1,1,1,1])
    #length = 10000
    k12 = k[0];
    k23 = k[1];
    k34 = k[2];
    k45 = k[3];
    k13 = k[4];
    k35 = k[5];
    k24 = k[6];


    N = np.zeros((5,length))
    ##Initial parameters
    N[0,0] = 1
    ##Solving differential equations numerically
    for i in range(length-1):
        dt = 1./(length*0.1) 
        N[0,i+1] = N[0,i] - k12*N[0,i]*dt - k13*N[0,i]*dt;
        N[1,i+1] = N[1,i] + k12*N[0,i]*dt - k23*N[1,i]*dt - k24*N[1,i]*dt;
        N[2,i+1] = N[2,i] + k23*N[1,i]*dt - k34*N[2,i]*dt + k13*N[0,i]*dt - k35*N[2,i]*dt;
        N[3,i+1] = N[3,i] + k34*N[2,i]*dt - k45*N[3,i]*dt + k24*N[1,i]*dt;
        N[4,i+1] = N[4,i] + k45*N[3,i]*dt + k35*N[2,i]*dt;
    return N

print str_success,'--->', 'Current version:',str_version, '---> last update',str_date
