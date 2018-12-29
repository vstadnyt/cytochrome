"""

#that should be it. YOu might need xlrd (pip install xlrd).
this is to read excel files but I couldn't make it work. I am doing something
wrong.


You will need to run 'run_csv()' functiop. It will find all .csv files in the
directory. You will need to get rid of everything but leave header, time and dA
columns, see files in the folder.

It will be more work to have this code to read excel files.


Date: July 10, 2018
Authors: By Jullian Ness, Valentyn Stadnytskyi

"""
__version__ = '1.1'



from numpy import loadtxt, transpose,argmax
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



def excel_to_csv():
    """
    converts excel to csv files.
    NOTE: DOESN'T WORK
    """
    import xlrd
    import csv
    import os
    current_folder = os.getcwd()
    lst_dir = os.listdir(current_folder)
    lst_xlsx = []
    for item in lst_dir:
        filename, file_extension = os.path.splitext(current_folder +'/'+ item)
        if file_extension == '.xlsx':
            lst_xlsx.append(item)
            wb = xlrd.open_workbook(filename+file_extension)
            sh = wb.sheet_by_name('Data')
            your_csv_file = open(filename, '.csv', 'w')
            wr = csv.writer(your_csv_file, quoting=csv.QUOTE_ALL)

            for rownum in range(sh.nrows):
                wr.writerow(sh.row_values(rownum))

            your_csv_file.close()

def get_list_csv(folder = ''):
    '''
     gets list of .csv file in the directory with this python code
    '''
    import os
    current_folder = os.getcwd()
    lst_dir = os.listdir(current_folder+'')
    lst_csv = []
    for item in lst_dir:
        filename, file_extension = os.path.splitext(current_folder +'/'+ item)
        if file_extension == '.csv':
            lst_csv.append(item)
    return lst_csv
            
        
def sin_fixed(x,A,t0,y0):
    """
    definition of sin function with fixed freqency
    where
    x - input x vector,
    A - a mplitude,
    t0 - phase,
    y0 - vertical displacement
    frequecny fixed to 60 Hz
    """
    from numpy import sin, pi
    return A*sin(pi*2*60*x+t0) + y0

def sin_free(x,A,f,t0,y0):
    """
    definition of sin function
    """
    from numpy import sin, pi
    return A*sin(pi*2*f*x+t0) + y0

def process(filename = '', plot = False):
    """
    - opens file with filename
    - runs fitting
    - creates a file in subdirectory named \processed with 2 extra columns:
        sin and corrected data
    """
    from numpy import loadtxt, savetxt, concatenate, zeros, mean
    data = loadtxt(filename, delimiter  = ',', skiprows = 1) # import data
    result_data = zeros((data.shape[0],data.shape[1]+2)) #create output array   

    y = data[:,1] #difference
    x = data[:,0] #time
    x_max = argmax(y) # maximum value which is flash lamp peak
    x_fit_to = x_max - 20 #fitted x range from 0 to x_fit_to
    initial_parameters = [max(y[:x_fit_to]),0.01,mean(y[:x_fit_to])] #initial parameters for
                              #A = max value in the region
                              #  phase s10 ms,
                              #vertical displacement = mean value in the region
    popt,pcov = curve_fit(sin_fixed,x[0:x_fit_to],y[0:x_fit_to],p0 = initial_parameters)
    y_sin = sin_fixed(x,*popt)
    y_corrected = y - y_sin
    result_data[:,:2] = data
    result_data[:,2] = y_sin
    result_data[:,3] = y_corrected

    import os
    if not os.path.exists('processed'):
        os.makedirs('processed')
                        
                        
    savetxt('processed/processed_' + filename, result_data, delimiter = ',')
    if plot:
        plt.figure()

        plt.plot(x,y, label = 'raw data')

        plt.plot(x,sin_fixed(x,*popt), linewidth = 1,label = "background")
        plt.plot(x,y-sin_fixed(x,*popt), linewidth = 2, label = 'Corrected data')
        plt.legend()
        plt.show()
    
    
def run_csv():
    """
    looks for .csv files makes a list
    goes through each entry in the list and processes it
    """
    lst = get_list_csv()
    for filename in lst:
        print('processing: %r' % filename)
        process(filename)

def plot():
    legend = ['time', '556nm' , '540nm','556nm-540nm']
    data = transpose(loadtxt('TRIAL 3 wt 0620.csv', delimiter  = ','))
    y = data[1,:] #dA data
    x = data[0,:] #time
    x_max = argmax(y)

    #fitting of sin function with initial parameters p0
    popt,pcov = curve_fit(sin_fixed,x[0:x_max-20],y[0:x_max-20],p0 = [0.0002,0.01,0])

    plt.figure(1)

    plt.plot(x,y, label = 'raw data')

    plt.plot(x,sin_fixed(x,*popt), linewidth = 1,label = "background")
    plt.plot(x,y-sin_fixed(x,*popt), linewidth = 2, label = 'Corrected data')
    plt.legend()

            

    plt.title('freq of noise oscillations fixed to 60Hz=')
                    
    print(popt)
    plt.show()

    """
    #exporting data to .csv file

    output_file = open('output.csv', 'w')#opens output.txt



    output_file.write()

    """

if __name__ == '__main__':
    print('run_csv() #just run me :)')
