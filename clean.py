#f = open('data_bwaves.txt','r');
#lines = f.readlines().strip().split(',');
#boards = [f.readline().strip().split(' ') for i in range(T)]



#import pandas as pd

#data = pd.read_csv('data_bwaves.txt')
#print data

import numpy
import matplotlib.pyplot as plt
#filename = 'data_bwaves.txt'

files = ['data_hmmer.txt' , 'data_mcf.txt', 'data_GemsFDTD.txt',
         'data_lbm.txt'   , 'data_milc.txt', 'data_bwaves.txt',
         'data_leslie3d.txt',  'data_stream.txt', 'data_gups.txt',
         'data_libquantum.txt', 'data_zeusmp.txt'];
for j in range(1,4):
    plt.figure(str(j));
    for i in range(0,11):
        filename = files[i];
        plt.subplot(str(341+i));
        my_data = numpy.genfromtxt(filename, delimiter=',',skip_header=1)
        plt.plot(my_data[:,j]/max(my_data[:,j]))
        plt.title(filename);
        


