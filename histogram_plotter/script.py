#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 13:04:40 2017

@author: bkmiller
"""

import os
import functions as my
import numpy as np
import subprocess
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
import datetime

t_stop = 800;
jstart = 75;
jstop = 150;
jincr = 5;
fA_ss = 0.2;
fA_ww = 0;

recalculate = 1;
if recalculate != 0:
	#change to
	#~ os.chdir('/home/bkmiller/stochastics/')
	#~ 
	#~ call_program = './muPillar_FB' 
	
	call_program =  '/home/bkmiller/stochastics/muPillar_FB' 
	call_program += ' -path ' + '/home/bkmiller/stochastics/data/'
	call_program += ' -srccode ' + ' /home/bkmiller/stochastics/muPillar_FB.cc'
	call_program += ' -sponE 1'
	call_program += ' -t_stop ' + str(t_stop)
	call_program += ' -jstop ' + str(jstop)
	call_program += ' -jstart ' + str(jstart)
	call_program += ' -jincr '  + str(jincr)
	call_program += ' -fA_ss '  + str(fA_ss)
	call_program += ' -fA_ww '  + str(fA_ww) 
	
	print(call_program)
	
	subprocess.call(call_program, shell=True)

# change to
#~ os.chdir('/home/bkmiller/stochastics/data/muPillar_FB/2017-5-23')
datadir = '/home/bkmiller/stochastics/data/muPillar_FB/2017-5-24/'

#~ current = '160'
#~ 
 #~ _, weak = my.read_data(str(current)+'.00_muA_w_ts.data')
 #~ _, strong = my.read_data(str(current)+'.00_muA_s_ts.data')
#~ 
 #~ my.histplot(strong,weak)

# Now all of em
cur = np.arange(jstart,jstop,jincr)

if not os.path.exists('/home/bkmiller/stochastics/histogram_plotter/plots'):
    os.makedirs('/home/bkmiller/stochastics/histogram_plotter/plots')
if not os.path.exists('/home/bkmiller/stochastics/histogram_plotter/plots/' + str(datetime.date.today()) + '/'):
    os.makedirs('/home/bkmiller/stochastics/histogram_plotter/plots/' + str(datetime.date.today()) + '/')
todaydir = '/home/bkmiller/stochastics/histogram_plotter/plots/' + str(datetime.date.today()) + '/'

for ind,val in enumerate(cur):
	
	if len(str(val)) < 3:
		val = '0'+str(val)
	
	_, weak = my.read_data(datadir + str(val)+'.00_muA_w_ts.data')
	_, strong = my.read_data(datadir + str(val)+'.00_muA_s_ts.data')
	
	#~ Let's create the plot and save it
	print('current: '+str(val))
	
	plt.figure(figsize=(11.69,8.27))
	plt.scatter(strong, weak, s=3.1415926*0.5**2, c = 'red', alpha=0.5)
	plt.gca().set_xlim(left=0)
	plt.gca().set_ylim(bottom=0)
	plt.title('Current: '+str(val))
	
	print(todaydir + str(val) + 'uA_histogram.png')
	plt.savefig(todaydir + str(val) + 'uA_histogram.png', bbox_inches='tight')
	
	plt.close()
	
	#~ Creates files which are way too big
	#~ with PdfPages('histogram.pdf') as pdf:
		#~ plt.figure(figsize=(11.69,8.27))
		#~ plt.scatter(strong, weak, s=3.1415926*0.5**2, c = 'red', alpha=0.5)
		#~ plt.gca().set_xlim(left=0)
		#~ plt.gca().set_ylim(bottom=0)
		#~ plt.title('Current: '+str(val))
		#~ 
		#~ pdf.savefig()
		#~ plt.close
