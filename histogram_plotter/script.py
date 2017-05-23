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

t_stop = 800;
jstart = 75;
jstop = 150;
jincr = 5;
fA_ss = 0.2;
fA_ww = 0;

recalculate = 1;
if recalculate != 0:
    #change to
    os.chdir('/home/bkmiller/stochastics/')
    
    call_program = './muPillar_FB' 
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
os.chdir('/home/bkmiller/stochastics/data/muPillar_FB/2017-5-23')


# current = '160'
#
# _, weak = my.read_data(str(current)+'.00_muA_w_ts.data')
# _, strong = my.read_data(str(current)+'.00_muA_s_ts.data')
#
# my.histplot(strong,weak)

# Now all of em

cur = np.arange(jstart,jstop,jincr)

for ind,val in enumerate(cur):
    
    if len(str(val)) < 3:
        val = '0'+str(val)
        
    _, weak = my.read_data(str(val)+'.00_muA_w_ts.data')
    _, strong = my.read_data(str(val)+'.00_muA_s_ts.data')
    
    print('current: '+str(val))
    my.histplot(strong,weak)