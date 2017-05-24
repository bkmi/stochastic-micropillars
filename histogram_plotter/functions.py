#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 12:08:26 2017

@author: bkmiller
"""

import csv
import matplotlib.pyplot as plt


def read_data(filename):
    """
    This function loads the data from a single file. It is designed to work 
    with time series data for the histogram.
    """
    
    time = []    
    intensity = []
    
    with open(filename) as tsv:
        next(tsv) # skip head
        reader = csv.reader(tsv,delimiter='\t')
        for t, i in reader:
            time.append(t)
            intensity.append(i)
            
    return time, intensity


def histplot(x, y):
    """
    Input the xdir intensity first and the ydir intensity second.
    """
    
    fig = plt.scatter(x, y, s=3.1415926*0.5**2, c = 'red', alpha=0.5)
    plt.gca().set_xlim(left=0)
    plt.gca().set_ylim(bottom=0)
    
    return fig
