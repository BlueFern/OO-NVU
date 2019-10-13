# -*- coding: utf-8 -*-
"""
Plotting functions 
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import parameters as p

'''Get suitable subplot dimensions based on the number of variables to be plotted'''
# Returns row & col, the number of rows and columns for the figure
# num_vars: number of variables to be plotted
def get_subplot_dims(num_vars):
    if np.sqrt(num_vars) > 5:
        col = 6                         # set maximum number of columns to 6
        row = np.ceil(num_vars/col)
    elif num_vars == 3:                 # if 3 subplots then plot all in a row
        col = 3
        row = 1
    elif num_vars == 4:                 # if 4 subplots then plot 2x2
        col = 2
        row = 2
    elif num_vars == 8:                 # if 8 subplots then plot 2x4
        col = 4
        row = 2
    else:                           
        col = np.ceil(np.sqrt(num_vars))
        row = np.ceil(num_vars/col)
    return row, col

'''Set a figure size dependent on the dimensions of the plot'''
# Returns width & height of the figure
# row, col: number of rows and columns in the figure
# single_fig_width: width of a 1x1 figure [optional, default is 8]
def get_figure_size(row, col, single_fig_width=8):
    
    if col > 3:
        width = 2*single_fig_width
        height = 4/3*single_fig_width*row/col
    else:
        width = single_fig_width
        height = 2/3*single_fig_width
    
    return width, height

'''Custom function for plotting a selection of variables on separate plots (or a single variable on one plot)'''
# Returns fig, the figure object
# t: time vector
# a: class containing the algebraic variables
# v: class containing the state variables
# xlim1 & xlim2: xlim limits
# plot_list: list containing the names of variables to be plotted

# alphabetical: 1 if want plots done alphabetically (state variables first then algebraic), 0 otherwise go by the list order [optional, default is 0]
# title: title of the whole figure [optional, default is blank]
# figsize = custom figure size if wanted [optional, default is calculated based on number of subplots]
# dpi: custom dpi [optional, default is 80dpi]
# num: figure number, set this if you want to keep adding to an existing figure rather than overwriting when running multiple simulations [optional, default is None]
# show_stimulation: 1 to show a grey box where the stimulation is, 0 to turn off [optional, 0:default]
def plot_variables_singles(t, a, v, xlim1, xlim2, plot_list, alphabetical=0, title='', figsize=(0,0), dpi=80, num = None, show_stimulation=0):
    
    # get the index corresponding to xlim1
    xlim1_idx = np.where(t == xlim1) 
    xlim1_idx = int(xlim1_idx[0][0])
     
    num_vars = len(plot_list)                               # calculate the number of variables to plot 
    row, col = get_subplot_dims(num_vars)                   # get the number of rows and columns for subplot
    
    if figsize == (0,0):
        fig_size = get_figure_size(row,col)                     # get the figure size based on rows and columns
    else:
        fig_size = figsize
        
    fig = plt.figure(figsize=fig_size, dpi=dpi, edgecolor='k', num=num)   # generate the main figure
    plt.tight_layout()
          
    i = 1
    # loop through and plot the state variables, checking if the variable is in plot_list
    loop = sorted(vars(v).items())  # order alphabetically
    for attr, value in loop:
        if attr in plot_list:
            
            if num_vars > 1: # if multiple plots then add subplot
                if alphabetical == 1:   # set position according to order of sorted v
                    position = i
                else:
                    position = plot_list.index(attr) + 1 # find position of variable in plot_list and set the subplot position accordingly (with +1 as indexing begins at zero) 
                ax = fig.add_subplot(row,col,position)
                ax.set_title(attr)
            else:
                ax = fig.gca()
                variable_name = attr
                
            ax.plot(t[xlim1_idx:], value[xlim1_idx:])               
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)   # stop weird exponential notation
            plt.xlabel('Time [s]')
            plt.xlim(xlim1,xlim2)
            plt.tight_layout()
            if abs((max(value[xlim1_idx:]) - min(value[xlim1_idx:]))) < 1e-6: # check if the array is pretty much constant and change ylim to just below & above the mean if so (otherwise gives a misleading figure)
                plt.ylim( 0.9 * np.mean(value[xlim1_idx:]), 1.1 * np.mean(value[xlim1_idx:]) )
            if abs(max(value[xlim1_idx:])) < 1e-6 and abs(min(value[xlim1_idx:])) < 1e-6 and abs(np.mean(value[xlim1_idx:])) < 1e-6: # check if the array is zero for all time and if so set ylim to -1,1
                plt.ylim(-1,1)
            if show_stimulation == 1:
                ymin, ymax = plt.gca().get_ylim()
                plt.gca().add_patch(Rectangle((0, ymin), p.lengthpulse/1e3, ymax-ymin, facecolor='k', alpha=0.1)) # plot an optional grey square where the stimulation is
            i = i + 1
            
    # loop through and plot the algebraic variables, checking if the variable is in plot_list  
    loop = sorted(vars(a).items()) # order alphabetically
    for attr, value in loop:
        if attr in plot_list:
            
            if num_vars > 1: # if multiple plots then add subplot
                if alphabetical == 1:   # set position according to order of sorted a 
                    position = i
                else:
                    position = plot_list.index(attr) + 1 # find position of variable in the plot_list list and set the subplot position accordingly (with +1 as indexing begins at zero)
                ax = fig.add_subplot(row,col,position)
                ax.set_title(attr)
            else:
                ax = fig.gca()  
                variable_name = attr
                                
            if isinstance(value, np.ndarray):   # check if the algebraic variable is an array
                ax.plot(t[xlim1_idx:], value[xlim1_idx:])
                if abs((max(value[xlim1_idx:]) - min(value[xlim1_idx:]))) < 1e-6: # check if the array is pretty much constant and change ylim to below & above the mean if so (otherwise gives a misleading figure)
                    plt.ylim( 0.9 * np.mean(value[xlim1_idx:]), 1.1 * np.mean(value[xlim1_idx:]) )
                if abs(max(value[xlim1_idx:])) < 1e-6 and abs(min(value[xlim1_idx:])) < 1e-6 and abs(np.mean(value[xlim1_idx:])) < 1e-6: # check if the array is zero for all time and if so set ylim to -1,1
                    plt.ylim(-1,1)
            else:                               # otherwise it's a constant (not array, must be handled differently)  
                ax.axhline(y=value)                     
                if abs(value) < 1e-6:
                    plt.ylim(-1,1)  # if value is zero then set ylim to -1,1
                else:
                    plt.ylim(0.9*value, 1.1*value) # otherwise set ylim to just below & above value
                    
            if show_stimulation == 1:
                ymin, ymax = plt.gca().get_ylim()
                plt.gca().add_patch(Rectangle((0, ymin), p.lengthpulse/1e3, ymax-ymin, facecolor='k', alpha=0.1)) # plot an optional grey square where the stimulation is
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
            plt.xlabel('Time [s]')
            plt.xlim(xlim1,xlim2)
            plt.tight_layout()
            i = i + 1

    if not title and num_vars == 1:             # if there is no user specified title and only one variable then make the title the name of the variable
        plt.title(variable_name, y=1)
    elif not title and num_vars > 1:            # if there is no user specified title but more than one variable then print with no title
        pass
    else:
       fig.suptitle(title, y=1)              # otherwise use the user specified title
       
    fig.subplots_adjust(hspace = 1)
    plt.show() 
    return fig

'''Custom function for plotting a selection of variables on the same graph'''
# Returns fig, the figure object
# t: time vector 
# a: class containing the algebraic variables
# v: class containing the state variables
# xlim1 & xlim2: xlim limits
# plot_list: list containing the names of variables to be plotted

# title: title of the whole figure [optional, default is blank]
# legend: custom legend if wanted, must be in order of the variables in plot_list [optional, default uses the standard variable names]
# figsize = custom figure size if wanted [optional, default is calculated based on number of subplots]
# dpi: custom dpi [optional, default is 80dpi]
# num: figure number, set this if you want to keep adding to an existing figure rather than overwriting when running multiple simulations [optional, default is None]
# show_stimulation: 1 to show a grey box where the stimulation is, 0 to turn off [optional, 0:default]
# colourstyle: specify the colour and style of the lines, e.g. ['k-','bo'] [optional, if omitted then just uses the default colour settings]
#   for a list of colour/style options see https://matplotlib.org/2.1.1/api/_as_gen/matplotlib.pyplot.plot.html
def plot_variables_samegraph(t, a, v, xlim1, xlim2, plot_list, title='',legend=[], figsize=(0,0), dpi=80, num = None, show_stimulation=0, colourstyle=[]):
    
    # get the index corresponding to xlim1
    xlim1_idx = np.where(t == xlim1) 
    xlim1_idx = int(xlim1_idx[0][0])
    
    num_vars = len(plot_list)                                   # calculate the number of variables to plot 
    
    if figsize == (0,0):
        fig_size = get_figure_size(1,1)                     # get the figure size based on rows and columns
    else:
        fig_size = figsize
     
    fig = plt.figure(figsize=fig_size, dpi=dpi, edgecolor='k', num=num)   # generate the main figure
    
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)   # stop weird exponential notation
    plt.xlabel('Time [s]')
    plt.xlim(xlim1,xlim2)   
    
    # loop through and plot the state variables, checking if the variable is in plot_list   
    loop = sorted(vars(v).items())  # order alphabetically
    for attr, value in loop:
        if attr in plot_list:
            
            if colourstyle == []:    # if no custom colour/style is specified then use the default
                style = ''
            else:
                style = colourstyle[plot_list.index(attr)]
                
            variable_name = attr
            if not legend:   # if there is no user specified legend then use the names of each variable   
                label = attr
            else:                   # otherwise use the entries from legend            
                label = legend[plot_list.index(attr)]
            plt.plot(t[xlim1_idx:], value[xlim1_idx:], style, label=label)

    # loop through and plot the algebraic variables, checking if the variable is in plot_list  
    loop = sorted(vars(a).items()) # order alphabetically
    for attr, value in loop:
        if attr in plot_list:

            if colourstyle == []:    # if no custom colour/style is specified then use the default
                style = ''
            else:
                style = colourstyle[plot_list.index(attr)]
                
            variable_name = attr
            if not legend:   # if there is no user specified legend then use the labels of each variable   
                label = attr
            else:                   # otherwise use the entries from legend  
                label = legend[plot_list.index(attr)]
                
            if isinstance(value, np.ndarray):   # check if the algebraic variable is an array
                plt.plot(t[xlim1_idx:], value[xlim1_idx:], style, label=label)
            else:
                plt.axhline(style, y=value, label=label)             # otherwise it's a constant (not array, must be handled differently)           
    
    if not title and num_vars == 1:             # if there is no user specified title and only one variable then make the title the name of the variable
        plt.title(variable_name)
    elif not title and num_vars > 1:            # if there is no user specified title but more than one variable then print with no title
        pass
    else:
       plt.title(title)              # otherwise use the user specified title
    
    if show_stimulation == 1:
        ymin, ymax = plt.gca().get_ylim()
        plt.gca().add_patch(Rectangle((0, ymin), p.lengthpulse/1e3, ymax-ymin, facecolor='k', alpha=0.1)) # plot an optional grey square where the stimulation is
    
    plt.legend(loc='best') 
    fig.subplots_adjust(hspace = 1)
    plt.tight_layout()
    plt.show() 
    return fig