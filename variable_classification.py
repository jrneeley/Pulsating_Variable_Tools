'''
    File name: variable_classification.py
    Author: Jillian Neeley
    Date last modified: 8/1/2021
    Python Version: 3.6
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
import sys
import glob
import plotting_utilities as putil
from config import *

global star
global results

# Allows you to switch classification modes
run = sys.argv[1]

# Allows you to look at one specific star
if len(sys.argv) == 3:
    star = int(sys.argv[2])
    results_file = data_dir+'fit_results.txt'
    lcv_dir = data_dir+'lcvs/'
    lcv_prefix = 'c'
    flags = np.loadtxt(results_file, usecols=(0), dtype=int)
    stars_left_to_do = np.argwhere(flags[:] == star)

# Use for first run through classification
# (Loops through all stars with visual flag == -1)
elif run == 'first pass':
    results_file = data_dir+'fit_results.txt'
    lcv_dir = data_dir+'lcvs/'
    lcv_prefix = 'c'
    flags = np.loadtxt(results_file, usecols=(13))
    stars_left_to_do = np.argwhere(flags[:] == -1)

# Use to make a second pass at classification
# (Loops through all stars that were not marked as non variable during first pass)
elif run == 'revise':
    results_file = data_dir+'fit_results.txt'
    lcv_dir = data_dir+'lcvs/'
    lcv_prefix = 'c'
    flags = np.loadtxt(results_file, usecols=(14), dtype='U3')
    stars_left_to_do = np.argwhere(flags[:] != 'NV')

# Use to classify eclipsing binary stars
# (Loops through all stars marked as EB)
elif run == 'EB':
    results_file = data_dir+'fit_results.txt'
    lcv_dir = data_dir+'lcvs/'
    lcv_prefix = 'c'
    flags = np.loadtxt(results_file, usecols=(14), dtype='U3')
    stars_left_to_do = np.argwhere(flags[:] == 'EB')

# Use to classify simulated stars
# (Different input file/output directory)
elif run == 'simulations':
    results_file = data_dir+'simulations/sim_results.txt'
    lcv_dir = data_dir+'simulations/lcvs/'
    lcv_prefix = 'sim_'
    flags = np.loadtxt(results_file, usecols=(21,22))
    stars_left_to_do = np.argwhere((flags[:,0] == 1) & (flags[:,1] == -1))

# Exit the program when we've looped through all stars
if len(stars_left_to_do) == 0:
    print('You\'ve already finished this dataset!')
    sys.exit()

next_star = stars_left_to_do[0][0]
iter = 0

# read in full catalog for cmd/image panels
# Required format for catalog:
#   col 2: x coordinate
#   col 3: y coordinate
#   col 4: magnitude in bluest filter
#   col 6: magnitude in reddest filter
#   col 9: sharpness value
# Format of other columns not important!

dt = np.dtype([('x', float), ('y', float), ('mag1', float), ('mag2', float),
    ('sharp', float)])
cat_data = np.loadtxt(catalog, dtype=dt, usecols=(1,2,3,5,8))
# assumes missing magnitudes are 99.999
cat_data['mag1'][cat_data['mag1'] > 90] = np.nan
cat_data['mag2'][cat_data['mag2'] > 90] = np.nan

def plot_star_data(star, fig, cat_data):

    if run != 'simulations':
        dt = np.dtype([('var_id', int), ('dao_id', int), ('x', float), ('y', float),
            ('fit_template', int), ('fit_period', float), ('fit_t0', float),
            ('fit_mag1', float), ('fit_amp1', float), ('fit_sig1', float),
            ('fit_mag2', float), ('fit_amp2', float), ('fit_sig2', float),
            ('pass_flag', int), ('fit_type', 'U3'), ('fit_mode', 'U2')])
    else:
        dt = np.dtype([('var_id', int), ('type', 'U3'), ('mode', 'U2'),
            ('template', int), ('period', float), ('t0', float),
            ('mag1', float), ('amp1', float), ('sig1', float),
            ('mag2', float), ('amp2', float), ('sig2', float),
            ('fit_template', int), ('fit_period', float), ('fit_t0', float),
            ('fit_mag1', float), ('fit_amp1', float), ('fit_sig1', float),
            ('fit_mag2', float), ('fit_amp2', float), ('fit_sig2', float),
            ('pass_index', int), ('pass_flag', int), ('fit_type', 'U3'),
            ('fit_mode', 'U2')])

    results = np.loadtxt(results_file, dtype=dt)
    nstars = len(results['var_id'])

    # exit program if we've reached the end of the file
    if star >= nstars:
        print('Reached end of file!')
        fig.close()
        sys.exit()

    # Setup plots for GUI window
    figure1 = fig
    figure1.clf()
    gs = GridSpec(3,4, figure=figure1)
    ax1 = figure1.add_subplot(gs[0,0]) # LCV1 (bluest)
    ax2 = figure1.add_subplot(gs[0,1]) # LCV2 (reddest)
    ax3 = figure1.add_subplot(gs[0,2]) # CMD
    ax4 = figure1.add_subplot(gs[0,3]) # Image
    ax5 = figure1.add_subplot(gs[1,2]) # period-color
    ax6 = figure1.add_subplot(gs[1,3]) # period-luminosity (reddest filter)
    ax7 = figure1.add_subplot(gs[2,2]) # Amp
    ax8 = figure1.add_subplot(gs[2,3]) # Amp ratio (bluest/reddest)
    axbig1 = figure1.add_subplot(gs[1,0:2]) # LCV1 raw (bluest)
    axbig2 = figure1.add_subplot(gs[2,0:2]) # LCV1 raw (reddest)

    f1 = FigureCanvasTkAgg(figure1, master)
    f1.get_tk_widget().place(relx=0.025, rely=0.025, relwidth=0.75, relheight=0.95)

    # Print some useful header information about current star
    figure1.suptitle('Star {}: Period={:.3f}    {} stars left\n Type: {} Mode: {} Flag: {}'.format(\
        results['var_id'][star], results['fit_period'][star],
        len(stars_left_to_do) - iter, results['fit_type'][star],
        results['fit_mode'][star], results['pass_flag'][star]))


    confirmed = (results['fit_type'] == 'RRL') | (results['fit_type'] == 'CEP') | \
        (results['fit_type'] == 'AC') | (results['fit_type'] == 'T2C')
    left = (results['pass_flag'] == -1)

    # set up some array short cuts
    logP = np.log10(results['fit_period'])
    color = results['fit_mag1'] - results['fit_mag2']
    amp_ratio = results['fit_amp1']/results['fit_amp2']

    # Plot CMD (NOT extinction corrected)
    sel = np.abs(cat_data['sharp']) < 0.1
    ax3.scatter(cat_data['mag1'][sel]-cat_data['mag2'][sel], cat_data['mag2'][sel],
        s=1, alpha=0.05, color='gray')
    ax3.scatter(color[confirmed], results['fit_mag2'][confirmed], s=1, color='black')
    ax3.scatter(color[star], results['fit_mag2'][star], s=5, color='red')
    ax3.invert_yaxis()
    ax3.set_xlim(-0.5, 3.0)
    ax3.set_xlabel(f'{filt1_name} - {filt2_name}')
    ax3.set_ylabel(f'{filt1_name}')

    # Image Plot
    # only plot if we are not doing simulations
    if run != 'simulations':
        putil.plot_region(results['x'][star]-montage_x_offset, results['y'][star]-montage_y_offset, image,
            axes=ax4, xall=cat_data['x'], yall=cat_data['y'], xoff=montage_x_offset, yoff=montage_y_offset,
            aperture=5, img_limits=image_limits)


    # plot PC relation
    # This needs to be edited!! currently only works for B & I
    ax5.scatter(logP[left], color[left] - (ext_mag1 - ext_mag2), s=1,
        color='xkcd:gray', alpha=0.5)
    ax5.scatter(logP[confirmed], color[confirmed] - (ext_mag1 - ext_mag2),
        s=1, color='xkcd:black', alpha=0.5)
    ax5.plot(logP[star], color[star] - (ext_mag1 - ext_mag2), marker='o',
        color='xkcd:red')
    if len(logP[left]) > 0:
        min_logP = np.nanmin(logP)
        max_logP = np.nanmax(logP)
    else:
        min_logP = np.nanmin(logP[confirmed])
        max_logP = np.nanmax(logP[confirmed])
    if plot_pc == True:
        x_cep = np.array([min_logP, max_logP])
        y_cep = cep_pc_relation_coeffs[0] + cep_pc_relation_coeffs[1]*x_cep
        x_rrl = np.array([min_logP, 0.0])
        y_rrl = rrl_pc_relation_coeffs[0] + rrl_pc_relation_coeffs[1]*x_rrl
        ax5.fill_between(x_cep, y_cep-cep_pc_relation_coeffs[2],
            y_cep+cep_pc_relation_coeffs[2], color='xkcd:steel blue', alpha=0.4)
        ax5.fill_between(x_rrl, y_rrl-rrl_pc_relation_coeffs[2],
            y_rrl+rrl_pc_relation_coeffs[2], color='xkcd:pale purple', alpha=0.4)
    ax5.set_xlabel('logP')
    ax5.set_ylabel(f'{filt1_name} - {filt2_name}')

    # Plot PL relation of redder band
    ax6.scatter(logP[left], results['fit_mag2'][left] - ext_mag2, s=1,
        color='xkcd:gray', alpha=0.5)
    ax6.scatter(logP[confirmed], results['fit_mag2'][confirmed] - ext_mag2, s=1,
        color='xkcd:black', alpha=0.5)
    ax6.plot(np.log10(results['fit_period'][star]), results['fit_mag2'][star] - ext_mag2,
        marker='o', color='xkcd:red')
    ax6.invert_yaxis()
    ax6.set_xlabel('logP')
    ax6.set_ylabel(f'{filt2_name}')
    # add on LMC PL relations if desired
    if plot_lmc == True:
        putil.plot_lmc_pl(axes=ax6, offset=DM, period_cutoff=np.max(results['fit_period'])+0.2)

    # Plot Bailey diagram for bluer band
    ax7.scatter(logP[left], results['fit_amp1'][left], s=1,
        color='xkcd:gray', alpha=0.5)
    ax7.scatter(logP[confirmed], results['fit_amp1'][confirmed], s=1,
        color='xkcd:black', alpha=0.5)
    ax7.plot(logP[star], results['fit_amp1'][star],
        marker='o', color='xkcd:red')
    ax7.set_xlabel('logP')
    ax7.set_ylabel(f'{filt1_name} amp')

    # Plot amplitude ratio
    ax8.scatter(logP[left], amp_ratio[left], s=1,
        color='xkcd:gray', alpha=0.5)
    ax8.scatter(logP[confirmed], amp_ratio[confirmed], s=1,
        color='xkcd:black', alpha=0.5)
    ax8.plot(logP[star], amp_ratio[star], marker='o', color='xkcd:red')
    ax8.set_xlabel('logP')
    ax8.set_ylabel(f'amp_{filt1_name}/amp_{filt2_name}')

    # get lcv information
    dt = np.dtype([('filt', int), ('mjd', float), ('mag', float), ('err', float)])
    lcv = np.loadtxt(lcv_dir+lcv_prefix+'{}'.format(results['var_id'][star])+lcv_suffix,
        dtype=dt, skiprows=3, usecols=(0,1,2,3))
    fit = np.loadtxt(lcv_dir+lcv_prefix+'{}'.format(results['var_id'][star])+fit_suffix,
        dtype=([('phase', float), ('mag1', float), ('mag2', float)]), skiprows=1)

    filt = lcv['filt'] == 0
    phase = np.mod((lcv['mjd'][filt]-results['fit_t0'][star])/results['fit_period'][star], 1)
    phase = np.concatenate((phase, phase+1))
    mag = np.tile(lcv['mag'][filt],2)
    err = np.tile(lcv['err'][filt],2)
    phase_fit = np.concatenate((fit['phase'], fit['phase']+1))
    mag_fit = np.tile(fit['mag1'], 2)

    # plot phased light curve for each filter
    ax1.errorbar(phase, mag, yerr=err, fmt='.', color='k')
    ax1.plot(phase_fit, mag_fit, color='xkcd:ocean blue')
    ax1.set_xlabel('Phase')
    ax1.set_ylabel(f'{filt1_name}')
    ax1.set_ylim(np.max(lcv['mag'][filt])+0.3,
        np.min(lcv['mag'][filt])-0.3)

    filt = lcv['filt'] == 1
    phase = np.mod((lcv['mjd'][filt]-results['fit_t0'][star])/results['fit_period'][star], 1)
    phase = np.concatenate((phase, phase+1))
    mag = np.tile(lcv['mag'][filt],2)
    err = np.tile(lcv['err'][filt],2)
    phase_fit = np.concatenate((fit['phase'], fit['phase']+1))
    mag_fit = np.tile(fit['mag2'], 2)

    ax2.errorbar(phase, mag, yerr=err, fmt='.', color='k')
    ax2.plot(phase_fit, mag_fit, color='xkcd:ocean blue')
    ax2.set_xlabel('Phase')
    ax2.set_ylabel(f'{filt2_name}')
    ax2.set_ylim(np.max(lcv['mag'][filt])+0.3,
        np.min(lcv['mag'][filt])-0.3)

    # plot unphased light curves
    mjd_order = np.argsort(lcv['mjd'])
    min_mjd = lcv['mjd'][mjd_order][0] - 0.1
    max_mjd = lcv['mjd'][mjd_order][-1] + 0.1
    mjd_window = max_mjd - min_mjd
    time_diff = np.diff(lcv['mjd'][mjd_order])

    axbig1.set_xlim(min_mjd, max_mjd)
    axbig2.set_xlim(min_mjd, max_mjd)
    axbig1.set_xlabel('MJD')
    axbig2.set_xlabel('MJD')
    axbig1.set_ylabel(f'{filt1_name}')
    axbig2.set_ylabel(f'{filt2_name}')

    filt = lcv['filt'] == 0
    axbig1.errorbar(lcv['mjd'][filt], lcv['mag'][filt], yerr=lcv['err'][filt],
        fmt='.', color='k')
    ncycles = int(np.ceil(mjd_window/results['fit_period'][star]))*2

    tt = fit['phase']*results['fit_period'][star] + results['fit_t0'][star]
    ttt = []
    for j in np.arange(0,ncycles):
        ttt = np.append(ttt, tt+(j-ncycles/2)*results['fit_period'][star])
    mm = np.tile(fit['mag1'], ncycles)
    axbig1.plot(ttt, mm, color='xkcd:ocean blue')
    axbig1.set_ylim(np.mean(fit['mag1'])+1.0, np.mean(fit['mag1'])-1.0)

    filt = lcv['filt'] == 1
    axbig2.errorbar(lcv['mjd'][filt], lcv['mag'][filt], yerr=lcv['err'][filt],
        fmt='.', color='k')
    ncycles = int(np.ceil(mjd_window/results['fit_period'][star]))*2

    tt = fit['phase']*results['fit_period'][star] + results['fit_t0'][star]
    ttt = []
    for j in np.arange(0,ncycles):
        ttt = np.append(ttt, tt+(j-ncycles/2)*results['fit_period'][star])
    mm = np.tile(fit['mag2'], ncycles)
    axbig2.plot(ttt, mm, color='xkcd:ocean blue')
    axbig2.set_ylim(np.mean(fit['mag2'])+1.0, np.mean(fit['mag2'])-1.0)

def var_states():

    # pulsation modes
    if fo.get() == 1:
        mode = 'FO'
    elif fu.get() == 1:
        mode = 'FU'
    # eclipsing binary types
    elif ew.get() == 1:
        mode = 'EW'
    elif ea.get() == 1:
        mode = 'EA'
    elif eb2.get() == 1:
        mode = 'EB'
    else:
        mode = 'XX'

    # determine variable type
    if rrl.get() == 1:
        type = 'RRL'
    elif cep.get() == 1:
        type = 'CEP'
    elif ac.get() == 1:
        type = 'AC'
    elif t2c.get() == 1:
        type = 'T2C'
    elif eb.get() == 1:
        type = 'EB'
    else:
        type = 'XXX'

    return type, mode

def reset_var_states():
    fo.set(0)
    fu.set(0)
    rrl.set(0)
    cep.set(0)
    ac.set(0)
    t2c.set(0)
    eb.set(0)
    ea.set(0)
    eb2.set(0)
    ew.set(0)
    unsure.set(0)

# star is likely a variable star
def pass_star():
    global next_star
    global iter
    if run != 'simulations':
        dt = np.dtype([('var_id', int), ('dao_id', int), ('x', float), ('y', float),
            ('fit_template', int), ('fit_period', float), ('fit_t0', float),
            ('fit_mag1', float), ('fit_amp1', float), ('fit_sig1', float),
            ('fit_mag2', float), ('fit_amp2', float), ('fit_sig2', float),
            ('pass_visual', int), ('fit_type', 'U3'), ('fit_mode', 'U2')])
        fmt1 = '%4d %7d %8.3f %8.3f %1d %9.6f %10.4f '
        fmt2 = '%6.3f %5.3f %4.2f '
        fmt3 = '%2d %3s %2s'
        fmt = fmt1+2*fmt2+fmt3

    else:
        dt = np.dtype([('var_id', int), ('type', 'U3'), ('mode', 'U2'),
            ('template', int), ('period', float), ('t0', float),
            ('mag1', float), ('amp1', float), ('sig1', float),
            ('mag2', float), ('amp2', float), ('sig2', float),
            ('fit_template', int), ('fit_period', float), ('fit_t0', float),
            ('fit_mag1', float), ('fit_amp1', float), ('fit_sig1', float),
            ('fit_mag2', float), ('fit_amp2', float), ('fit_sig2', float),
            ('pass_index', int), ('pass_visual', int), ('fit_type', 'U3'),
            ('fit_mode', 'U2')])
        fmt1 = '%4d %3s %2s %1d %9.6f %10.4f '
        fmt2 = '%6.3f %5.3f %4.2f '
        fmt3 = '%1d %9.6f %10.4f '
        fmt4 = '%1d %1d %3s %2s'
        fmt = fmt1+2*fmt2+fmt3+2*fmt2+fmt4

    results = np.loadtxt(results_file, dtype=dt)
    if unsure.get() == 0:
        results['pass_visual'][next_star] = 1
    if unsure.get() == 1:
        results['pass_visual'][next_star] = 2
    type, mode = var_states()
    results['fit_type'][next_star] = type
    results['fit_mode'][next_star] = mode

    np.savetxt(results_file, results, fmt=fmt)
    iter += 1
    if iter == len(stars_left_to_do):
        print('Finished!')
        sys.exit()
    next_star = stars_left_to_do[iter][0]
    reset_var_states()
    plot_star_data(next_star, figure1, cat_data)

# star does not pass visual inspection
def fail_star():
    global next_star
    global iter
    if run != 'simulations':
        dt = np.dtype([('var_id', int), ('dao_id', int), ('x', float), ('y', float),
            ('fit_template', int), ('fit_period', float), ('fit_t0', float),
            ('fit_mag1', float), ('fit_amp1', float), ('fit_sig1', float),
            ('fit_mag2', float), ('fit_amp2', float), ('fit_sig2', float),
            ('pass_visual', int), ('fit_type', 'U3'), ('fit_mode', 'U2')])
        fmt1 = '%4d %7d %8.3f %8.3f %1d %9.6f %10.4f '
        fmt2 = '%6.3f %5.3f %4.2f '
        fmt3 = '%2d %3s %2s'
        fmt = fmt1+2*fmt2+fmt3

    else:
        dt = np.dtype([('var_id', int), ('type', 'U3'), ('mode', 'U2'),
            ('template', int), ('period', float), ('t0', float),
            ('mag1', float), ('amp1', float), ('sig1', float),
            ('mag2', float), ('amp2', float), ('sig2', float),
            ('fit_template', int), ('fit_period', float), ('fit_t0', float),
            ('fit_mag1', float), ('fit_amp1', float), ('fit_sig1', float),
            ('fit_mag2', float), ('fit_amp2', float), ('fit_sig2', float),
            ('pass_index', int), ('pass_visual', int), ('fit_type', 'U3'),
            ('fit_mode', 'U2')])
        fmt1 = '%4d %3s %2s %1d %9.6f %10.4f '
        fmt2 = '%6.3f %5.3f %4.2f '
        fmt3 = '%1d %9.6f %10.4f '
        fmt4 = '%1d %1d %3s %2s'
        fmt = fmt1+2*fmt2+fmt3+2*fmt2+fmt4

    results = np.loadtxt(results_file, dtype=dt)
    results['pass_visual'][next_star] = 0
    results['fit_type'][next_star] = 'NV'
    results['fit_mode'][next_star] = 'NV'

    np.savetxt(results_file, results, fmt=fmt)
    iter += 1
    if iter == len(stars_left_to_do):
        print('Finished!')
        sys.exit()
    next_star = stars_left_to_do[iter][0]
    reset_var_states()
    plot_star_data(next_star, figure1, cat_data)

# go back to the star you just did
def go_back():
    global next_star
    global iter
    iter -= 1
    next_star = stars_left_to_do[iter][0]
    reset_var_states()
    plot_star_data(next_star, figure1, cat_data)

# move on without changing file
def move_on():
    global next_star
    global iter
    iter += 1
    if iter == len(stars_left_to_do):
        print('Finished!')
        sys.exit()
    next_star = stars_left_to_do[iter][0]
    reset_var_states()
    plot_star_data(next_star, figure1, cat_data)
    #print('Star index {}'.format(star))

def lpv_star():
    global next_star
    global iter
    if run != 'simulations':
        dt = np.dtype([('var_id', int), ('dao_id', int), ('x', float), ('y', float),
            ('fit_template', int), ('fit_period', float), ('fit_t0', float),
            ('fit_mag1', float), ('fit_amp1', float), ('fit_sig1', float),
            ('fit_mag2', float), ('fit_amp2', float), ('fit_sig2', float),
            ('pass_visual', int), ('fit_type', 'U3'), ('fit_mode', 'U2')])
        fmt1 = '%4d %7d %8.3f %8.3f %1d %9.6f %10.4f '
        fmt2 = '%6.3f %5.3f %4.2f '
        fmt3 = '%2d %3s %2s'
        fmt = fmt1+2*fmt2+fmt3

    else:
        dt = np.dtype([('var_id', int), ('type', 'U3'), ('mode', 'U2'),
            ('template', int), ('period', float), ('t0', float),
            ('mag1', float), ('amp1', float), ('sig1', float),
            ('mag2', float), ('amp2', float), ('sig2', float),
            ('fit_template', int), ('fit_period', float), ('fit_t0', float),
            ('fit_mag1', float), ('fit_amp1', float), ('fit_sig1', float),
            ('fit_mag2', float), ('fit_amp2', float), ('fit_sig2', float),
            ('pass_index', int), ('pass_visual', int), ('fit_type', 'U3'),
            ('fit_mode', 'U2')])
        fmt1 = '%4d %3s %2s %1d %9.6f %10.4f '
        fmt2 = '%6.3f %5.3f %4.2f '
        fmt3 = '%1d %9.6f %10.4f '
        fmt4 = '%1d %1d %3s %2s'
        fmt = fmt1+2*fmt2+fmt3+2*fmt2+fmt4

    results = np.loadtxt(results_file, dtype=dt)
    results['pass_visual'][next_star] = 3
    results['fit_type'][next_star] = 'LP'
    results['fit_mode'][next_star] = 'LP'

    np.savetxt(results_file, results, fmt=fmt)
    iter += 1
    if iter == len(stars_left_to_do):
        print('Finished!')
        sys.exit()
    next_star = stars_left_to_do[iter][0]
    reset_var_states()
    plot_star_data(next_star, figure1, cat_data)

master = tk.Tk()

canvas = tk.Canvas(master, height=HEIGHT, width=WIDTH)
canvas.pack()

frame1 = tk.Frame(master, bg='#abc7d1')
frame1.place(relx=0.775, rely=0.025, relwidth=0.2, relheight=0.95)

# Variable types
rrl = tk.IntVar()
t1 = tk.Checkbutton(frame1, text='RRL', variable=rrl)
cep = tk.IntVar()
t2 = tk.Checkbutton(frame1, text='CEP', variable=cep)
ac = tk.IntVar()
t3 = tk.Checkbutton(frame1, text='AC', variable=ac)
t2c = tk.IntVar()
t4 = tk.Checkbutton(frame1, text='T2C', variable=t2c)
eb = tk.IntVar()
t5 = tk.Checkbutton(frame1, text='EB', variable=eb)
# Modes
fo = tk.IntVar()
t6 = tk.Checkbutton(frame1, text='FO', variable=fo)
fu = tk.IntVar()
t7 = tk.Checkbutton(frame1, text='FU', variable=fu)
ew = tk.IntVar()
t8 = tk.Checkbutton(frame1, text='EW', variable=ew)
ea = tk.IntVar()
t9 = tk.Checkbutton(frame1, text='EA', variable=ea)
eb2 = tk.IntVar()
t10 = tk.Checkbutton(frame1, text='EB', variable=eb2)

# Quality Flags
unsure = tk.IntVar()
t11 = tk.Checkbutton(frame1, text='Questionable', variable=unsure)


b1 = tk.Button(frame1, text='Variable', command=pass_star)
b2 = tk.Button(frame1, text='Not Variable', command=fail_star)
b4 = tk.Button(frame1, text='Previous Star', command=go_back)
b3 = tk.Button(frame1, text='Long Period', command=lpv_star)
b5 = tk.Button(frame1, text='Next Star', command=move_on)

l1 = tk.Label(frame1, text='Variable Type')
l1.place(relx=0.05, rely=0.02)
t1.place(relx=0.05, rely=0.07)
t2.place(relx=0.25, rely=0.07)
t3.place(relx=0.5, rely=0.07)
t4.place(relx=0.05, rely=0.12)
t5.place(relx=0.25, rely=0.12)

l2 = tk.Label(frame1, text='Mode/Subtype')
l2.place(relx=0.05, rely=0.2)
t6.place(relx=0.05, rely=0.25)
t7.place(relx=0.25, rely=0.25)
t8.place(relx=0.05, rely=0.3)
t9.place(relx=0.25, rely=0.3)
t10.place(relx=0.5, rely=0.3)

l3 = tk.Label(frame1, text='Quality Flags')
l3.place(relx=0.05, rely=0.38)
t11.place(relx=0.05, rely=0.43)

b1.place(relx=0.5, rely=0.55, relwidth=0.95, relheight=0.08, anchor='n')
b2.place(relx=0.5, rely=0.64, relwidth=0.95, relheight=0.08, anchor='n')
b3.place(relx=0.5, rely=0.73, relwidth=0.95, relheight=0.08, anchor='n')
b4.place(relx=0.5, rely=0.82, relwidth=0.95, relheight=0.08, anchor='n')
b5.place(relx=0.5, rely=0.91, relwidth=0.95, relheight=0.08, anchor='n')

figure1 = plt.figure(constrained_layout=True, figsize=(10,6))

plot_star_data(next_star, figure1, cat_data)


master.mainloop()
