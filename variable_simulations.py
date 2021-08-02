import numpy as np
import matplotlib.pyplot as plt
#import statsmodels.api as sm
#from astropy.stats import sigma_clip
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.filterwarnings('ignore')
#import sys
#import os
import random
import seaborn as sns
sns.set_theme()

# Function that describes the period distribution of a particular type of
# variable star in the LMC
def get_period_distribution(type, ogle_dir):

    # Read in OGLE catalog of variables stars
    dt = np.dtype([('id', 'U40'), ('I', float), ('V', float), ('P', float),
        ('eP', float), ('t0', float), ('Iamp', float), ('R21', float),
        ('phi21', float), ('R31', float), ('phi31', float)])

    if type == 'CEP':
        fo_vars = np.loadtxt(ogle_dir+'ccep/cep1O.dat', dtype=dt)
        fu_vars = np.loadtxt(ogle_dir+'ccep/cepF.dat', dtype=dt)

        bins = np.linspace(0.2, 40.2, 401)

    elif type == 'AC':
        fo_vars = np.loadtxt(ogle_dir+'acep/acep1O.dat', dtype=dt)
        fu_vars = np.loadtxt(ogle_dir+'acep/acepF.dat', dtype=dt)

        bins = np.linspace(0.3, 3.0, 28)

    elif type == 'T2C':
        t2c = np.loadtxt(ogle_dir+'t2cep/t2cep.dat', dtype=dt)
        dt2 = np.dtype([('id', 'U40'), ('type', 'U10')])
        t2c2 = np.loadtxt(ogle_dir+'t2cep/ident.dat', dtype=dt2, usecols=(0,1))
        fo_vars = t2c[t2c2['type'] == 'BLHer']
        fu_vars = t2c[t2c2['type'] == 'WVir']

        bins = np.linspace(1.0, 40.0, 79)

    elif type == 'RRL':
        fu_vars = np.loadtxt(ogle_dir+'rrl/RRab.dat', dtype=dt)
        fo_vars = np.loadtxt(ogle_dir+'rrl/RRc.dat', dtype=dt)

        bins = np.linspace(0.2, 1.0, 81)


    bins_centered = (bins[:-1] + bins[1:]) / 2
    bin_size = bins[1] - bins[0]

    fig, ax = plt.subplots(1,1)

    sns.histplot(fu_vars['P'], bins=bins, color='xkcd:rose')
    sns.histplot(fo_vars['P'], bins=bins, color='xkcd:steel blue')

    n_fu, _ = np.histogram(fu_vars['P'], bins=bins,
        density=True)
    n_fo, _ = np.histogram(fo_vars['P'], bins=bins,
        density=True)
    ax.set_xlabel('Period')
    ax.set_ylabel('Normalized N')
    plt.title('{} simulated period distribution'.format(type))
    plt.show()

    # probability of each period bin
    p_prob_fu = n_fu*bin_size
    p_prob_fo = n_fo*bin_size

    return bins_centered, p_prob_fo, p_prob_fu


# Function that describes the amplitude distribution of a particular type of
# variable star in the LMC
def get_amp_distribution(type, ogle_dir):

    # Read in OGLE catalog of variables stars
    dt = np.dtype([('id', 'U40'), ('I', float), ('V', float), ('P', float),
        ('eP', float), ('t0', float), ('Iamp', float), ('R21', float),
        ('phi21', float), ('R31', float), ('phi31', float)])

    if type == 'CEP':
        fo_vars = np.loadtxt(ogle_dir+'ccep/cep1O.dat', dtype=dt)
        fu_vars = np.loadtxt(ogle_dir+'ccep/cepF.dat', dtype=dt)

    elif type == 'AC':
        fo_vars = np.loadtxt(ogle_dir+'acep/acep1O.dat', dtype=dt)
        fu_vars = np.loadtxt(ogle_dir+'acep/acepF.dat', dtype=dt)

    elif type == 'T2C':
        t2c = np.loadtxt(ogle_dir+'t2cep/t2cep.dat', dtype=dt)
        dt2 = np.dtype([('id', 'U40'), ('type', 'U10')])
        t2c2 = np.loadtxt(ogle_dir+'t2cep/ident.dat', dtype=dt2, usecols=(0,1))
        fo_vars = t2c[t2c2['type'] == 'BLHer']
        fu_vars = t2c[t2c2['type'] == 'WVir']

    elif type == 'RRL':
        fu_vars = np.loadtxt(ogle_dir+'rrl/RRab.dat', dtype=dt)
        fo_vars = np.loadtxt(ogle_dir+'rrl/RRc.dat', dtype=dt)


    if (type == 'RRL') | (type == 'CEP'):
        bins = np.linspace(0.0, 1.0, 101)
    else:
        bins = np.linspace(0.0, 1.0, 11)
    bins_centered = (bins[:-1] + bins[1:]) / 2
    bin_size = bins[1] - bins[0]

    fig, ax = plt.subplots(1,1)
    sns.histplot(fu_vars['Iamp'], bins=bins, color='xkcd:rose')
    sns.histplot(fo_vars['Iamp'], bins=bins, color='xkcd:steel blue')

    n_fu, _ = np.histogram(fu_vars['Iamp'], bins=bins,
        density=True)
    n_fo, _ = np.histogram(fo_vars['Iamp'], bins=bins,
        density=True)
    ax.set_xlabel('Amplitude')
    ax.set_ylabel('Normalized N')
    plt.title('{} simulated amp distribution'.format(type))
    plt.show()

    # probability of each amplitude bin
    amp_prob_fu = n_fu*bin_size
    amp_prob_fo = n_fo*bin_size

    return bins_centered, amp_prob_fo, amp_prob_fu

# Generates the parameters of simulated stars
def simulate_params(var_types, num_stars, ogle_dir, id_offset=0, seed=0):

    dt = np.dtype([('id', int), ('type', 'U3'), ('mode', 'U2'),
        ('period', float), ('mag', float), ('amp', float)])
    sim_data = np.zeros(num_stars*2*len(var_types), dtype=dt)

    for i in range(len(var_types)):
        type = var_types[i]

        ### First, pick period
        bins, prob_fo, prob_fu = get_period_distribution(type, ogle_dir)
        # pick approximate periods from distribution
        if seed != 0:
            np.random.seed(seed)
        periods_fu = np.random.choice(bins, num_stars, p=prob_fu)
        periods_fo = np.random.choice(bins, num_stars, p=prob_fo)
        bin_size = bins[1] - bins[0]
        periods_fu += np.random.normal(0, bin_size, num_stars)
        periods_fo += np.random.normal(0, bin_size, num_stars)

        # put in hard limits for periods of certain types
        # NOTE: for T2C, fo/fu is actually BL Her/W Virginis
        if type == 'T2C':
            for j in range(num_stars):
                # this causes a bit of weirdness, but ok for now
                periods_fo[j] = np.max([periods_fo[j], 1.0])
                periods_fo[j] = np.min([periods_fo[j], 4.0])
                periods_fu[j] = np.max([periods_fu[j], 4.0])

        ### Now do amplitudes
        if type == 'RRL':
            # Use period-amplitude relation to assign amplitudes
            # period-amplitude relaiton is defined by OGLE IV LMC RRL
            PA_coeff_fu = np.array([-1.58, -2.99, -0.08])
            PA_coeff_fo = np.array([-4.11, -3.94, -0.66])
            PA_fu_sig = 0.14
            PA_fo_sig = 0.06

            fu_amp = PA_coeff_fu[0]*np.log10(periods_fu)**2 + \
                PA_coeff_fu[1]*np.log10(periods_fu) + PA_coeff_fu[2]
            fo_amp = PA_coeff_fo[0]*np.log10(periods_fo)**2 + \
                PA_coeff_fo[1]*np.log10(periods_fo) + PA_coeff_fo[2]

            fu_amp += np.random.normal(0, PA_fu_sig, num_stars)
            fo_amp += np.random.normal(0, PA_fo_sig, num_stars)

        # no other variables follow a period-amplitude relation, so
        # randomly assign an amplitude
        else:
            # Pick amplitudes from probability distribution
            bins, prob_fo, prob_fu = get_amp_distribution(type, ogle_dir)

            fu_amp = np.random.choice(bins, num_stars, p=prob_fu)
            fo_amp = np.random.choice(bins, num_stars, p=prob_fo)

            # perturb the amplitudes
            bin_size = bins[1] - bins[0]
            fu_amp += np.random.normal(0, bin_size, num_stars)
            fo_amp += np.random.normal(0, bin_size, num_stars)

        ### Now assign mean magnitudes

        if type == 'RRL':
            # Use I band PL relation to assign magnitudes from periods
            fu_pl_sig = 0.15
            fo_pl_sig = 0.16

            fo_mag = -2.014*np.log10(periods_fo) + 17.743
            fu_mag = -1.889*np.log10(periods_fu) + 18.164

            # add some noise
            fo_mag += np.random.normal(0, fo_pl_sig, num_stars)
            fu_mag += np.random.normal(0, fu_pl_sig, num_stars)

        if type == 'CEP':

            fu_pl_sig = 0.15
            fo_pl_sig = 0.16

            fo_mag = -3.311*(np.log10(periods_fo)-1.0) + 12.897
            fu_mag = -2.912*(np.log10(periods_fu)-1.0) + 13.741

            fo_mag += np.random.normal(0, fo_pl_sig, num_stars)
            fu_mag += np.random.normal(0, fu_pl_sig, num_stars)

        if type == 'AC':

            fu_pl_sig = 0.23
            fo_pl_sig = 0.16

            fo_mag = -3.302*np.log10(periods_fo) + 16.656
            fu_mag = -2.962*np.log10(periods_fu) + 17.368

            fo_mag += np.random.normal(0, fo_pl_sig, num_stars)
            fu_mag += np.random.normal(0, fu_pl_sig, num_stars)

        if type == 'T2C':

            fo_pl_sig = 0.40
            fu_pl_sig = 0.40

            fo_mag = -2.033*np.log10(periods_fo) + 18.015
            fu_mag = -2.033*np.log10(periods_fu) + 18.015

            fo_mag += np.random.normal(0, fo_pl_sig, num_stars)
            fu_mag += np.random.normal(0, fu_pl_sig, num_stars)


        start_fo = i*num_stars*2
        start_fu = i*num_stars + (i+1)*num_stars
        stop_fu = (i+1)*num_stars*2
        sim_data['id'] = np.arange(num_stars*2*len(var_types))+1+id_offset
        sim_data['type'][start_fo:stop_fu] = type
        sim_data['mode'][start_fo:start_fu] = 'FO'
        sim_data['mode'][start_fu:stop_fu] = 'FU'
        sim_data['period'][start_fo:start_fu] = periods_fo
        sim_data['period'][start_fu:stop_fu] = periods_fu
        sim_data['mag'][start_fo:start_fu] = fo_mag
        sim_data['mag'][start_fu:stop_fu] = fu_mag
        sim_data['amp'][start_fo:start_fu] = fo_amp
        sim_data['amp'][start_fu:stop_fu] = fu_amp

    return sim_data

# Generates the lcv files of simulated stars, which can be run through your
# fitting algorithm.
# Will create true_params.txt, as well as .fitlc files for each of your simulated
# stars.
def make_simulated_lcv(sim_data, mjds, filters, exptimes, mu, scatter_array, error_key,
    error_array, output_dir, period_search_range=[0.1, 3.0], append=False, rrl_cutoff=99,
    template_file = 'var_templates.txt', seed=0):

    if seed != 0:
        np.random.seed(seed)

    num_stars = len(sim_data['id'])

    if append == False:
        params = open(output_dir+'true_params.txt', 'w')
        params.write('#id  type mode temp period   t0         mag1   amp1 sig1 mag2   amp2 sig2\n')
    else:
        params = open(output_dir+'true_params.txt', 'a')

    # transform parameters to desired filters
    # first define shift in zero point according to the approximate distance modulus
    zp_shift = mu - 18.477 # estimated target distance modulus - LMC distance modulus

    # Now transform magnitudes into desired filters using period color relations
    # and amplitudes using amplitude ratios
    filters_unique = np.unique(filters)
    num_filters = len(filters_unique)
    mags = np.zeros((num_stars, num_filters))
    amps = np.zeros((num_stars, num_filters))
    for i in range(num_filters):

        # need to add in other colors!
        filter = filters_unique[i]

        if filter == 'B':
            rrl = sim_data['type'] == 'RRL'
            if len(sim_data['id'][rrl]) > 0:
                logP = np.log10(sim_data['period'][rrl])
                fo = sim_data['mode'][rrl] == 'FO'
                logP[fo] += 0.127
                color = 1.072 + 1.325*logP # CRRP (B-I)
                color += np.random.normal(0, 0.10, len(logP))
                mags[rrl,i] = sim_data['mag'][rrl] + zp_shift + color
                amps[rrl,i] = sim_data['amp'][rrl]*2.06 # amp ratio from CRRP RRL

            cep = ~rrl
            if len(sim_data['id'][cep]) > 0:
                P = sim_data['period'][cep]
                fo = sim_data['mode'][cep] == 'FO'
                logP = np.log10(sim_data['period'][cep])
                color = 0.636*logP + 0.685 # Sandage 2009
                color += np.random.normal(0, 0.10, len(logP))
                mags[cep,i] = sim_data['mag'][cep] + zp_shift + color
                amps[cep,i] = sim_data['amp'][cep]*2.06 # amp ratio from CRRP RRL

        if filter == 'V':
            rrl = sim_data['type'] == 'RRL'
            if len(sim_data['id'][rrl]) > 0:
                logP = np.log10(sim_data['period'][rrl])
                fo = sim_data['mode'][rrl] == 'FO'
                logP[fo] += 0.127
                color = 0.606 + 0.676*logP # CRRP (V-I)
                color += np.random.normal(0, 0.05, len(logP))
                mags[rrl,i] = sim_data['mag'][rrl] + zp_shift + color
                amps[rrl,i] = sim_data['amp'][rrl]*1.57 # amp ratio from CRRP RRL

            cep = ~rrl
            if len(sim_data['id'][cep]) > 0:
                P = sim_data['period'][cep]
                fo = sim_data['mode'][cep] == 'FO'
                logP = np.log10(sim_data['period'][cep])
                color = 0.276*logP + 0.450 # Sandage 2009 - apply to all types of cepheids
                color += np.random.normal(0, 0.08, len(logP))
                mags[cep,i] = sim_data['mag'][cep] + zp_shift + color
                amps[cep,i] = sim_data['amp'][cep]*1.57 # amp ratio from CRRP RRL

    #    if filter == 'R':

        if filter == 'I':
            mags[:,i] = sim_data['mag'] + zp_shift
            amps[:,i] = sim_data['amp']

    #    if filter == 'J':
    #    if filter == 'H':
    #    if filter == 'K':
    #    if filter == '[3.6]':
    #    if filter == '[4.5]':

    # force amplitude to be positive
    amps = np.abs(amps)

    fo_temps = np.array([6])
    fu_temps = np.array([0, 1, 2,3, 4, 5, 7, 8])
    template_number = np.zeros(num_stars)
    template_number[sim_data['mode'] == 'FO'] = np.random.choice(fo_temps, int(num_stars/2))
    template_number[sim_data['mode'] == 'FU'] = np.random.choice(fu_temps, int(num_stars/2))

    mjd_range = np.array([np.min(mjds), np.max(mjds)])
    t0 = np.random.uniform(low=mjd_range[0], high=mjd_range[1], size=num_stars)

    # Get template light curves
    dt = np.dtype([('phase', float), ('ab1', float), ('ab2', float),
        ('ab3', float), ('ab4', float), ('ab5', float), ('ab6', float),
        ('c', float), ('ab7', float), ('ab8', float)])
    templates = np.loadtxt(template_file, dtype=dt)

    templates_new = np.c_[templates['ab1'], templates['ab2'], templates['ab3'],
        templates['ab4'], templates['ab5'], templates['ab6'], templates['c'],
        templates['ab7'], templates['ab8']]

    templates_mean = np.mean(templates_new, axis=0)

    for i in range(num_stars):

        ph = templates['phase']# + np.random.uniform()
        #ncycles = 50 ### make smarter
        # calculate how many pulsation cycles we need for this star
        mjd_window = mjd_range[1] - mjd_range[0]
        ncycles = int(np.ceil(mjd_window/sim_data['period'][i]))*2

        mags_all = np.zeros(len(mjds))
        errs_all = np.zeros(len(mjds))

        xx = int(template_number[i])
        lcv_scatter = np.zeros(num_filters)
        for j in range(num_filters):

            # build template with correct period, amp, mean mag
            mag_temp = (templates_new[:,xx]-templates_mean[xx])*amps[i,j] + mags[i,j]

            t_one_cycle = ph*sim_data['period'][i] + t0[i]
            t_repeated = []
            for k in np.arange(0,ncycles):
                t_repeated = np.append(t_repeated, t_one_cycle+(k-ncycles/2)*sim_data['period'][i])
            mag_repeated = np.tile(mag_temp, ncycles)

            # sample template at proper mjds
            mjds_in_filter = mjds[filters == filters_unique[j]]
            num_obs_in_filter = len(mjds_in_filter)
            mag_sim = np.interp(mjds_in_filter, t_repeated, mag_repeated)


            # add in scatter to light curve
            # scatter is a function of magnitude, as described in inputs
            # to this functiion
            if scatter_array.ndim == 2:
                lcv_scatter[j] = np.random.normal(scatter_array[j,0], scatter_array[j,1])
                lcv_scatter[j] = np.abs(lcv_scatter[j])
            elif scatter_array.ndim == 3:
                temp_scatter = np.interp(mags[i,j], scatter_array[j,:,0], scatter_array[j,:,1])
                temp_sig = np.interp(mags[i,j], scatter_array[j,:,0], scatter_array[j,:,2])
                lcv_scatter[j] = np.random.normal(temp_scatter, temp_sig, 1)
                lcv_scatter[j] = np.abs(lcv_scatter[j])

            mag_sim += np.random.normal(0, lcv_scatter[j], num_obs_in_filter)

            mags_all[filters == filters_unique[j]] = mag_sim


            # assign photometric uncertainty to each point
            # photometric uncertainty is also a function of magnitude
            if error_array.ndim == 3:

                err_sim = np.zeros(len(mag_sim))

                for k in range(num_obs_in_filter):

                    this_exptime = exptimes[filters == filters_unique[j]][k]

                    # find right row in error_array for this filter and exptime
                    row = (error_key['filter'] == filters_unique[j]) & \
                        (error_key['exptime'] == this_exptime)

                    temp_err = np.interp(mag_sim[k], error_array[row,:,0][0], error_array[row,:,1][0])
                    temp_sig = np.interp(mag_sim[k], error_array[row,:,0][0], error_array[row,:,2][0])

                    # put in safeguard to avoid negative numbers
                    err_sim[k] = np.random.normal(temp_err, temp_sig, 1)
                    err_sim[k] = np.abs(err_sim[k])
                    # put in error floor so we don't get unrealisticaly low numbers
                    if err_sim[k] < 0.005:
                        err_sim[k] = 0.005
            elif error_array.ndim == 2:
                err_sim = np.random.normal(error_array[j,0], error_array[j,1],
                    len(mag_sim))
                err_sim = np.abs(err_sim)

            errs_all[filters == filters_unique[j]] = err_sim

        # save simulated light curve into file
        filters_int = np.zeros(len(filters))
        for j in range(num_filters):
            filters_int[filters == filters_unique[j]] = j
        data_save = np.c_[filters_int, mjds, mags_all, errs_all]

        f = open(output_dir+'sim_{}.fitlc'.format(sim_data['id'][i]), 'w')
        if rrl_cutoff < 99:
            if mags[i,-1] > rrl_cutoff:
                f.write('amoeba\nperst=0.1\npered=1.0\n')
            else:
                f.write('amoeba\nperst={}\npered={}\n'.format(period_search_range[0],
                    period_search_range[1]))
        else:
            f.write('amoeba\nperst={}\npered={}\n'.format(period_search_range[0],
                period_search_range[1]))
        np.savetxt(f, data_save, fmt='%2i %12.6f %7.4f %6.4f')
        f.close()


        format_1 = '{:4} {:4} {:4} {:4} {:8.5f} {:10.4f}'
        string_1 = format_1.format(int(sim_data['id'][i]), sim_data['type'][i],
            sim_data['mode'][i], int(template_number[i]), sim_data['period'][i],
            t0[i])
        string_2 = ''
        for j in range(num_filters):
            format_2 = ' {:6.3f} {:4.2f} {:4.2f}'
            string_2 += format_2.format(mags[i,j], amps[i,j], lcv_scatter[j])

        line = string_1 + string_2 + '\n'
        params.write(line)
    params.close()


def make_simulated_lcv_single_snr(sim_data, mjds, filters, exptimes, mu,
    scatter_mean_snr, scatter_std_snr, error_mean_snr, error_std_snr, output_dir,
    period_search_range=[0.1, 1.0], append=False, rrl_cutoff=99,
    template_file='var_templates.txt'):

    num_stars = len(sim_data['id'])

    if append == False:
        params = open(output_dir+'true_params.txt', 'w')
        params.write('#id  type mode temp period   t0         mag1   amp1 sig1 mag2   amp2 sig2\n')
    else:
        params = open(output_dir+'true_params.txt', 'a')

    # transform parameters to desired filters
    # first transform from I band LMC magnitude to I band target magnitude
    zp_shift = mu - 18.477 # estimated target distance modulus - LMC distance modulus

    # Now transform magnitudes into desired filters using period color relations
    # and amplitudes using amplitude ratios
    filters_unique = np.unique(filters)
    num_filters = len(filters_unique)
    mags = np.zeros((num_stars, num_filters))
    amps = np.zeros((num_stars, num_filters))
    for i in range(num_filters):

        # need to add in other colors!
        filter = filters_unique[i]

        if filter == 'B':
            rrl = sim_data['type'] == 'RRL'
            if len(sim_data['id'][rrl]) > 0:
                logP = np.log10(sim_data['period'][rrl])
                fo = sim_data['mode'][rrl] == 'FO'
                logP[fo] += 0.127
                color = 1.072 + 1.325*logP # CRRP (B-I)
                color += np.random.normal(0, 0.10, len(logP))
                mags[rrl,i] = sim_data['mag'][rrl] + zp_shift + color
                amps[rrl,i] = sim_data['amp'][rrl]*2.06

            cep = ~rrl
            if len(sim_data['id'][cep]) > 0:
                P = sim_data['period'][cep]
                fo = sim_data['mode'][cep] == 'FO'
                logP = np.log10(sim_data['period'][cep])
                color = 0.636*logP + 0.685 # Sandage 2009
                color += np.random.normal(0, 0.10, len(logP))
                mags[cep,i] = sim_data['mag'][cep] + zp_shift + color
                amps[cep,i] = sim_data['amp'][cep]*2.06

        if filter == 'V':
            rrl = sim_data['type'] == 'RRL'
            if len(sim_data['id'][rrl]) > 0:
                logP = np.log10(sim_data['period'][rrl])
                fo = sim_data['mode'][rrl] == 'FO'
                logP[fo] += 0.127
                color = 0.606 + 0.676*logP # CRRP (V-I)   ### something is off here
                color += np.random.normal(0, 0.05, len(logP))
                mags[rrl,i] = sim_data['mag'][rrl] + zp_shift + color
                amps[rrl,i] = sim_data['amp'][rrl]*1.57

            cep = ~rrl
            if len(sim_data['id'][cep]) > 0:
                P = sim_data['period'][cep]
                fo = sim_data['mode'][cep] == 'FO'
                logP = np.log10(sim_data['period'][cep])
                color = 0.276*logP + 0.450 # Sandage 2009 - apply to all types of cepheids
                color += np.random.normal(0, 0.08, len(logP))
                mags[cep,i] = sim_data['mag'][cep] + zp_shift + color
                amps[cep,i] = sim_data['amp'][cep]*1.57

    #    if filter == 'R':

        if filter == 'I':
            mags[:,i] = sim_data['mag'] + zp_shift
            amps[:,i] = sim_data['amp']

    #    if filter == 'J':
    #    if filter == 'H':
    #    if filter == 'K':
    #    if filter == '[3.6]':
    #    if filter == '[4.5]':

    # force amplitude to be positive
    amps = np.abs(amps)

    fo_temps = np.array([6])
    fu_temps = np.array([0, 1, 2,3, 4, 5, 7, 8])
    template_number = np.zeros(num_stars)
    template_number[sim_data['mode'] == 'FO'] = np.random.choice(fo_temps, int(num_stars/2))
    template_number[sim_data['mode'] == 'FU'] = np.random.choice(fu_temps, int(num_stars/2))

    mjd_range = np.array([np.min(mjds), np.max(mjds)])
    t0 = np.random.uniform(low=mjd_range[0], high=mjd_range[1], size=num_stars)

    # Get template light curves
    dt = np.dtype([('phase', float), ('ab1', float), ('ab2', float),
        ('ab3', float), ('ab4', float), ('ab5', float), ('ab6', float),
        ('c', float), ('ab7', float), ('ab8', float)])
    templates = np.loadtxt(template_file, dtype=dt)

    templates_new = np.c_[templates['ab1'], templates['ab2'], templates['ab3'],
        templates['ab4'], templates['ab5'], templates['ab6'], templates['c'],
        templates['ab7'], templates['ab8']]

    templates_mean = np.mean(templates_new, axis=0)

    for i in range(num_stars):

        ph = templates['phase']# + np.random.uniform()
        #ncycles = 50 ### make smarter
        # calculate how many pulsation cycles we need for this star
        mjd_window = mjd_range[1] - mjd_range[0]
        ncycles = int(np.ceil(mjd_window/sim_data['period'][i]))*2

        mags_all = np.zeros(len(mjds))
        errs_all = np.zeros(len(mjds))


        xx = int(template_number[i])
        lcv_scatter = np.zeros(num_filters)
        for j in range(num_filters):

            # build template with correct period, amp, mean mag
            mag_temp = (templates_new[:,xx]-templates_mean[xx])*amps[i,j] + mags[i,j]

            t_one_cycle = ph*sim_data['period'][i] + t0[i]
            t_repeated = []
            for k in np.arange(0,ncycles):
                t_repeated = np.append(t_repeated, t_one_cycle+(k-ncycles/2)*sim_data['period'][i])
            mag_repeated = np.tile(mag_temp, ncycles)

            # sample template at proper mjds
            mjds_in_filter = mjds[filters == filters_unique[j]]
            num_obs_in_filter = len(mjds_in_filter)
            mag_sim = np.interp(mjds_in_filter, t_repeated, mag_repeated)


            # add in scatter to light curve
            temp_scatter = -99
            while temp_scatter < 0.01:
                temp_scatter = 1.09/np.random.normal(scatter_mean_snr, scatter_std_snr)
                if temp_scatter > 0.4:
                    temp_scatter = -99
            lcv_scatter[j] = temp_scatter

            # check that it is positive
            mag_sim += np.random.normal(0, lcv_scatter[j], num_obs_in_filter)

            mags_all[filters == filters_unique[j]] = mag_sim


            # assign photometric uncertainty to each point
            err_sim = np.zeros(len(mag_sim))
            for k in range(num_obs_in_filter):
                temp_err = -99
                while temp_err < 0.005:
                    temp_err = 1.09/np.random.normal(error_mean_snr, error_std_snr)
                    if temp_err > 0.4:
                        temp_err = -99
                err_sim[k] = temp_err

            errs_all[filters == filters_unique[j]] = err_sim

        # save simulated light curve into file
        filters_int = np.zeros(len(filters))
        for j in range(num_filters):
            filters_int[filters == filters_unique[j]] = j
        data_save = np.c_[filters_int, mjds, mags_all, errs_all]

        f = open(output_dir+'sim_{}.fitlc'.format(sim_data['id'][i]), 'w')
        if rrl_cutoff < 99:
            if mags[i,-1] > rrl_cutoff:
                f.write('amoeba\nperst=0.1\npered=1.0\n')
            else:
                f.write('amoeba\nperst={}\npered={}\n'.format(period_search_range[0],
                    period_search_range[1]))
        else:
            f.write('amoeba\nperst={}\npered={}\n'.format(period_search_range[0],
                period_search_range[1]))
        np.savetxt(f, data_save, fmt='%2i %12.6f %7.4f %6.4f')
        f.close()


        format_1 = '{:4} {:4} {:4} {:4} {:8.5f} {:10.4f}'
        string_1 = format_1.format(int(sim_data['id'][i]), sim_data['type'][i],
            sim_data['mode'][i], int(template_number[i]), sim_data['period'][i],
            t0[i])
        string_2 = ''
        for j in range(num_filters):
            format_2 = ' {:6.3f} {:4.2f} {:4.2f}'
            string_2 += format_2.format(mags[i,j], amps[i,j], lcv_scatter[j])

        line = string_1 + string_2 + '\n'
        params.write(line)
    params.close()
    print(np.mean(mags[:,0]), np.mean(mags[:,1]))


def compile_sims():
    f = open('sim_results.txt', 'w')
    #f = open('test.txt', 'w')

    dt = np.dtype([('id', int), ('type', 'U3'), ('mode', 'U2'),
        ('template', int), ('period', float), ('t0', float),
        ('mag1', float), ('amp1', float), ('sig1', float),
        ('mag2', float), ('amp2', float), ('sig2', float)])

    true_params = np.loadtxt('lcvs/true_params.txt', dtype=dt, skiprows=1)

    nsims = len(true_params['id'])

    dt2 = np.dtype([('id', int), ('type', 'U3'), ('mode', 'U2'),
        ('template', int), ('period', float), ('t0', float),
        ('mag1', float), ('amp1', float), ('sig1', float),
        ('mag2', float), ('amp2', float), ('sig2', float),
        ('fit_template', int), ('fit_period', float), ('fit_t0', float),
        ('fit_mag1', float), ('fit_amp1', float), ('fit_sig1', float),
        ('fit_mag2', float), ('fit_amp2', float), ('fit_sig2', float),
        ('pass1', int), ('pass2', int), ('fit_type', 'U3'), ('fit_mode', 'U2')])
    fit_params = np.zeros(nsims, dtype=dt2)

    for i in range(nsims):

        star = true_params['id'][i]
        dt = np.dtype([('template', int), ('period', float), ('epoch', float),
            ('amp1', float), ('mag1', float), ('amp2', float), ('mag2', float)])
        star_props = np.loadtxt('lcvs/sim_{}.fitlc_props'.format(star),
            dtype=dt, skiprows=1, usecols=(0,1,3,4,6,7,9))

        fit_params['id'][i] = true_params['id'][i]
        fit_params['type'][i] = true_params['type'][i]
        fit_params['mode'][i] = true_params['mode'][i]
        fit_params['template'][i] = true_params['template'][i]
        fit_params['period'][i] = true_params['period'][i]
        fit_params['t0'][i] = true_params['t0'][i]
        fit_params['mag1'][i] = true_params['mag1'][i]
        fit_params['amp1'][i] = true_params['amp1'][i]
        fit_params['sig1'][i] = true_params['sig1'][i]
        fit_params['mag2'][i] = true_params['mag2'][i]
        fit_params['amp2'][i] = true_params['amp2'][i]
        fit_params['sig2'][i] = true_params['sig2'][i]
        fit_params['fit_template'][i] = star_props['template']-1
        fit_params['fit_period'][i] = star_props['period']
        fit_params['fit_t0'][i] = star_props['epoch']
        #fit_params['fit_mag1'][i] = star_props['mag1'] # Changed to mean of template!
        fit_params['fit_amp1'][i] = star_props['amp1']
        #fit_params['fit_mag2'][i] = star_props['mag2'] # Changed to mean of template!
        fit_params['fit_amp2'][i] = star_props['amp2']
        fit_params['pass1'][i] = -1
        fit_params['pass2'][i] = -1
        fit_params['fit_type'][i] = 'XXX'
        fit_params['fit_mode'][i] = 'XX'

        star_data = np.loadtxt('lcvs/sim_{}.fitlc_phase'.format(star),
            skiprows=1, usecols=(0,1,2))

        star_fit = np.loadtxt('lcvs/sim_{}.fitlc_fit'.format(star),
            skiprows=1, usecols=(0,1,2))

        f1 = star_data[:,0] == 0
        calc = np.interp(star_data[f1,1], star_fit[:,0], star_fit[:,1])
        resid = star_data[f1,2] - calc
        fit_params['fit_sig1'][i] = np.std(resid)

        f2 = star_data[:,0] == 1
        calc = np.interp(star_data[f2,1], star_fit[:,0], star_fit[:,2])
        resid = star_data[f2,2] - calc
        fit_params['fit_sig2'][i] = np.std(resid)

        # this isn't intensity mean for simulations!
        fit_params['fit_mag1'][i] = np.mean(star_fit[:,1])
        fit_params['fit_mag2'][i] = np.mean(star_fit[:,2])

    fmt1 = '%4d %3s %2s %1d %9.6f %10.4f '
    fmt2 = '%6.3f %5.3f %4.2f '
    fmt3 = '%1d %9.6f %10.4f '
    fmt4 = '%1d %1d %3s %2s'
    fmt = fmt1+2*fmt2+fmt3+2*fmt2+fmt4
    np.savetxt(f, fit_params, fmt=fmt)

    f.close()

    # shuffle the results, to prevent bias in visual check.
    lines = open('sim_results.txt').readlines()
    random.shuffle(lines)
    open('sim_results.txt', 'w').writelines(lines)
