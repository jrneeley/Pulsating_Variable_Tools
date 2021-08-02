'''
    File name: variable_identification.py
    Author: Jillian Neeley
    Date last modified: 8/1/2021
    Python Version: 3.6
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from progress.bar import Bar
from astropy.stats import sigma_clip
#import AstroTools as at
import os

# Calculate the weighted mean, as defined by Peter Stetson
def stetson_robust_mean(mags, errs):

    if len(mags) == 1:
        return mags

    # As a first guess, calculate the weighted mean the usual way
    weights = 1/errs**2
    initial_guess = np.sum(mags*weights)/np.sum(weights)
    n = len(mags)

    # Iteratively reweight points based on their difference from the weighted
    # mean, and recalculate the weighed mean. Stop when this converges.
    diff = 99
    old_mean = initial_guess
    for i in range(1000):

        delta = np.sqrt(n/(n-1))*(mags-old_mean)/errs
        weight_factor = 1/(1+(np.abs(delta)/2)**2)
        weights = weights*weight_factor

        new_mean = np.sum(mags*weights)/np.sum(weights)

        diff = np.abs(old_mean - new_mean)
        # break the loop if the weighted mean has converged
        if diff < 0.00001:
            break
        old_mean = new_mean

    return new_mean

# Simple function for the standard weighted mean.
def weighted_mean(mags, errs):

    # Select out indices where magnitude and error are both finite (not nan)
    finite = (~np.isnan(mags)) & (~np.isnan(errs))
    weights = 1./errs[finite]**2
    sum_weights = np.sum(weights)

    mean = np.sum(mags[finite]*weights)/sum_weights

    return mean

# A function that computes a variety of variability indices for a single star.
def compute_variability_index(filters, mjds, mags, errs, max_time=0.02):
# INPUTS:
# * filters - array of filter of each observation for this star
# * mjds - Modified Julian Dates for each observation of this star
# * mags - array of relative magnitudes for this star (output from calc_relative_phot)
# * errs - array of photometric uncertainties for this star
# * max_time - parameter to control the maximum separation in MJD to be counted
# as a pair of observations
# OUTPUTS:
# * var_indices - array of all computed variability indices: Stetson J, chi
# squared, weighted standard deviation (change?), median absolute deviation, and
# robust median statistic

    # first filter out nan mags
    filters = filters[~np.isnan(mags)]
    mjds = mjds[~np.isnan(mags)]
    errs = errs[~np.isnan(mags)]
    mags = mags[~np.isnan(mags)]

    if len(mags) <= 1:
        return np.repeat(np.nan, 5)


    # Determine unique filters in dataset
    filter_list = np.unique(filters)
    n_filts = len(filter_list)

    # How many total observations (in all filters) do we have?
    num_obs_total = len(mags)

    # set up output array
    var_indices = np.zeros(5)

    ####### Compute the Stetson J index #####

    # set up arrays for the weighted mean, the number of observations in each
    # filter, and the filter number of each observation (integer for indexing)
    weighted_means = np.zeros(n_filts)
    num_obs = np.zeros(n_filts)
    filter_num = np.zeros(num_obs_total)

    for i in range(n_filts):

        # select out only the observations in this filter.
        # NOTE: This creates a boolean
        # array of the same length as filters, where the elements are true
        # if it satisfies the condition. It can be used to index another array -
        # I do this a lot in python, but don't remember if you can do the same
        # in IDL without using the where() function
        f = filters == filter_list[i]
        # compute the weighted mean in each filter
        weighted_means[i] = stetson_robust_mean(mags[f], errs[f])
        num_obs[i] = float(len(mags[f]))
        filter_num[f] = i


    # This index requires the data to be sorted by time, so do that
    order = np.argsort(mjds)
    mags_temp = mags[order]
    errs_temp = errs[order]
    mjds_temp = mjds[order]
    filt_temp = filter_num[order]

    P = 0
    n_pairs = 0
    skip_next = False
    for i in range(num_obs_total-1):

        # If skip_next == True, then this observation has already been counted
        # in a pair, so change it back to False and move on to the next
        # iteration of the loop
        if skip_next == True:
            skip_next = False
            continue

        # Check if the current observation and the next one were taken close
        # together in time. If they are within your maximum time difference,
        # count them as a pair
        if mjds_temp[i+1] - mjds_temp[i] <= max_time:

            # Check which filters the observations in our pair were taken in, so
            # we compare them to the appropriate weighted mean.
            # This allows for the possibility that these two observations are
            # from the same or different filters
            fnum1 = int(filt_temp[i])
            fnum2 = int(filt_temp[i+1])

            temp1 = (mags_temp[i] - weighted_means[fnum1])/errs_temp[i]
            delta1 = np.sqrt(num_obs[fnum1]/(num_obs[fnum1]-1))*temp1

            temp2 = (mags_temp[i+1] - weighted_means[fnum2])/errs_temp[i+1]
            delta2 = np.sqrt(num_obs[fnum2]/(num_obs[fnum2]-1))*temp2
            # Stetson math
            P += np.sign(delta1*delta2)*np.sqrt(np.abs(delta1*delta2))
            # We paired observation i and i+1, so we will need to skip the
            # next iteration
            skip_next = True


        # This observation is not part of a pair, (could be an isolated
        # observation or part of a grouping of an odd nubmer of observations)
        # and is now treated as a single observation.
        else:
            fnum = int(filt_temp[i])
            temp = (mags_temp[i] - weighted_means[fnum])/errs_temp[i]
            delta = np.sqrt(num_obs[fnum]/(num_obs[fnum]-1))*temp

            P += np.sign(delta*delta-1)*np.sqrt(np.abs(delta*delta-1))
            skip_next = False

        # We don't actually need this variable anymore, but it was useful
        # for testing
        n_pairs += 1

    # Set first variability index to the Stetson J value
    var_indices[0] = P

    ####### Compute chi squared index ###########

    sum = 0

    # Loop through all filters, but a single chi squared using observations
    # in all filters is computed in the end
    for i in range(n_filts):
        f = filters == filter_list[i]
        # Let's use traditional weighted mean this time.
        weighted_mean_mag = weighted_mean(mags[f], errs[f])
        # Use += so we can combine information from different filters
        sum += np.sum((mags[f]-weighted_mean_mag)**2/errs[f]**2)

    chi_squared = 1./(float(num_obs_total)-1)*sum
    # Set second variability index to the
    var_indices[1] = chi_squared

    ######## Compute weighted standard deviation? #####

    # Weighted standard deviation is already computed in the CCHP catalog, but
    # independently for differnt filters. Need to decide if I want to give a
    # combined estimate here, or to replace this with a different index.
    # For now, lets just set it to nan, so we know we haven't done anything yet.

    var_indices[2] = np.nan


    ######## Calculate median absolute deviation (MAD) #########

    # set up empty array for median magnitude. This array has same length as
    # mags, but each element will be the median of all magnitudes of the
    # corresponding filter.
    median_mags = np.zeros(num_obs_total)

    for i in range(n_filts):
        f = filters == filter_list[i]
        # get the median magnitude in this filter, and copy it into an array,
        # whose corresponding elements in mags are the same filter.
        median_mags[f] = np.nanmedian(mags[f])

    absolute_deviation = np.abs(mags - median_mags)
    mad = np.nanmedian(absolute_deviation)

    # Set 4th variability index to median absolute devation
    var_indices[3] = mad

    ######## Calculate Robust median statistic  (RoMS) ########

    sum = 0
    for i in range(n_filts):
        # get all observations in this filter
        f = filters == filter_list[i]
        # Use += so we can combine observations from different filters.
        sum += np.sum(np.abs(mags[f] - np.median(mags[f]))/errs[f])

    # normalize by the total number of observations
    RoMS = sum/(float(num_obs_total)-1)

    # set 5th variability index to robust median statistic
    var_indices[4] = RoMS

    return var_indices


# Helper function to do a binned sigma clip
def binned_sigma_clip(xdata, ydata, bins=10, sigma=3, iters=5):

    # check for and remove nan values

    clipped = np.zeros(len(ydata), dtype=int)
    ind = np.arange(len(ydata))

    # check for and remove nan values
    good = (~np.isnan(ydata))

    std, edges, num = stats.binned_statistic(xdata[good],
        ydata[good], 'std', bins=bins)

    for i in range(bins):

        in_bin = (xdata[good] >= edges[i]) & (xdata[good] < edges[i+1])
        filtered_data = sigma_clip(ydata[good][in_bin], sigma=sigma, maxiters=iters)
        s = ind[good][in_bin]
        clipped[s] = filtered_data.mask*1

    return clipped, edges, std

# Function to identify variable star candidates, and make some diagnostic plots.
# Calls the function compute_variability_index (for each star), and then
# identifies potential variables from the variability index of all stars in the
# catalog. Each index has their own set of criteria to flag potential variables.
def find_vars(ids, chip, data, mags_relative, frame_flags, clean=False, plot=True):
# INPUTS:
# * ids, chip, and data are outputs from read_catalog
# * mags_relative, frame_flags are outputs from calc_relative_phot
# * clean is a flag used to determine if you want to remove flagged frames or not
# * plot is a flag to determine if diagnostic plots are made or not
# OUTPUTS:
# variable_flag array that identifies potential variable candidates based on
# each individual variability index (0=not variable 1=variable)

    # Make directory for variable search plots, if it doesn't already exist
    # This is for my setup - will change based on pipeline architechture
    if not os.path.exists('var_search'):
        os.makedirs('var_search')

    # How many stars do we have?
    nstars = len(data['id'])
    # What filters are we working with?
    filters = data['filter'][0]
    n_filters = len(filters)

    # check if we have already computed the variability index before - time
    # saver for me, but won't be necessary in pipeline
    if os.path.exists('var_search/vi.txt'):
        print('Variability index already calculated.')
        vi = np.loadtxt('var_search/vi.txt')

    else:
        # Set up array for variability index of all stars
        vi = np.zeros((nstars,5))
        # progress bar for my sanity
        bar = Bar('Calculating variability', max=nstars)
        # Loop through all stars in catalog - this makes it pretty slow, but
        # I couldn't think of an obvious alternative
        for i in range(nstars):
            # use relative photometry
            mag_all = mags_relative[i]

            # First check if the star was detected in enough frames, if
            # not, don't bother computing the variability index.
            fail = 0
            for j in range(n_filters):
                f = data['filters'][i] == filters[j]
                mags_in_filter = mag_all[f]
                # how many magnitudes should there be for this star/filter?
                num_possible_obs = len(mags_in_filter)
                # get number of finite magnitudes
                num_actual_obs = len(mags_in_filter[~np.isnan(mags_in_filter)])
                # If this star was not detected in at least half the number of
                # frames (in any filter), move on and don't compute variability index
                if num_actual_obs < num_possible_obs/2:
                    fail = 1
            if fail == 1:
                vi[i,0] = np.nan
                bar.next()
            else:
                # If we are not removing flagged frames, only remove nan values
                if clean == False:
                    keep = ~np.isnan(mag_all)
                # If we are removing flagged frames, do that here
                else:
                    keep = (~np.isnan(mag_all)) & (frame_flags[i] == 0)

                # calculate various variability indices for this star
                vi[i] = compute_variability_index(data['filters'][i][keep],
                    data['mjd_obs'][i][keep], mag_all[keep],
                    data['emag'][i][keep], max_time=0.02)

                bar.next()
        bar.finish()

        # write to file - convenience for me, in CCHP pipeline, results
        # will be stored in fits catalog
        np.savetxt('var_search/vi.txt', vi, fmt='%10f %10f %10f %10f %10f')


    # Pick candidates based on variability index. Originally I did this
    # independently on each chip, but after more testing I think it is OK to
    # do all at once

    # set up variable flag arrays
    variable_flag = np.zeros((nstars, 5))

    vi_copy = np.copy(vi[:,0])

    ################ STETSON J CRITERIA (index 0) ####################
    # First, accept variable candidates with Stetson index > 15
    group1 = vi_copy > 15

    # mask the accepted candidates for the next step
    vi_copy[group1] = np.nan

    # Use sigma clipping to determine variability cutoff
    clipped, bins, std = binned_sigma_clip(data['mag_sw'][:,1], vi_copy,
        bins=10, sigma=5, iters=5)
    bins_centered = (bins[:-1] + bins[1:])/2

    # identify positive values
    pos = vi_copy > 0
    # add to candidates list stars that are clipped and positive
    npass = clipped  + pos*1
    group2 = npass == 2

    # combine the two criteria for list of all candidates
    candidates = group1 | group2

    variable_flag[candidates,0] = 1

    # plotting stuff
    if plot == True:

        fig, ax = plt.subplots(1,1)
        ax.scatter(data['mag_sw'][:,1], vi[:,0],
            marker='.', color='black', alpha=0.3, s=1, rasterized=True)
        ax.errorbar(bins_centered, np.zeros(len(bins_centered)), yerr=std,
            fmt='.', color='blue')
        ax.scatter(data['mag_sw'][:,1][group1], vi[:,0][group1],
            marker='.', color='r', alpha=0.5, s=5)
        ax.scatter(data['mag_sw'][:,1][group2], vi[:,0][group2],
            marker='.', color='green', alpha=0.5, s=5)
        ax.set_xlabel('{} magnitude'.format(data['filter'][0,1]))
        ax.set_ylabel('Stetson J index')
        plt.savefig('var_search/var_index_stetson_j.pdf', format='pdf',
            rasterized=True, dpi=400)

        fig, ax = plt.subplots(1,1)
        color = data['mag_sw'][:,0] - data['mag_sw'][:,1]
        mag = data['mag_sw'][:,1]
        #at.AstroPlots.plot_cmd(color, mag, ylim=[np.nanmin(mag), np.nanmax(mag)],
        #    plt_axes=ax, cmap=plt.cm.Greys, cbar_scale='log', cbar_min=-1)
        ax.scatter(color, mag, s=1, alpha=0.3, color='k')
        ax.scatter(color[candidates], mag[candidates], s=4, color='r', alpha=0.7)
        ax.set_xlabel('{} - {}'.format(filters[0], filters[1]))
        ax.set_ylabel('{}'.format(filters[1]))
        plt.savefig('var_search/cmd_candidates_StetsonJ.pdf', rasterized=True,
            format='pdf')


    # Print a few useful things out to the terminal - mostly for testing
    n_candidates = len(data['id'][candidates])
    c1_candidates = candidates & (chip == 1)
    c2_candidates = candidates & (chip == 2)
    n_candidates_chip1 = len(data['id'][c1_candidates])
    n_candidates_chip2 = len(data['id'][c2_candidates])
    print('{} candidates: {} on chip 1, {} on chip 2'.format(n_candidates,
        n_candidates_chip1, n_candidates_chip2))

    ##### TO DO: Flag based on chi squared

    if plot == True:
        fig, ax = plt.subplots(1,1)
        ax.scatter(data['mag_sw'][:,1], vi[:,1],
            marker='.', color='black', alpha=0.3, s=1, rasterized=True)
        ax.set_xlabel('{} magnitude'.format(data['filter'][0,1]))
        ax.set_ylabel('$\chi^2$')
        plt.savefig('var_search/var_index_chisq.pdf', format='pdf',
            rasterized=True, dpi=400)

    ##### TO DO: Flag based on weighted standard deviation

#    if plot == True:
#        fig, ax = plt.subplots(1,1)
#        ax.scatter(data['mag_sw'][:,1], vi[:,2],
#            marker='.', color='black', alpha=0.3, s=1, rasterized=True)
#        ax.set_xlabel('{} magnitude'.format(data['filter'][0,1]))
#        ax.set_ylabel('$\sigma$')
#        plt.savefig('var_search/var_index_stddev.pdf', format='pdf',
#            rasterized=True, dpi=400)

    ##### TO DO: Flag based on MAD

    if plot == True:
        fig, ax = plt.subplots(1,1)
        ax.scatter(data['mag_sw'][:,1], vi[:,3],
            marker='.', color='black', alpha=0.3, s=1, rasterized=True)
        ax.set_xlabel('{} magnitude'.format(data['filter'][0,1]))
        ax.set_ylabel('MAD')
        plt.savefig('var_search/var_index_mad.pdf', format='pdf',
            rasterized=True, dpi=400)

    ##### TO DO: Flag based on RoMS

    if plot == True:
        fig, ax = plt.subplots(1,1)
        ax.scatter(data['mag_sw'][:,1], vi[:,4],
            marker='.', color='black', alpha=0.3, s=1, rasterized=True)
        ax.set_xlabel('{} magnitude'.format(data['filter'][0,1]))
        ax.set_ylabel('RoMS')
        plt.savefig('var_search/var_index_roms.pdf', format='pdf',
            rasterized=True, dpi=400)

    return variable_flag
