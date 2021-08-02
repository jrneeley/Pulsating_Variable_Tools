import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic

# currently only works with two filters, with one exptime in each.
def get_properties(target, n_epochs, filters=['F475W', 'F814W'], exptimes=[1000, 1000]):

    n_filters = len(np.unique(filters))
    #n_exptimes = len(exptimes)

    multiepoch_data = np.loadtxt('{}_allcal_12_1line.raw.tran'.format(target))


    b_mags = multiepoch_data[:,3:3+n_epochs*2:2]
    b_emags = multiepoch_data[:,4:4+n_epochs*2:2]

    i_mags = multiepoch_data[:,3+n_epochs*2:-2:2]
    i_emags = multiepoch_data[:,4+n_epochs*2:-2:2]

    b_mags[b_mags > 50] = np.nan
    b_emags[b_emags > 5] = np.nan
    i_mags[i_mags > 50] = np.nan
    i_emags[i_emags > 5] = np.nan

    b_mags[b_mags == 0.0] = np.nan
    b_emags[b_emags == 0.0] = np.nan
    i_mags[i_mags == 0.0] = np.nan
    i_emags[i_emags == 0.0] = np.nan

    sharp = multiepoch_data[:,-1]

    sharp_repeat = np.zeros((len(sharp), len(b_mags[0])))
    for i in range(len(b_mags[0])):
        sharp_repeat[:,i] = sharp
    sharp_flat = sharp_repeat.flatten()

    low_sharp = np.abs(sharp) < 0.1
    low_sharp_flat = np.abs(sharp_flat) < 0.1

    dt = np.dtype([('filter', 'U5'), ('exptime', float)])
    nrows = len(filters)
    nbins = 15

    error_key = np.zeros(nrows, dtype=dt)
    for i in range(len(filters)):
        error_key['filter'][i] = filters[i]
        error_key['exptime'][i] = exptimes[i]

    error_array = np.zeros((nrows,nbins,3))
    scatter_array = np.zeros((n_filters,nbins,3))

    b_mags_flat = b_mags.flatten()
    b_emags_flat = b_emags.flatten()
    i_mags_flat = i_mags.flatten()
    i_emags_flat = i_emags.flatten()

    notnan = ~np.isnan(b_mags_flat) & ~np.isnan(b_emags_flat)
    select = notnan & low_sharp_flat

    def mad(x):
        med_x = np.nanmedian(x)
        return np.nanmedian(np.abs(x - med_x))

    phot_unc, edges, _ = binned_statistic(b_mags_flat[select], b_emags_flat[select],
                    'median', bins=nbins)
    sig_phot_unc, edges2, _ = binned_statistic(b_mags_flat[select], b_emags_flat[select],
                    statistic=mad, bins=nbins)

    centers1 = (edges[:-1] + edges[1:])/2.0

    error_array[0,:,0] = centers1
    error_array[0,:,1] = phot_unc
    error_array[0,:,2] = sig_phot_unc


    fig, ax = plt.subplots(1,1)

    ax.scatter(b_mags_flat[select], b_emags_flat[select], s=0.1, color='k', alpha=0.1)
    ax.errorbar(centers1, phot_unc, yerr=sig_phot_unc, fmt='o', color='orange')
    #ax.set_ylim(0,0.4)
    ax.set_xlabel('{} mag'.format(filters[0]))
    ax.set_ylabel('$\sigma_i$ {}'.format(filters[0]))
    ax.text(0.1, 0.9, 'Exp = {} s'.format(exptimes[0]), transform=ax.transAxes)

    plt.show()

    # repeat for i band
    notnan = ~np.isnan(i_mags_flat) & ~np.isnan(i_emags_flat)
    select = notnan & low_sharp_flat

    phot_unc, edges, _ = binned_statistic(i_mags_flat[select], i_emags_flat[select],
                    'median', bins=nbins)
    sig_phot_unc, edges2, _ = binned_statistic(i_mags_flat[select], i_emags_flat[select],
                    statistic=mad, bins=nbins)

    centers1 = (edges[:-1] + edges[1:])/2.0

    error_array[1,:,0] = centers1
    error_array[1,:,1] = phot_unc
    error_array[1,:,2] = sig_phot_unc


    fig, ax = plt.subplots(1,1)

    ax.scatter(i_mags_flat[select], i_emags_flat[select], s=0.1, color='k', alpha=0.1)
    ax.errorbar(centers1, phot_unc, yerr=sig_phot_unc, fmt='o', color='orange')
    #ax.set_ylim(0,0.4)
    ax.set_xlabel('{} mag'.format(filters[1]))
    ax.set_ylabel('$\sigma_i$ {}'.format(filters[1]))
    ax.text(0.1, 0.9, 'Exp = {} s'.format(exptimes[1]), transform=ax.transAxes)

    plt.show()

    b_mean_mags = np.nanmean(b_mags, axis=1)
    b_stds = np.nanstd(b_mags, axis=1)
    i_mean_mags = np.nanmean(i_mags, axis=1)
    i_stds = np.nanstd(i_mags, axis=1)

    notnan = ~np.isnan(b_mean_mags) & ~np.isnan(b_stds)
    select = notnan & low_sharp
    scatter, edges, _ = binned_statistic(b_mean_mags[select], b_stds[select],
        'median', bins=nbins)
    sig_scatter, edges2, _ = binned_statistic(b_mean_mags[select], b_stds[select],
        statistic=mad, bins=nbins)

    centers2 = (edges[:-1] + edges[1:])/2.0

    scatter_array[0,:,0] = centers2
    scatter_array[0,:,1] = scatter
    scatter_array[0,:,2] = sig_scatter

    fig, ax = plt.subplots(1,1)

    ax.scatter(b_mean_mags[select], b_stds[select], s=0.1, color='k', alpha=0.1)
    ax.errorbar(centers2, scatter, yerr=sig_scatter, fmt='o', color='orange')
    #ax.set_ylim(0,0.4)
    ax.set_xlabel('{} mag'.format(filters[0]))
    ax.set_ylabel('$\sigma$ {}'.format(filters[0]))

    plt.show()

    notnan = ~np.isnan(i_mean_mags) & ~np.isnan(i_stds)
    select = notnan & low_sharp
    scatter, edges, _ = binned_statistic(i_mean_mags[select], i_stds[select],
        'median', bins=nbins)
    sig_scatter, edges2, _ = binned_statistic(i_mean_mags[select], i_stds[select],
        statistic=mad, bins=nbins)

    centers2 = (edges[:-1] + edges[1:])/2.0

    scatter_array[1,:,0] = centers2
    scatter_array[1,:,1] = scatter
    scatter_array[1,:,2] = sig_scatter

    fig, ax = plt.subplots(1,1)

    ax.scatter(i_mean_mags[select], i_stds[select], s=0.1, color='k', alpha=0.1)
    ax.errorbar(centers2, scatter, yerr=sig_scatter, fmt='o', color='orange')
    #ax.set_ylim(0,0.4)
    ax.set_xlabel('{} mag'.format(filters[1]))
    ax.set_ylabel('$\sigma$ {}'.format(filters[1]))

    plt.show()

    return error_key, error_array, scatter_array
