# Height and Width of GUI window
HEIGHT = 800
WIDTH = 1500

# Where does your data live? (best to give absolute path in practice)
data_dir = 'test_data/vv124/'

# What is the file extension for your light curve and fit files?
lcv_suffix = '.fitlc'
fit_suffix = '.fitlc_fit'

# Do you want to plot the LMC PL relations as a guide for classification?
plot_lmc = True
# What is the approximate distance modulus for galaxy you are working on?
DM = 25.71
# What is the extinction in the two filters (bluest first)?
ext_mag1 = 0.052
ext_mag2 = 0.024
# Name of star catalog with xy positions
catalog = data_dir+'vv124_bi_cal.raw.tran.radec'
# What is the name of your montage image?
image = data_dir+'images/montage.fits'
# Does is have xy offsets from the catalog? (likely if created using Montage2)
montage_x_offset = -10
montage_y_offset = -12
# Define appropriate limits for your image
image_limits = [0,1000]

# Names of filters
filt1_name = 'F475W'
filt2_name = 'F814W'


# zp, slope, and scatter of period-color relation for cep and rrl in relevant bands
plot_pc = True
cep_pc_relation_coeffs = [0.685, 0.636, 0.1] # B-I relation from Sandage et al 2009
rrl_pc_relation_coeffs = [1.072, 1.325, 0.1] # B-I relation from ??
