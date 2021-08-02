# Variable-Star-Analysis
Routines for the identification, classification, and simulation of variable stars (Neeley et al. 2021). This software writes and reads light 
curve files that are compatible with FITLC (Sarajedini et al. 2009), but can be modified fairly easily to work with other 
period fitting programs. In order to use this software without modification, you need to replicate the directory structure 
and files seen in test_data/vv124. 


![Classification example](classification_example.png?raw=true "variable_classification.py window")

## Identification 
`variable_identification.py` contains useful functions for the automated detection of potential variable stars, including the calculation of 
several variability indices. 

## Classification

`variable_classification.py` Runs an interactive window (shown above) that allows you to quickly classify candidate variable stars 
into their respective types (e.g. Classical Cepheids, RR Lyrae stars, etc.). The program plots the light curve (phased and raw) as
well as the color-magnitude diagram, period-luminosity relation, period-color relation, period-amplitude relation, period-amplitude ratio, 
and a postage stamp of the star on the image to help with classification. The user chooses the appropriate classification using the buttons 
on the right panel, and the text file with the list of candidates is automatically updated. We have included test_data 

**BEFORE RUNNING**: Create a file called `fit_results.txt` that contains the best fit parameters for all of your candidate variable stars.
See file in test_data folder for an example, but columns should be: 
1. candidate number (must also be the number in candidate light curve file name) 
2. variable id number in catalog (e.g. DAOPHOT ID number) 
3. X coordinate in pixels on master image 
4. Y coordinate in pixels on master image
5. Variable template number 
6. period in days 
7. epoch of maximum in MJD 
8. mean magnitude in bluest band 
9. amplitude in bluest band 
10. scatter around fit in bluest band 
11. mean magnitude in reddest band 
12. amplitude in reddest band 
13. scatter around fit in reddest band 
14. flag (should be set to -1 before classification) 
15. variable type (should be set to XXX before classification) 
16. variable subtype/pulsation mode (should be set to XX before classification) 

Copy the file `config_example.py` to `config.py`, and update the relevant fields for your system / the data you are working on.

**RUNNING** 
Classification program has 5 operating modes. 
mode | Purpose
------------ | -------------
first pass | Use for going through your candidates for the first time. Loops through all stars in fit_results.txt with flag = -1 
revise | Use to verify previous classifications. Loops through all stars in fit_results.txt with type != NV 
simulations | Use to classify simulated stars. No image will be shown, but all other diagnositic plots available. 
EB | Use to classify eclipsing binary stars (they use a different set of templates) 
specific star | Use to look at one specific star. Give a second argument with the candidate number (col 1 in fit_results.txt). 


To start the classification script, run the .py file and the mode you want to use. For example:
`python variable_classification.py 'first pass'` 
For the specific star mode, give both the mode and the star number you want: 
`python variable_classification.py 'star' 123`



## Simulations
`variable_simulations.py` This file contains functions to generate simulated variable stars that mimic the properties of your dataset. These simulations are 
used to understand potential bias (such as incompleteness) stemming from the limits of your observations. They are also used to estimate 
realistic uncertainties on the periods, mean magnitudes, and amplitudes derived for your dataset. In order to run these simulations, you 
must have the OGLEIV catalog of variable stars downloaded locally (http://ogle.astrouw.edu.pl). 

See the jupyter notebook included in this repository for a detailed example of how to generate a set of simulated variable stars for our test data. 
