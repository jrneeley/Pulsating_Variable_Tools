# Variable-Star-Analysis
Routines for the identification, classification, and simulation of variable stars. 


![Classification example](classification_example.png?raw=true "variable_classification.py window")

## Classification

variable_classification.py runs an interactive window (shown above) that allows you to quickly classify candidate variable stars 
into their respective types (e.g. Classical Cepheids, RR Lyrae stars, etc.). The program plots the light curve (phased and raw) as
well as the color-magnitude diagram, period-luminosity relation, period-color relation, period-amplitude relation, and a postage stamp 
of the star on the image to help with classification. The user chooses the appropriate classification using the buttons on the right panel,
and the text file with the list of candidates is automatically updated. 

## Identification 
variable_identification.py contains useful functions for the automated detection of potential variable stars, including the calculation of 
several variability indices. 

## Simulations
This repository also contains code to generate simulated variable stars that mimic the properties of your dataset. These simulations are 
used to understand potential bias (such as incompleteness) stemming from the limits of your observations. They are also used to estimate 
realistic uncertainties on the periods, mean magnitudes, and amplitudes derived for your dataset. 

See the notebook Create_Simulated_stars.ipynb for an example of how to generate a set of simulated variable stars in the galaxy M33. 

Full API to come....
