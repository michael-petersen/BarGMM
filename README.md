# BarGMM
### Dissecting the inner galaxy with Gaussian Mixture Models
#### MSP May 2023

The analysis is currently set up for the data, which is coded in the file `models/apogee.py`. The variables in this file are passed to the other scripts via the `from models.apogee import *` line in
 - `classifystars.py`: run the classification, print the table of the components
 - `plotellipses.py`: demonstration of how to display ellipses from the fits
 - `readclassifications.py`: demonstration of how to read the classification files

You may also look at the mock files by changing the relevant model call to `from models.bulgemock import *`, which will import the variables from `models/bulgemock.py`.
