# BarGMM
### Dissecting the inner galaxy with Gaussian Mixture Models
#### MSP (started May 2023)

The main directory has three scripts to run classifications, and one example script:
 - `classifyapogee.py`: run the classification for APOGEE stars, print the table of the components.
 - `classifybarmock.py`: run the classification for a mock with a disc + bar, print the table of the components. This is described further in `models/barmock.py` and `data/barmock/README.md`.
 - `classifybulgemock.py`: run the classification for a mock with a disc + bar + bulge, print the table of the components. This is described further in `models/bulgemock.py` and `data/bulgemock/README.md`.
 - `readclassifications.py`: demonstration of how to read the classification files.


 Additionally, the `diagnostics` folder includes
 - `plotellipses.py`: demonstration of how to display ellipses from the fits.
 - `testclassification.py`: for the mocks, which have known memberships, validate the classifications to estimate per-star uncertainties.
