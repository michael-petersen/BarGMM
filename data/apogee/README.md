## APOGEE data

This data has been pre-processed to match the APOGEE catalog against
 - Gaia DR3
 - Distance catalogues (Bailer-Jones, AstroNN, StarHorse)

The data files themselves are not included at present (they are relatively large for GitHub), but are available upon request.

### Validation Notes

The script `apogeevisualization.py` creates a view of the nine radial bins that are fit to the data.
For each radial bin, we show the $L_x$, $L_y$, and $L_z$ histograms of stars in the bin.
We also overlay the projected Gaussian fits to the data as blue curves, and the sum of the Gaussians as a grey curve.
The Gaussian fits are projected because the $x$ and $y$ Gaussians are allowed to rotate, which creates covariance between the dimensions.
The curves are all 500 realisations from the posterior chains, which gives a measure of uncertainty on individual components.
The effect of non-Gaussianity is visible in the quality of the fits in some bins; in particular the large amplitude of the $L_y$ component at smaller radii.
In general the fits are higher quality where the angular momentum uncertainties are smaller (larger Galactocentric radii).
However, all fits appear to be tolerable approximations for the distributions.

Angular momentum is a good proxy for orbital content. 
Particularly at a given radius, different orbital families will live in different (distinct) regions of angular momentum space.
If the potential is changing (e.g. from the effect of a rotating bar), these families may be smeared in angular momentum space around a locus.
Families that form a sequence (e.g. $x_1$ orbits that comprise the bar, or near-circular orbits) will have a smoothly changing locus in angular momentum space.
We always expect $\langle L_x \rangle = \langle L_y\rangle = 0$ for planar orbits.
For sufficiently different orbital families, we expect that angular momentum can provide a path to differentiate types of orbits. 
For the purposes of this work, we are interested in decomposing the inner galaxy into disc-like and bar-like orbits.
This decomposition revealed a third, unanticipated, distinct orbital family: with $\langle L_z \rangle$, independent of radius.
Such a population suggests a spherical distribution (i.e. with no preferential axis).

We show the fit validations in the figure below.
![Check angular momentum fits.](../../figures/fitcheck.png "Angular Momemtum reconstruction")