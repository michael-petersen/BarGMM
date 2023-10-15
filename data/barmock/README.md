### Bar Mock

This mock file is the bar 'growth phase' snapshot from the cusp halo model of Petersen, Weinberg, Katz (2021)[^1]. This test measures our ability to recover known self-consistent components -- and not identify any extra components (e.g., a knot). Given that we know the salient values for the bar and untrapped component, can we recover those from our GMM?

 - To convert from virial units, we assume a virial radius of $r_{\rm vir}$=220 kpc, and a virial velocity of $v_{\rm vir}$=150 km/s. This corresponds to a frequency scale of omegafac = $v_{\rm vir}/r_{\rm vir}$= 150/220 = 0.68 (i.e. multiply the virial frequencies by this value), and an angular momentum scale of Lfac = $v_{\rm vir}\times r_{\rm vir} = 150\times220 = 3.3e5$ kpc km/s.
 - The pattern speed of the model is 37.5 in virial units, or 25.5 km/s/kpc. This is modestly lower than estimates for the MW.
 - To estimate bar membership, I took 10 snapshots on either side of the exact time and used the classification from PWK21 to provide an estimate from 0-20 of whether the orbit is a member of the bar. A reasonable threshold might be 10.  
 - The bar was rotated by -12 degrees. This corresponds to an expected alpha fit of 78 degrees. We also expect then to see much larger x velocity dispersions compared to y velocity dispersions.
 - The trapped component of the bar ends at approximately 3 kpc in physical units (though the apparent length of the bar is much longer, see Petersen, Weinberg, Katz [2023]).
 - We cannot reliably classify any individual star (that is, we cannot provide a $p>0.8$ estimate for most stars). Instead, our validation asks how well we can recover the mean of the distributions.

![Figure 1](../../figures/angularmomentum_diagnostic.png "Lz Diagnostic"), shows angular momentum vs the square of the radius. The model for the bar is shown in black, while the angular momentum of the circular orbit at each radius is shown in blue. We estimated the angular momentum two different ways: the solid curve extracts the true potential in the $z=0$ plane to measure the circular velocity, while the dashed curve is an estimate for the circular velocity from a spherical mass enclosed measurement.

[^1]: For internal consumption, this is snapshot 900 from the run001 simulation.
