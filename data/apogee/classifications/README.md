### Files of classifications for the APOGEE data

 - All classifications have now been put into a single file, `AllClassifications.h5`.
 - For each star, there are now a list of 20 different classifications, each from a different posterior chain.
 - This file must be read with `h5py`.
 - When opening, each star is a `key` in the dataset, indexed by APOGEE ID. If a classification comes from a half-integer bin, the APOGEE ID has a `*` appended to the end.
 - Each dataset contains a 20x4 array. `[:,0]` gives all disc probabilities, `[:,1]` gives all bar probabilities, `[:,2]` gives all knot probabilites. The last row of the array is a diagnostic row, with `[R,x,y,z,Lx,Ly,Lz]` (not needed).
 - As an example to read and compute the median bar probability given the 20 classifications:
 ```
 import h5py
 f = h5py.File("data/apogee/classifications/AllClassifications.h5","r")
 # f.keys() # this would list all the APOGEE IDs, including doubles for half-integer bins
 median_bar_probability = []
 for key in f.keys():
   median_bar_probability.append(np.nanmedian(f[key][:,1]))

f.close() # don't forget to close the file!
```
- Analogous calls could get the median disc or knot probabilities.
- This loop setup isn't exactly the fastest, as the data structure creates some overhead, but it is the most compact way to represent the data.

*For FEH8_SNRC classifications*
For different bins, the mean radii and star numbers are given below.

| minimum radius | maximum radius | median radius | number of stars     | mean Lz uncertainty |
| -------------- | -------------- | ------------- | ------------------- | ------------------  |
| 0.0            | 1.0            | 0.699         | 1429                | 115.6               |
| 0.5            | 1.5            | 1.017         | 2061                | 112.3               |
| 1.0            | 2.0            | 1.525         | 2245                | 112.8               |
| 1.5            | 2.5            | 2.035         | 2522                | 114.5               |
| 2.0            | 3.0            | 2.572         | 3070                | 115.5               |
| 2.5            | 3.5            | 3.049         | 3850                | 113.9               |
| 3.0            | 4.0            | 3.556         | 4848                | 106.7               |
| 3.5            | 4.5            | 4.033         | 5785                |  98.9               |
| 4.0            | 5.0            | 4.545         | 6679                |  88.5               |


*For FEH6_SNRC classifications*
For different bins, the mean radii and star numbers are given below.

| minimum radius | maximum radius | median radius | number of stars     | mean Lz uncertainty |
| -------------- | -------------- | ------------- | ------------------- | ------------------  |
| 0.0            | 1.0            | 0.708         | 2196                | 136.6               |
| 0.5            | 1.5            | 1.000         | 3206                | 133.9               |
| 1.0            | 2.0            | 1.541         | 3537                | 132.4               |
| 1.5            | 2.5            | 1.996         | 3839                | 131.5               |
| 2.0            | 3.0            | 2.530         | 4007                | 124.8               |
| 2.5            | 3.5            | 3.021         | 4398                | 117.6               |
| 3.0            | 4.0            | 3.571         | 5493                | 109.1               |
| 3.5            | 4.5            | 4.037         | 6812                | 101.5               |
