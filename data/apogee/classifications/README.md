### Files of classifications for the APOGEE data

 - All files start with `3Component_AllFeHCutMembership_Percentiles_reduceSNR_`.
 - There are two sets of unique bins: [[0,1],[1,2],[2,3],[3,4]] and [[0.5,1.5],[1.5,2.5],[2.5,3.5],[3.5,4.5]].
 - The integer bins are labelled as `r{}R{}_cyl.csv`, where the `{}` are replaced by the minimum and maximum radii, respectively.
 - The half-integer bins follow the same nomenclature, except the minimum and maximum radii are multiplied by 10. That is, for the bin that is [1.5,2.5], the filename is `3Component_AllFeHCutMembership_Percentiles_reduceSNR_r15R25_cyl.csv`.

For different bins, the mean radii and star numbers are given below.

| minimum radius | maximum radius | median radius | number of stars     | mean Lz uncertainty |
| -------------- | -------------- | ------------- | ------------------- | ------------------  |
| 0.0            | 1.0            | 0.708         | 1867                | 313.6               |
| 0.5            | 1.5            | 0.992         | 2677                | 280.1               |
| 1.0            | 2.0            | 1.536         | 2894                | 229.1               |
| 1.5            | 2.5            | 1.992         | 3108                | 188.9               |
| 2.0            | 3.0            | 2.541         | 3289                | 164.1               |
| 2.5            | 3.5            | 3.026         | 3695                | 135.9               |
| 3.0            | 4.0            | 3.580         | 4739                | 95.9                |
| 3.5            | 4.5            | 4.033         | 5930                | 75.0                |
