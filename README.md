# msdialreadr
<!-- badges: start -->
  [![msdialreadr status badge](https://ethanbass.r-universe.dev/badges/msdialreadr)](https://ethanbass.r-universe.dev)
  [![stability-experimental](https://img.shields.io/badge/stability-experimental-orange.svg)](https://github.com/emersion/stability-badges#experimental)
<!-- badges: end -->
   
A few functions for analyzing [MS-DIAL](http://prime.psc.riken.jp/compms/msdial/main.html) alignments in R.

**Note**: This package is not created or endorsed by the creators of MS-DIAL.

### Installation

Currently it's recommended to install msdialreadr directly from GitHub.

```
install.packages("devtools")
devtools::install_github("https://github.com/ethanbass/msdialreadr/")
```

It can also be installed via [R Universe](https://ethanbass.r-universe.dev/):

```
install.packages("msdialreadr", repos="https://ethanbass.r-universe.dev/", type="source")
```


### Exporting your alignment

Export your alignment from MS-DIAL by selecting the `Export:Alignment Result` menu.

### Reading alignment into R

Read your alignment into R using the `read_alignment` function, providing the path to your MS-DIAL alignment file. This will produce a list of matrices. The first element will be the actual peak table (`tab`) with the areas or heights of the features detected by MS-DIAL. The second element will contain the peak metadata returned by MS-DIAL such as retention indices, mass spectra, and any identifications done in MS-DIAL (`peak_meta`). The third element is intended to hold sample metadata (`sample_meta`) provided by the user, which can be added using the `attach_metadata` function.

The `ms_attach_metadata` function takes three arguments: an MS-DIAL alignment object (`x`), a dataframe or matrix containing the sample metadata (`meta`), and a string specifying the column in `meta` to be matched with the names of the samples (`col`).

###  Normalization

There are currently several different options for feature normalization, including normalization by an internal standard (`ms_normalize_itsd`), total sum normalization (`ms_normalize_tsn`) and probabilistic quotient normalization (`ms_normalize_pqn`). You can also subtract the mean or median value from a set of blanks from each peak using the `ms_subtract_blanks` function.

### Peak identification

There is preliminary support for peak identification by searching a user-provided mass-spectral database through the `ms_search_spectra` function. It takes several parameters, including an ms-dial alignment object (`x`), a spectral database (`db`), the column or columns to identify (`cols`), the maximum retention index shift to exclude a match from consideration (`ri_thresh`), the relative weight to give spectral similarity versus retention index similarity (`spectral_weight`), the number of results to return (`n_results`), and the number of cores to use for parallel processing (`mc.cores`).

To compile a mass spectral database, I recommend using [mspcompiler](https://github.com/QizhiSu/mspcompiler).

### Visualization

Spectra can be plotted using either "base R" graphics or "plotly" graphics using the `ms_plot_spectrum` function, which takes an ms-dial alignment object (`x`) and a column index (`col`) as argument.

### Citation

If you use msdialreadr in published work, please cite it as follows:

Bass, E. (2023). msdialreadr: Read and Analyze MS-DIAL Alignment Files (version 0.3.0).
