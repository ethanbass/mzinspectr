# mzinspectr
<!-- badges: start -->
  [![mzinspectr status badge](https://ethanbass.r-universe.dev/badges/mzinspectr)](https://ethanbass.r-universe.dev)
  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10426253.svg)](https://doi.org/10.5281/zenodo.10426253)
  [![stability-experimental](https://img.shields.io/badge/stability-experimental-orange.svg)](https://github.com/emersion/stability-badges#experimental)
<!-- badges: end -->
   
A few functions for analyzing mass spectrometry alignments in R. Currently, MS-DIAL `.txt` alignment files are supported.

### Installation

It's recommended to install mzinspectr directly from GitHub.

```
install.packages("devtools")
devtools::install_github("https://github.com/ethanbass/mzinspectr/")
```

It can also be installed via [R Universe](https://ethanbass.r-universe.dev/):

```
install.packages("mzinspectr", repos="https://ethanbass.r-universe.dev/", type="source")
```

### Exporting alignments from MS-DIAL.

Export your alignment from MS-DIAL by selecting the `Export:Alignment Result` menu. Select Area or Height as appropriate.

### Reading alignment into R

Read your alignment into R using the `ms_read_alignment` function, providing the path to your MS-DIAL alignment file. This will produce a list of matrices. The first element will be the actual peak table (`tab`) with the areas or heights of the features detected by MS-DIAL. The second element will contain the peak metadata returned by MS-DIAL such as retention indices, mass spectra, and any identifications done in MS-DIAL (`peak_meta`). The third element is intended to hold sample metadata (`sample_meta`) provided by the user, which can be added using the `ms_attach_metadata` function.

The `ms_attach_metadata` function takes three arguments: an MS-DIAL alignment object (`x`), a dataframe or matrix containing the sample metadata (`meta`), and a string specifying the column in `meta` to be matched with the names of the samples (`col`).

###  Normalization

There are currently several different options for feature normalization, including normalization by an internal standard (`ms_normalize_itsd`), total sum normalization (`ms_normalize_tsn`) and probabilistic quotient normalization (`ms_normalize_pqn`). You can also subtract the mean or median value from a set of blanks from each peak using the `ms_subtract_blanks` function.

### Peak identification

There is preliminary support for peak identification by searching a user-provided mass-spectral database through the `ms_search_spectra` function. The database can be loaded into R from an MSP file using the It takes several parameters, including an ms-dial alignment object (`x`), a spectral database (`db`), the column or columns to identify (`cols`), the maximum retention index shift to exclude a match from consideration (`ri_thresh`), the relative weight to give spectral similarity versus retention index similarity (`spectral_weight`), the number of results to return (`n_results`), and the number of cores to use for parallel processing (`mc.cores`).

To compile a mass spectral database, I recommend using [mspcompiler](https://github.com/QizhiSu/mspcompiler).

### Visualization

Spectra can be plotted using either "base R" graphics or "plotly" graphics using the `ms_plot_spectrum` function, which takes an ms-dial alignment object (`x`) and a column index (`col`) as argument.

### Citation

If you use mzinspectr in published work, please cite it as follows:

Bass, E. (2023). mzinspectr: Read and Analyze Mass Spectrometry Alignment Files (version 0.4.0). https://doi.org/10.5281/zenodo.10426253.


