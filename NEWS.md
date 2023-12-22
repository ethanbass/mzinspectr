# mzinspectr 0.4.0

* Changed name of package to "mzinspectr".
* Updated `ms_filter_alignment` so it filters matches when removing columns from the peak table.
* Added titles to mass spectra plots with peak name and retention index or time.
* Fixed bug in retention index matching algorithm for `ms_search_spectra` function.
* Added `fixed_levels` argument to `ms_reshape_peaktable` so features can be plotted in the order they're provided by the user.
* Changed the name of the `mzinspectr` alignment object from `msdial_alignment` to `ms_alignment`.

# msdialreadr 0.3.3

* Fixed bug `ms_tidy_msdial` when renaming peaks that are out of order.
* Renamed `ms_tidy_msdial` to `ms_reshape_peaktable`. The old function name is now deprecated.
* Made some improvements to documentation.

# msdialreadr 0.3.2

* Added option to for renaming peaks via `ms_tidy_msdial` by providing a named character vector.
* Deprecated `treatment` argument in `ms_tidy_msdial`, replacing it with new `metadata` argument.

# msdialreadr 0.3.1

* Made rcdk (for rendering structures in `ms_search_gadget`) suggested instead of
a dependency.
* Simplified `ms_read_alignment` function and made corrections so that correct
types are read in for peak metadata fields. Missing retention indices are now
correctly entered as NAs.
* Added error message to `search_spectra` when peak retention indices are not
present and added option to search entire database by setting `ri_thresh` to `NULL`.
* Corrected normalization functions so they now return the peak table as a
data.frame instead of a matrix.

# msdialreadr 0.3.0

* Added support for calling MSDIAL console app directly through the `ms_call_msdial` function.
* Improved speed of `ms_search_spectra` function and added progress bar.
* Removed `search_msp` from NAMESPACE.
* Added `ms_search_gadget` for interactively viewing `search_spectra` matches.
* Changed names of normalization functions: `tsn` to `ms_normalize_tsn` and
`pqn` to `ms_normalize_pqn`.
* Added `ms_rt_to_ri` function for conversion of retention times to retention indices.
* Added `ms_` prefix to all exported functions to prevent NAMESPACE conflicts.

# msdialreadr 0.2.0

* Update license to GPL3.
* Added a `NEWS.md` file to track changes to the package.
* Added `search_spectra` & `search_msp` functions for peak identification.
