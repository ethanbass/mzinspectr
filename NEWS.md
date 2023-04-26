# msdialreadr 0.3.1 dev

* Made rcdk (for rendering structures in `ms_search_gadget`) suggested instead of
a dependency.
* Simplified `ms_read_alignment` function and made corrections so that correct
types are read in for peak metadata fields. Missing retention indices are now
correctly entered as NAs.
* Added error message to `search_spectra` when peak retention indices are not
present and added option to search entire database by setting `ri_thresh` to `NULL`.
* Corrected normalizatio functions so they now return the peak table as a
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
