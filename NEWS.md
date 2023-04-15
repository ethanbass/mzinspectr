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
