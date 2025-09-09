# fpcountr v1.0.1

* Prepped package for CRAN submission
* Updated installation instructions
* Added bug reports link

# fpcountr v1.0.0

### New functionality

* Added ability to process single timepoint data to `process_plate()`, `calc_fppercell()` and `calc_fpconc()` with `timecourse = FALSE` argument. These produce 96-well heatmap plots in place of timecourse plots.
* Added ability to process fluorescence data without OD (cell free fluorescence data) to `process_plate()`. This is triggered by `od_name = NULL`.
* Added parsers for Tecan Spark plate reader running SparkControl software: `parse_sparkcontrol()` for standard/endpoint or timecourse/kinetic data, and `parse_sparkcontrol_spectrum()` for spectrum data.
* Added new plate formatting functions: `find_wells()`, `find_wells_in_row()`, `find_wells_in_column()` and `find_wells_in_rowcolumn()`. These enable plate reader experiment analysis using plates of any size from 6- to 384-well sizes. The functions can also be used as helpers as a shorthand for referring to all wells in a given row, e.g. for specifying blank or empty wells in `process_plate()`, `calc_fppercell()` or `calc_fpconc()`. They also help behind the scenes to fit plots produced by `process_plate()`, `calc_fppercell()` and `calc_fpconc()` to the plate size automatically.
* Added new function to generate metadata template files, `save_metadata_template()`.

### Updated functions

* Standardised function names by: 
  * adopting verb-first function names, i.e. `parse_[instrument_software_name]()`, e.g. `magellan_parse()` -> `parse_magellan()`, `magellan_spectrum_parse()` -> `parse_magellan_spectrum()`
  * clarifying others e.g. `get_properties()` -> `get_fpbase_properties()`, `plot_absorbance_spectrum()` -> `process_absorbance_spectrum()`
  * standardising them as all lowercase e.g. `get_conc_A280()` -> `get_conc_a280()`, `get_conc_ECmax()` -> `get_conc_ecmax()`.
* Standardised argument names: `layout_csv` updated to `metadata_csv` in parsers and `generate_cfs()`.
* Standardised saved file names from `process_plate()`, `calc_fppercell()` and `calc_fpconc()`.
* Standardised plot styles in `process_plate()` and `calc_fppercell()`.
* Corrected plot labels in `process_plate()`.
* Removed hardcoded requirement in certain functions (`process_absorbance_spectrum()`, `get_conc_ecmax()`, `generate_cfs()`) for specific column names.
* Removed the ability of any function to silently save files. Most functions will now only save files if `outfolder` is specified as a valid path. Parsers only save files if `save_file` is set to `TRUE`.

### Breaking changes

* Parsing functions now only save files if `save_file` is set to `TRUE`. Without this, they just return the parsed data as a dataframe.

### New and updated vignettes

* The **Get Started** vignette has been updated for clarity. Every relevant function and argument is now explained in the plate reader data parsing-calibration-quantification workflow, along with the details of all output CSV files. Updated plots at the end to clarify information gained from calibration. Code updated to improve concentration estimation in example.
* **Data Parsing and Metadata** vignette explains what the data parsing functions do, and the metadata they expect. It also shows you other ways to parse your data if the package functions don't cover your plate reader.
* **Path Lengths** vignette describes why path length determination is important for calculating absolute concentrations from absorbance data, and the options for path length calculation within `plot_absorbance_spectrum()`.
* **Experiments I: Timecourse cellular data** describes how to use the experimental data processing functions (`process_data()` etc.) for data collected from microbial growth-linked fluorescent protein expression.
* **Experiments II: Endpoint data** describes how to use the experimental data processing functions (`process_data()` etc.) for data collected at a single timepoint, and the different types of plots it produces.
* **Experiments III: Cell-free data** describes how to use the experimental data processing functions (`process_data()` etc.) for fluorescence data collected without OD readings (e.g. for cell-free translation or purified FP quantification).
* **Plate formatting functions** vignette describes the use of functions such as `find_wells_in_column()`.

### Updated documentation

* Updated package Title and Description.
* Updated README.
* Added CITATION file.
* Updated function documentation.
* Updated website to Bootstrap 5, and to include new functions in Reference list.

### Minor bug fixes and improvements

* Numerous minor bug fixes.
* The OD name is now case insensitive in `process_plate()`.
* Improved error checks in `process_plate()`, `get_properties()` and `get_conc_ECmax()`.
* Fixed function examples in `process_plate()` and `calc_fppercell()`.
* Updated `parse_magellan()` to follow tidyverse style code.
* Compressed image files on website to reduce package size.

# fpcountr v0.2.0

- Improvements to functionality:
  - `get_conc_ECmax()` updated to correctly display plots and to save plots that report on success of normalisation.
  - fpbase data updated, alongside `get_properties()` data message
- Documentation updates:
  - README introduction expanded
  - README protocols.io link updated
  - README updated with citation links, badges and logo

# fpcountr v0.1.0

- Initial release.
