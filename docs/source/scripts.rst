Scripts
=====

The scripts are divided into five sets:

 1. Common variables and functions
 2. Data setup
 3. Main analyses
 4. Simulation validation
 5. Plotting

Script names have a three-digit prefix indicating their set (first digit) and running order.

Within each set, the x00 and x01 script contain common objects and functions respectively. This helps keep the code separate and clean.

.. 1xx scripts: common
1xx scripts: common
------------

Common functions and a re-write of several gamlss functions (that had issues within our pipeline).

* 100.common-variables.r
* 101.common-functions.r
    * **`Create.Folders()`:** Create folder structure for a specific subset
    * **`Check.Attributes()`:** Utility function for subsets of data
* 102.gamlss-recode.r
    * **`bfpNA()`:** Custom version of `bfp()` that can handle `NA`s
    * **`GGalt()`:** Custom version of `GG()` (GAMLSS family) with robust `GGalt$mean()` and `GGalt$variance()`

These scripts are sourced in later scripts, they define common variables/objects and functions.

Importantly, they also include a re-write of several gamlss functions to address numerical instability (these may not be necessary in a novel replication, however they are required for using our output fitted objects).


.. 2xx scripts: data setup
2xx scripts: data setup
------------

Import and clean the data ready for the gamlss fitting. Also, generate simulated dataset (called omega) used for validation.

* 200.variables.r
* 201.functions.r
    * None
* 211.data-setup.r
    * Custom script to load and clean the raw data, outputs are `SUBSET.rds` and model objects
* 220.simulation-omega-setup.r
    * Custom script to create simulated data, outputs are `SUBSET.rds` and model objects

Each data script is a custom

.. 3xx scripts: 
3xx scripts: 
------------

These scripts perform the substantial calculations and model fitting.

* 300.variables.r
* 301.functions.r
    * **`Fit.Function()`:** Calls `Extract.Wrapper()`
        * **`gamlssWrapper()`:** Simple wrapper around `gamlss()` to ensure consistent calls
        * **`Extract.Summary()`:** Generate consistent summary of subset/dataset
        * **`Extract.Param()`:** Create custom `ParamObj` (new class of object) from `gamlss()` output
        * **`Extract.Wrapper()`:** Combined call of `gamlssWrapper()`, `Extract.Summary()` and `Extract.Param()`
        * **`Save.Extracted()`:** Save `gamlss()` output objects, called from within `Fit.Function()`
    * **`Make.bfpNA.model.from.extract()`:** Convert `fp()` into `bfpNA()`, i.e. define fractional polynomial explicitly
    * **`Find.Models.To.Fit()`:** Find candidate models under MODEL folder
    * **`Make.Longitudinal()`:** For longitudinal follow-up, generate subject-specific summary measures, e.g. IQR
    * **`Boot.Function()`:** Generate a bootstrap replicate dataset and call `Extract.Wrapper()`
    * **`ValidateCleanInput()`:** Check dataset conforms to `ParamObj`
    * **`Apply.Param()`:** Use `ParamObj` to generate predictions for a dataset (calls `ValidateCleanInput()`)
    * **`Apply.FitAndBoot()`:** Calls `Apply.Param()` on all bootstrap replicate `ParamObj`s
    * **`Load.Subset.Wrapper()`:** Load multiple elements into `HOLDER` object
    * **`Calc.Expanded()`:** Wrapper calling `Ranef.MLE.Func()` and `Add.New.Ranefs()` 
        * **`Find.Fitted.Levels()`:** Compare `ParamObj` with dataset to find studies with fitted random-effects
        * **`Find.Missing.Levels()`:** Compare `ParamObj` with dataset to find studies with missing/unknown random-effects
        * **`Ranef.MLE.Func()`:** Estimate random-effects using maximum likelihood (using dXX from GAMLSS family)
        * **`Add.New.Ranefs()`:** Expand a `ParamObj` with new study random-effects
* 310.fitting.r
    * Uses `Fit.Function()`
* 320.best-fit.r
    * Extracts BIC values from fitted models and selects the best (makes a copy or symlink as `MODEL.rds`)
* 330.bootstrapping.r
    * Uses `Boot.Function()`
* 340.bootstrap-merge.r
    * Merges separate bootstrap outputs into `BOOT.EXTRACT.rds`
* 350.calc-derived.r
    * Uses `Apply.Param()` and `Apply.FitandBoot()` to create all derived curves and outputs (for lifespan and study-specific curves) saved as `DERIVED.rds`
* 350.calc-novel.r
    * Uses `Calc.Expanded()` to estimate random-effects for novel data saved as `FIT.EXPANDED.rds` (for fit and bootstrap replicates)

Main scripts, these fit the gamlss model(s), select the best (via
BIC), perform bootstrapping, and calculate all necessary derived values.


.. 5xx scripts: plotting
5xx scripts: plotting
------------

Some example plotting scripts using GNU R's base graphics. The article uses "nicer" plots generated using ggplot2 (not included in this repository).

* 500.plotting-variables.r
* 501.plotting-functions.r
    * None
* 510.plotting.r
    * Example plots using `DERIVED.rds` object

Plotting functions, these *only* use the `DERIVED.rds` and the fitted objects from novel data (`FIT.EXPANDED.rds`).
