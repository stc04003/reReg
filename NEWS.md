# `reReg` 1.4.3
	* Added an option to specify the heuristic constant (numAdj).
	* Move duplicated Rd files for Recur(), %2% (%to%), and mcf() into reexport.Rd
# `reReg` 1.4.2
  	* Added an option to perform GMM or EL procedure to improve cox.LWYY.
	* Added a vignette on CPPL.
g; typing ?Recur will now be directed to the help page resid in reda.
# `reReg` 1.4.1
  	* printCoefmat now all has the argument `has.Pvalue = TRUE`
	* Improved print methods for summary.reReg
	* Improved matrix algebra in calculation
	* Replace exceeded cat() and print() with message() and warning()
	* Export reReg.control
# `reReg` 1.4.0
  	* `simGSC()` now allows users to specifies design matrix and censoring distribution.
	* Some name changes (mostly arguments) to match with the JSS submission.
	* Change names "SC" to "GSC"; this includes the `model` argument and `simSC()` to `simGSC()`.
	* Added a basebind() function to combine baseline rate/hazard plots.
	* Temporary disabled sandwich variance estimators for further investigation.
	* Import `mcf()` from reda.
	* Event plot can now plot in calendar times.
	* Updated examples and online vignettes.
# `reReg` 1.3.1
  	* simSC() now returns times in t.start and t.stop.
	* Added origin in simSC()
	* Updated examples
	* Import summary.Recur from reda
# `reReg` 1.3.0
  	* Re-organized `reReg()`; it now provides general model assumptions.
	* Improved speed via RCpp
# `reReg` 1.2.0
  	* Fixed memory errors (checked with valgrind 3.15)
# `reReg` 1.2.0
  	* Adopt `Recur()` from package `reda`
	* Added a draft for regression vignettes
	* Changed function name `simDat` to `simSC`
	* Fixed bug with only 1 covariate
# `reReg` 1.1.6
	* Added sandwish variance estimations to most implementations
	* Added vignettes on simulation and plots
# `reReg` 1.1.5
	* Cleaned reSurv and am.GL codes
	* updated event plot, baseline function plots
	* added CMF plot
# `reReg` 1.1.4
	* Rebuild with oxygen
# `reReg` 1.1.3
	* Removed plotEvent function and add its features to plot.reSurv
	* Updated reSurv so it split out tibble df
# `reReg` 1.1.2
	* Fixed bugs reported in CRAN regarding r-devel
# `reReg` 1.1.0
	* Major make over
	* Cleaned codes
	* New data structure required for `reSurv`: Time, id, event, status
# `reReg` 1.0.0
	* Version updated to published version.   
 