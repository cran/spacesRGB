# NEWS for **spacesRGB** package


### Version 1.7-0  [2025-01-16]

* logging is now done with the **logger** package, which is imported
* the initial dictionary of 8 color spaces is now computed during `.onLoad()`, and `sysdata.rda` has been removed


### Version 1.6-1  [2024-11-16]

* added functions `RGBfromLab()` and `LabfromRGB()`
* improved accuracy for the breakpoint and slope of the linear part of the sRGB transfer function, and its inverse
* improved some man pages


### Version 1.5-0  [2024-01-24]

* improved clarity of the User Guide
* improved some function man pages
* removed `exportClasses` directive


### Version 1.4-0  [2021-12-06]

* fixed `bibliography.bib` to be compatible with `pandoc` v. 2.16.2
* fixed some stale URLs in `bibliography.bib` and man pages


### Version 1.3-0  [2019-12-10]

* fixed tolerance for ATLAS alternative BLAS/LAPACK implementation


### Version 1.2-2  [2019-01-30]

* added partial support for ACES Color workflows
* added S3 class `TransferFunction` with many associated methods
* added a User Guide


### Version 1.1-1  [2018-07-19]

* added 'scene RGB' - linear color data
* fixed some conceptual and terminology errors
* renamed functions and altered argument list appropriately
* allow non-trivial OOTF
* added a few more RGB spaces
* added function `plotPatchesRGB()`


### Version 1.0-4  [2018-05-31]

initial version on CRAN
