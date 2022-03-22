# bistablehistory 1.0.0
## First CRAN Release
* Initial CRAN Release

# bistablehistory 1.1.0
## Improvements
* Custom prior values for history parameters, intercept terms, history effect, and fixed effects.
* Simplified Stan code.
* predict() computes values from history, reducing fit object size.
* predict() returns a vector of length that matches original table.

## Bug Fixes
* Change to difference in history values instead of the weighted mean.
* Use for scale instead of rate for prediction for Gamma family.
* Spelling in documentation.
* Additional tests.

# bistablehistory 1.1.1
## Improvements
* Posterior samples for cumulative history
