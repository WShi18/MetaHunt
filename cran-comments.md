## Test environments

* Local: macOS Tahoe 26.2, R 4.6.0 (arm64)
* GitHub Actions: ubuntu-latest, windows-latest, macOS-latest (R-release and R-devel)
* win-builder: R-devel, R-release, R-oldrelease — all clean
* R-hub v2: ubuntu R-devel, atlas (Fedora-clang), donttest — all clean

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

* The "Possibly misspelled words" flagged by the CRAN incoming-feasibility
  check are author surnames (Imai, Shi, Zhang) and proper nouns / algorithm
  names (MetaHunt, d-fSPA, denoised) used throughout the methodology paper
  cited in the Description field.

## Downstream dependencies

There are currently no downstream dependencies for this package.
