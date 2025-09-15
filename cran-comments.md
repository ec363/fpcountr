## Resubmission
This is a resubmission. In this version I have:

* Fixed the issue of the example in `get_fpbase_properties()` function running too slowly on CRAN.
* Re-checked packaged with `devtools::check(remote = TRUE)`.

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

❯ checking CRAN incoming feasibility ... [3s/20s] NOTE
New submission

Size of tarball: 7730891 bytes

❯ checking installed package size ... NOTE
installed size is  7.2Mb
sub-directories of 1Mb or more:
  doc    4.1Mb
help   2.4Mb

Reasons for package size explained below.

## Resubmission
This is a resubmission. In this version I have:

* Put single quotes around software names in the Description.

Response to comment on size: The size is 7MB+ primarily due to documentation (help and doc files). Vignettes and function docs contain extensive guidance for new users, including users new to R. This includes datasets and figures to illustrate what the package can do. The figures have been compressed. I believe the guidance is necessary to support the usability of this package for the intended audience.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

Installed size of package is 7.2Mb, primarily due to documentation (help and doc files). Vignettes and function docs contain extensive guidance for new users, including users new to R. This includes datasets and figures to illustrate what the package can do. The figures have been compressed. I believe this will support the reproducibility and usability of this package.
