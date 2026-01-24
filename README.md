
# dsMatchingClient

[![codecov](https://codecov.io/gh/roygusinow/dsMatchingClient/graph/badge.svg?token=19CXOIZ0I9)](https://codecov.io/gh/roygusinow/dsMatchingClient)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![License:
MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

------------------------------------------------------------------------

<figure>
<img
src="https://raw.githubusercontent.com/roygusinow/dsMatchingClient/master/vignettes/resources/images/matching_noise.png"
alt="Differentially Private Propensity Score Matching" />
<figcaption aria-hidden="true">Differentially Private Propensity Score
Matching</figcaption>
</figure>

## Overview

[dsMatchingClient](https://github.com/roygusinow/dsMatchingClient)
implements privacy-preserving covariate balancing and causal effect
estimation for observational studies conducted in federated,
non-poolable data environments. It is designed for use within the
DataSHIELD infrastructure and provides client-side orchestration and
aggregation for federated matching and subclassification workflows.

The package accompanies the server-side package
[dsMatching](https://github.com/roygusinow/dsMatching), which performs
all operations involving individual-level data. Together, the two
packages enable causal analyses across multiple data providers without
sharing patient-level records.

## Functionality

- Federated propensity score subclassification and matching  
- Differentially private, non-identifying distance evaluation for
  matching  
- Estimation of ATE, ATT, ATU, risk ratios, and odds ratios  
- Cluster-robust sandwich standard errors via federated aggregation  
- Federated balance diagnostics and privacy-preserving plots (Love
  plots, QQ, eCDF)

All matched sets, subclass assignments, and outcomes remain local to
each server. For comparisons to the central estimation approach, please
see the
[methods](https://github.com/roygusinow/dsMatchingClient/tree/master/methods)
folder and unit test cases

## Quick Start

To get started with `dsMatchingClient`, you will need to have a working
DataSHIELD environment with the `dsMatching` package installed on the
server side. Try exploring the package
[vignette](https://github.com/roygusinow/dsMatchingClient/tree/master/vignettes)
for detailed examples and workflows.

## Installation

``` r
# Install the development version from GitHub
# install.packages("remotes")
remotes::install_github("roygusinow/dsMatchingClient")
```

Server-side nodes must also install:

``` r
remotes::install_github("roygusinow/dsMatching")
```
