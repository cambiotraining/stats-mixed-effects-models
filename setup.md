---
title: "Data & Setup"
number-sections: false
---

<!-- 
Note for Training Developers:
We provide instructions for commonly-used software as commented sections below.
Uncomment the sections relevant for your materials, and add additional instructions where needed (e.g. specific packages used).
Note that we use tabsets to provide instructions for all three major operating systems.
-->

## Data

The data used in these materials are provided as a zip file.

First, create a folder on your computer for the course (e.g. `mixed-models`).
Then, download and unzip the the following data to the project folder you created.

<!-- Note for Training Developers: add the link to 'href' -->
<a href="https://github.com/cambiotraining/stats-mixed-effects-models/raw/main/materials/data/data.zip">
  <button class="btn"><i class="fa fa-download"></i> Download (Zip file)</button>
</a>

## R and RStudio

::: {.tabset group="os"}

#### Windows

Download and install all these using default options:

- [R](https://cran.r-project.org/bin/windows/base/release.html)
- [RTools](https://cran.r-project.org/bin/windows/Rtools/)
- [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

#### Mac OS

Download and install all these using default options:

- [R](https://cran.r-project.org/bin/macosx/)
- [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

#### Linux

- Go to the [R installation](https://cran.r-project.org/bin/linux/) folder and look at the instructions for your distribution.
- Download the [RStudio](https://www.rstudio.com/products/rstudio/download/#download) installer for your distribution and install it using your package manager.

:::

## R packages

From an R console, you can run the following command to install all the packages used in this course: 

```r
install.packages(c("broom", 
                   "corrr", 
                   "distributional", 
                   "downlit", 
                   "downloadthis", 
                   "ggdist", 
                   "ggResidpanel", 
                   "janitor", 
                   "kableExtra", 
                   "knitr", 
                   "lme4", 
                   "patchwork", 
                   "pwr", 
                   "reticulate", 
                   "rstatix", 
                   "tidyverse", 
                   "lmerTest", 
                   "broom.mixed", 
                   "performance", 
                   "see",
                   "glmmTMB", 
                   "brms", 
                   "MASS", 
                   "RLRsim", 
                   "pbkrtest",
                   "DHARMa"))
```