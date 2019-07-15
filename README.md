## CVEPJointModel: An R Package for Joint Modelling Procedure in Analyzing Highly Unbalanced CVEP Data

Authors: [Qian M. Zhou](http://qianmichellezhou.net), Jingjing Pan,

Vignette: Joint Modeling in Analyzing Highly Unbalanced Multi-Environment Trial Data (under review)

Please email all comments/questions to qz70@msstate.edu or panjj1125@outlook.com

### Installation Instructions

You can load the package and use the function install_github

```
library(devtools)
install_github("panjj1125/CVEPJointModel",dependencies=TRUE)
```

Note that this will install all the packages suggested and required to run our package.  It may take a few minutes the first time, but this only needs to be done on the first use.  In the future you can update to the most recent development version using the same code. 

### Getting Started
The main function to estimate the model is `CVEP_JM()` but there are a host of other useful functions. Here is one demo:

```
library(CVEPJointModel)
dat <- read.csv('3year_NCVT.csv')
Model <- CVEP_JM(dat, factors = c("Year","Loc","Rep", "Variety"),
                 TT_mm=list(fixed=c("Year","Loc"),
                 random=c("Variety","Variety:Year","Variety:Loc")),
                 DS_mm=list(fixed=c("Year","Loc"),random=c("Variety")),
                 converg_control=list(nsamp=5,max.iter=5,err=10^(-5),err1=10^(-5)))
```
When there is controls in fixed effects, we can use the following code:

```
Model <- CVEP_JM(dat, factors = c("Year","Loc","Rep", "Variety", "Checks"), TT_mm=list(fixed=c("Year","Loc", "Checks"),random=c("Variety","Variety:Year","Variety:Loc")),DS_mm=list(fixed=c("Year","Loc"),random=c("Variety")),converg_control=list(nsamp=5,max.iter=5,err=10^(-5),err1=10^(-5)))

```

To predict the random effects, we use
```
r <- raneff(Model)
```

### Citation information

Please cite the code using following formula:

    @misc{Qian2019,
        author = {Qian M. Zhou, Jingjing Pan},
        title = {CVEPJointModel: An R Package for Joint Modelling Procedure in Analyzing Highly Unbalanced CVEP Data},
        year = {2019},
        publisher = {GitHub},
        journal = {GitHub repository},
        howpublished = {\url{https://github.com/panjj1125/CVEPJointModel}},
    }
