### This file explains how to reproduce the figures of the article
# The files of the inst/reproduce folder contain the scripts to reproduce the figures.
#
# The files finishing in .run.R are the ones creating the result files;
# these scripts take several hours to run on a desktop computer.
#
# The files finishing in .plots.R create the plots, assuming that the .run.R files
# have been run already. Indeed the plot files require specific .RData files that are
# created by the .run.R files.

# folders (these should be changed to point to existing directories)
# scriptfolder <- "~/path to this package/inst/reproduce"
scriptfolder <- "~/Dropbox/Modularization/code/bettertogether/inst/reproduce"
# resultsfolder <- "~/path to a folder where large RData files will be created"
resultsfolder <- "~/Dropbox/ModularizationResults/test/"
#
setwd(resultsfolder)

## Biased data (warning: can take a few hours to run!)

# source(file.path(scriptfolder, "biaseddata.generatedata.run.R"))
# source(file.path(scriptfolder, "biaseddata.posterior1.run.R"))
# source(file.path(scriptfolder, "biaseddata.posterior2alone.run.R"))
# source(file.path(scriptfolder, "biaseddata.posterior1given2.run.R"))
# source(file.path(scriptfolder, "biaseddata.posterior2given1.run.R"))
source(file.path(scriptfolder, "biaseddata.posterior2given1plugin.run.R"))
source(file.path(scriptfolder, "biaseddata.posterior2given1cut.run.R"))
## plots and tables
source(file.path(scriptfolder, "biaseddata.plots.R"))

## Epidemiology data
# source(file.path(scriptfolder, "epidemiology.posterior1.run.R"))
# source(file.path(scriptfolder, "epidemiology.posterior2alone.run.R"))
# source(file.path(scriptfolder, "epidemiology.posterior1given2.run.R"))
# source(file.path(scriptfolder, "epidemiology.posterior2given1.run.R"))
# source(file.path(scriptfolder, "epidemiology.posterior2given1plugin.run.R"))
# source(file.path(scriptfolder, "epidemiology.posterior2given1cut.run.R"))
## plots and tables
# source(file.path(scriptfolder, "epidemiology.plots.R"))

