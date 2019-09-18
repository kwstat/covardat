# alz.r
# Time-stamp: c:/x/rpack/covardat/data-raw/alz.r

# Note, this is NOT the same data as in the books
# Givens & Hoeting, Computational Statistics
# Der & Everitt, Handbook of Statistical Analyses using SAS

# Source: http://ftp.sas.com/samples/A58086

library(asreml)
library(kw)
library(Hmisc)
library(lattice)
library(rio)

datw <- read.table("c:/x/rpack/covardat/data-raw/unialz.txt",
                   header=TRUE, sep="", na.strings=".")
datw$trt <- factor(datw$trt, labels=c("High", "Low", "Placebo"))
head(datw)
require(reshape2)
dat <- melt(datw, id=c("trt","pat"))
names(dat) <- c("trt", "patient", "month", "score")
dat$trt <- relevel(dat$trt,"Placebo")
dat <- transform(dat, patient=factor(patient))
dat <- transform(dat, month=substring(month,5))
dat <- transform(dat, month=as.numeric(month))
head(dat)
library(rio)
export(dat, "c:/x/rpack/covardat/data/alzheimers.txt")
