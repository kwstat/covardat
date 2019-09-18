# mouse.r

The data originally appear in :
A Distribution-Free Test for Tumor-Growth Curve Analyses with Application to an Animal Tumor Immunotherapy Experiment
James A. Koziol, Donna A. Maxwell, Matsuro Fukushima, M. E. M. Colmerauer and Yosef H. Pilch
Biometrics
Vol. 37, No. 2 (Jun., 1981), pp. 383-390
http://www.jstor.org/stable/2530427

Also appear in Dawson 1997 (with different mouse numbering).
# ----------------------------------------------------------------------------

dat <- read.table("c:/x/rpack/covardat/data/koziol.mouse.txt", header=TRUE)

dat <- tumor
require(lattice)
xyplot(size~day, data=dat, groups=mouse, type='l')
xyplot(size~day|group, data=dat, groups=mouse, type='l',
       main="Mouse tumor size",
       layout=c(3,1))

require(reshape2)
d2 <- melt(dat, id.vars=c('group','mouse','day'))
d2 <- acast(d2, mouse~day)
splom(d2)
#pairs(scale(d2))
round(cor(scale(d2), use="pair"),2)

lib(dplyr)
d3 <- dat
d3 <- mutate(d3, logsize=log(size+1))
d3 <- group_by(d3, group, day)
d3 <- mutate(d3, ssize=scale(logsize, scale=FALSE)) %>% filter(!is.na(size))
d3 <- acast(melt(d3, measure.vars=c('ssize')), mouse~day)
#pairs(d3)
lib(corrgram)
# Doesn't look quite right
corrgram(d3, upper.panel=panel.pts, lower.panel=panel.conf)
