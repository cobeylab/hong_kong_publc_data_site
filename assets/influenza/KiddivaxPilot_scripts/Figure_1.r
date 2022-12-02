#
# R syntax to reproduce information for Figure 1 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective efficacy of seasonal influenza vaccination against seasonal and
# pandemic influenza virus infection during 2009 in Hong Kong
# CID, 2010 (in press).
#
# Last updated by Fang VJ, Ng S and Cowling BJ.
# November 8, 2010

dir <- "http://sph.hku.hk/data/KiddivaxPilotV1/"
hchar <- read.csv(paste(dir, "housechar_h.csv", sep=""))
random <- read.csv(paste(dir, "randomcode_h.csv", sep=""))

# randomized households
nrow(hchar)

# Allocation
hchar <- merge(hchar,random,all.x=T)
table(hchar$intervention)

# Analysis
vax <- hchar[hchar$intervention=="TIV",]
placebo <- hchar[hchar$intervention=="placebo",]
c(mean(vax$familysize),min(vax$familysize),max(vax$familysize))
c(mean(placebo$familysize),min(placebo$familysize),max(placebo$familysize))

#
# End of script.
#

