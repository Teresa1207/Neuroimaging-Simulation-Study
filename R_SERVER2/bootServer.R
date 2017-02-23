##########################################
############ Mediation Models ############
##########################################
### Create bootstroop sample data ###
library(nlme)
library(plyr)
source("medFunctions.R")
B = 1000
at = 1

dataB = readRDS('dataB.RDA')
source("bootMed.R")
atrophy = c("amygst","entorst","hippst","insst",
            "ltemst","mtemst","phippst","pcingst","slfst","splenst","thalst")
cog = c("memst","execst")
dx = list(NL = 'NL',EMCI = 'EMCI',LMCI = 'LMCI',TOTAL  = c('NL','EMCI','LMCI'))

bootRun = bootMed(B = B,atVar = atrophy[at],
                  cogVar = cog, dxVar = dx)

saveRDS(bootRun,file = paste(atrophy[at],'_bootRun.RDA',sep =''))


