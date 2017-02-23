### EVAN FLETCHER
### MEDIATION ANALYSIS
library(nlme)
library(plyr)
fdata = readRDS(file = 'fdata.RDA')
fdata_bl = fdata[which(fdata$VISCODE2=='bl'),]
set.seed(12071983)
### Bootstrap Setup ###
nl = nrow(fdata_bl[which(fdata_bl$DX=='NL'),])
emci = nrow(fdata_bl[which(fdata_bl$DX=='EMCI'),])
lmci = nrow(fdata_bl[which(fdata_bl$DX=='LMCI'),])
### Create Data
NL.B = data.frame(DX = 'NL',ID = 1:nl,
                  sapply(1:B,function(B){sample(fdata_bl$RID[which(fdata_bl$DX=='NL')],
                                                nl,replace = TRUE)}))
EMCI.B = data.frame(DX = 'EMCI',ID = nl+1:emci,
                    sapply(1:B,function(B){sample(fdata_bl$RID[which(fdata_bl$DX=='EMCI')],
                                                  emci,replace = TRUE)}))
LMCI.B = data.frame(DX = 'LMCI',ID = (nl+emci)+1:lmci,
                    sapply(1:B,function(B){sample(fdata_bl$RID[which(fdata_bl$DX=='LMCI')],
                                                  lmci,replace = TRUE)}))
TOTAL.B = rbind(NL.B,EMCI.B,LMCI.B)




dataB = lapply(1:B,function(x){medBoot(B =x,data = fdata,dBoot = TOTAL.B)})
saveRDS(dataB,file = 'dataB.RDA')