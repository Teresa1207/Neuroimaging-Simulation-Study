### EVAN ANALYSIS
setwd("~/Documents/Davis 2016-2017/ADC_Evan/R_FOLDER")
require(dplyr)
#load data
edata = read.csv('~/Documents/Davis 2016-2017/ADC_Evan/R_FOLDER/SequencingAnalysisRedo.csv',sep = ',',na.strings = "")
dandata = read.csv('~/Documents/Davis 2016-2017/ADC_Evan/R_FOLDER/adni_mediation.csv',sep = ',',na.strings = "")
#cog measures
cdata = read.csv('~/Documents/Davis 2016-2017/ADC_Evan/R_FOLDER/UWNPSYCHSUM_04_22_16.ccsv.csv',sep = ',',na.strings = "")
edata$ID = as.character(edata$PTID)

eRID = strsplit(edata$ID,"-")
edata$RID = unlist(lapply(eRID,function(x){as.integer(x[3])}))
cdata$edata = cdata$RID%in%edata$RID
edata = edata[order(edata$RID),]
vis = data.frame(vis = c("bl","m06","m12","m18","m24"),time = c(0,6,12,18,24))
cdata.sub = subset(cdata,edata==TRUE & VISCODE2%in%vis$vis)
cdata.sub = merge(cdata.sub,vis,by.x = 'VISCODE2',by.y = 'vis')
cdata.sub = cdata.sub[order(cdata.sub$RID),]

##Dates
cdata.sub$edate <- as.Date(as.character(cdata.sub$EXAMDATE),format="%m/%d/%Y")
cdata.sub$edate <- as.Date(as.character(cdata.sub$EXAMDATE),format="%m/%d/%Y")
bltime = cdata.sub[which(cdata.sub$VISCODE2=='bl'),c("RID",'edate')]
colnames(bltime) = c('RID','vtime')
cdata.sub = merge(cdata.sub,bltime,by = 'RID')
cdata.sub$time2 = with(cdata.sub, (edate - vtime)/365.25)


##Dataset for analysis
adata = cdata.sub[,c('RID',"VISCODE2",'time','time2',"ADNI_MEM","ADNI_EF")]
adata = merge(adata,edata[,c("RID","AGE","PTEDUCAT","PTGENDER",'DX_bl')],by = 'RID')
adata$Age75 = adata$AGE-75
adata$Educ16 = adata$PTEDUCAT-16
adata$Sex = as.factor(adata$PTGENDER)
adata$Sex = factor(adata$Sex,levels = c("0","1"),labels = c('Male',"Female"))
adata$DX = as.factor(adata$DX_bl)
adata$DX = factor(adata$DX, levels = c("1","2","3"),labels = c("NL","EMCI","LMCI"))
### Combine Cognitive data with AB, Tau, Atrophy ###
eVars = c("RID","ABETA_bl","TAU_bl",colnames(edata)[grep('Atrophy',x = colnames(edata))])
fdata = merge(adata,edata[,eVars],by = 'RID')

View(fdata)
### DAN'S BLOM TRANSFORM
#   New code - full sample for adni_mem and adni_ef, bl for csf & brain

blomvarin <- c("ADNI_MEM","ADNI_EF")
blomvarout <- c("memst","execst")
med <- recodeBlom("fdata",varlist_orig=blomvarin,varlist_tr=blomvarout)

blomvarin <- c(colnames(fdata)[grep('Atrophy',colnames(fdata))],
               "ABETA_bl","TAU_bl")
blomvarout <- c("amygst","entorst","hippst","insst",
                "ltemst","mtemst","phippst","pcingst","slfst","splenst","thalst",
                "abetast","taust")


med4 <- med[med$VISCODE2=='bl',c("RID",blomvarin)]
med4 <- recodeBlom("med4",varlist_orig=blomvarin,varlist_tr=blomvarout)

med <- merge(med,med4[,c("RID",blomvarout)],by="RID")

saveRDS(fdata,file = '~/Documents/Davis 2016-2017/ADC_Evan/R_SERVER/fdata1.RDA')
saveRDS(med,file = '~/Documents/Davis 2016-2017/ADC_Evan/R_SERVER/fdata.RDA')
