#########################################
#### EVAN MEDIATION ANALYSIS RESULTS ####
#########################################
setwd("~/Documents/Davis 2016-2017/ADC_Evan/Neuroimaging-Simulation-Study/R_RESULTS2")
library(DT)
library(data.table)
library(Rmisc)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
`%ni%` <- Negate(`%in%`) 


temp = list.files(pattern="*.RDA")
filenames <- list.files(pattern="*.RDA")
ldf <- lapply(filenames, readRDS)
res <- lapply(ldf, function(r){
  rbind(data.frame(r$ResultsBoot[[1]][[1]][[1]]),
        data.frame(r$ResultsBoot[[1]][[1]][[2]]),
        data.frame(r$ResultsBoot[[1]][[1]][[3]]),
        data.frame(r$ResultsBoot[[1]][[1]][[4]]),
        data.frame(r$ResultsBoot[[1]][[2]][[1]]),
        data.frame(r$ResultsBoot[[1]][[2]][[2]]),
        data.frame(r$ResultsBoot[[1]][[2]][[3]]),
        data.frame(r$ResultsBoot[[1]][[2]][[4]]))
})

for(i in 1:length(temp)){assign(paste('dt',i,sep = ''),data.table(data.frame(res[i])))}
temp1 = lapply(paste0("dt", 1:length(temp)),get)

for (i in 1:length(temp1)) {
  
  # read data:
  sample <- data.frame(temp1[i])
  sample$C1 = with(sample, ab2-abc)
  sample$C2 = with(sample, ab2-a2c)
  sample$C3 = with(sample, abc-a2c)
  sample$dPrime = with(sample, d - (abc+a2c+ab2))
  sample$ab = with(sample, a * b)
  assign(paste('dt',i,sep = ''),sample)
  
}
temp1 = lapply(paste0("dt", 1:length(temp)),get)


medSum = function(x){list(mean = round(mean(x),6),
                          se = round(sd(x),6),
                          LCI = round(quantile(x,.025),6),
                          UCI = round(quantile(x,.975),6))}
medPaths =c('a','b','c','d','dP','a2','b2','ab','bc','abc','ab2','a2c',
            'C1','C2','C3','dPrime')

mP = rep(unlist(lapply(medPaths,function(x){rep(x,4)})),8)
tP = rep(c('mean','se','LCI','UCI'),8*length(medPaths))

dataSum1 = function(dt,med){
  dt[,unlist(lapply(.SD, medSum)), .SDcols = med,by = .(Atrophy,Cog,DX)]
}

dataSum2 = function(dt,med){
  dtNew = dataSum1(dt,med)
  dtNew[,MedPath:=mP][,Type:=tP]
  dcast(dtNew,Atrophy+Cog + DX + MedPath ~ Type, value.var = "V1")
}

for(i in 1:length(temp1)){assign(paste('sum',i,sep = ''),dataSum2(data.table(temp1[[i]]),medPaths))}

totalSum <- do.call(rbind, lapply(paste0("sum", 1:length(temp1)), get))
mPnames = data.frame(med = medPaths,
                     Paths = c('Abeta to Tau','Tau to Atrophy','Atrophy to Cog',
                               'Abeta to Cog','dP = Abeta to Cog',
                               'Abeta to Atrophy','Tau to Cog','Abeta to Tau to Atrophy',
                               'Tau to Atrophy to Cog',
                               'B2 = Abeta to Tau to Atrophy to Cog',
                               'B1 = Abeta to Tau to Cog',
                               'B3 = Abeta to Atrophy to Cog',
                               'C1 = ab2-abc',
                               'C2 = ab2-a2c',
                               'C3 = abc-a2c',
                               'dPrime = d-(abc+a2c+ab2)'))
totalSum = merge(totalSum, mPnames,by.x = 'MedPath',by.y = 'med')


totalSum$Sig = ifelse(sign(totalSum$LCI)==sign(totalSum$UCI),"1","0")
totalSum$Sig = factor(totalSum$Sig)
totalSum$MedPath = factor(totalSum$MedPath)
totalDF = data.frame(totalSum)
totalDF = totalDF[,c(2:4,1,9,10,5:8)]
totalSumDT = datatable(totalDF,filter = 'top',
                       options = list(pageLength = 50,autoWidth = TRUE
                                      #columnDefs = list(list(visible=FALSE, targets=c(5:6)))
                                      ),
                       caption = 'Table 1: This table summarizes results of 1000 bootstrapped samples by Atrophy, Cognitive Measure, Diagnosistic Category and Mediation Path (medPath).
                       mean = average of 1000 estimates, se = standard error of 1000 estimates, Sig = 1: the value 0 was not in the 95% Conf band, 0: the value zero was in the 95%CI band'
                       
)

totalSumDT

write.csv(totalDF,file = 'resultsSum.csv')
tr = read.csv('resultsSum.csv')


### OUTPUT TABLES
library(plyr)

### Normals

### EMCI
### LMCI
### TOTAL
options(scipen=999)
sumFunc = function(beta,mem){
  sdf = subset(totalDF,MedPath==beta&Cog==mem)
  ddply(sdf, .(DX,Atrophy),summarise,
        Path = MedPath,
        cog = Cog,
        Beta = mean,
        CI = paste('[',LCI,', ',UCI,']',sep = ''),
        sig = Sig
        )
}
mpaths = c('ab2','a2c','abc','ab')
dpaths = c('a2','d','a','dP','dPrime')
cogpaths = c('ADNI_MEM','ADNI_EF')
Diags = unique(sample$DX)

Fig34 = lapply(cogpaths,
                     function(i){
                       lapply(mpaths,function(x){sumFunc(x,i)})})

for(d in 1:length(Diags)){
for(i in 1:length(cogpaths)){
for(j in 1:length(mpaths)){
    sample = data.frame(Fig34[[i]][[j]])
    subdf = subset(sample,DX==Diags[d])
    colnames(subdf)[which(colnames(subdf)=="Beta")] = if(unique(subdf$Path)=='ab2'){'B1'}else{if(unique(subdf$Path)=='abc'){'B2'}else{'B3'}}
    subdf$Atrophy = gsub("^.*?_","",subdf$Atrophy)
    assign(paste(Diags[d],cogpaths[i],mpaths[j],sep = '_'),subdf)
  }}}



Fig2 =lapply(cogpaths,
                      function(i){
                        lapply(dpaths,function(x){sumFunc(x,i)})})

for(d in 1:length(Diags)){
for(i in 1:length(cogpaths)){
  for(j in 1:length(dpaths)){
    sample = data.frame(Fig2[[i]][[j]])
    subdf = subset(sample,DX==Diags[d])
    subdf$Atrophy = gsub("^.*?_","",subdf$Atrophy)
    assign(paste(Diags[d],cogpaths[i],dpaths[j],sep = '_'),subdf)
  }}}

### AB correlations

fig2CN = cbind(" " = c('ROI',rep(" ",11),'Cognition',rep(" ",2),'Other'," "), 
               "Variable" = c(" ",NL_ADNI_EF_a2$Atrophy,
                 c(" ",'EF',"MEM", " "),
                 'Tau'),
               rbind(c(" ", " ", " "), data.frame(NL_ADNI_EF_a2[,c("Beta","CI",'sig')]), # controlling for tau
             c(" ", " ", " "),
             data.frame(NL_ADNI_EF_dP[1,c("Beta","CI",'sig')]),
             data.frame(NL_ADNI_MEM_dP[1,c("Beta","CI",'sig')]),
             c(" ", " ", " "),
             data.frame(NL_ADNI_EF_a[1,c("Beta","CI",'sig')])))
fig2EMCI = cbind(" " = c('ROI',rep(" ",11),'Cognition',rep(" ",2),'Other'," "), 
               "Variable" = c(" ",EMCI_ADNI_EF_a2$Atrophy,
                              c(" ",'EF',"MEM", " "),
                              'Tau'),
               rbind(c(" ", " ", " "), data.frame(EMCI_ADNI_EF_a2[,c("Beta","CI",'sig')]), # controlling for tau
                     c(" ", " ", " "),
                     data.frame(EMCI_ADNI_EF_dP[1,c("Beta","CI",'sig')]),
                     data.frame(EMCI_ADNI_MEM_dP[1,c("Beta","CI",'sig')]),
                     c(" ", " ", " "),
                     data.frame(EMCI_ADNI_EF_a[1,c("Beta","CI",'sig')])))

fig2LMCI = cbind(" " = c('ROI',rep(" ",11),'Cognition',rep(" ",2),'Other'," "), 
                 "Variable" = c(" ",LMCI_ADNI_EF_a2$Atrophy,
                                c(" ",'EF',"MEM", " "),
                                'Tau'),
                 rbind(c(" ", " ", " "), data.frame(LMCI_ADNI_EF_a2[,c("Beta","CI",'sig')]), # controlling for tau
                       c(" ", " ", " "),
                       data.frame(LMCI_ADNI_EF_dP[1,c("Beta","CI",'sig')]),
                       data.frame(LMCI_ADNI_MEM_dP[1,c("Beta","CI",'sig')]),
                       c(" ", " ", " "),
                       data.frame(LMCI_ADNI_EF_a[1,c("Beta","CI",'sig')])))

fig2TOTAL = cbind(" " = c('ROI',rep(" ",11),'Cognition',rep(" ",2),'Other'," "), 
                 "Variable" = c(" ",TOTAL_ADNI_EF_a2$Atrophy,
                                c(" ",'EF',"MEM", " "),
                                'Tau'),
                 rbind(c(" ", " ", " "), data.frame(TOTAL_ADNI_EF_a2[,c("Beta","CI",'sig')]), # controlling for tau
                       c(" ", " ", " "),
                       data.frame(TOTAL_ADNI_EF_dP[1,c("Beta","CI",'sig')]),
                       data.frame(TOTAL_ADNI_MEM_dP[1,c("Beta","CI",'sig')]),
                       c(" ", " ", " "),
                       data.frame(TOTAL_ADNI_EF_a[1,c("Beta","CI",'sig')])))



library(xtable)
figure2Func = function(fig2df,caption){
sigCol1 = which(fig2df$sig==1)
fig2df$Variable = as.character(fig2df$Variable)
fig2df$Variable[which(fig2df$Variable=='SLF_PT')]='SLF.PT'
print(xtable(fig2df, 
       caption = caption,
       digits = 4), hline.after=c(12,15),booktabs = TRUE,include.rownames = FALSE, caption.placement = 'top',
      add.to.row=list(
        pos=list(as.list(sigCol1-1))[[1]],
        command=rep("\\rowcolor{green!10}",
                    length(seq(from=1,to=length(sigCol1),by=1)))),
      sanitize.text.function=identity,table.placement = 'H'
)
}


## B2 and B3 
## Mediations and correlations
fig3EFCN = cbind(data.frame(NL_ADNI_EF_abc[,c("B2","CI",'sig','Atrophy')]),
                   data.frame(NL_ADNI_EF_a2c[,c("B3","CI",'sig')]))
fig3MEMCN = cbind(data.frame(NL_ADNI_MEM_abc[,c("B2","CI",'sig','Atrophy')]),
                    data.frame(NL_ADNI_MEM_a2c[,c("B3","CI",'sig')]))


fig3EFEMCI = cbind(data.frame(EMCI_ADNI_EF_abc[,c("B2","CI",'sig','Atrophy')]),
               data.frame(EMCI_ADNI_EF_a2c[,c("B3","CI",'sig')]))
fig3MEMEMCI = cbind(data.frame(EMCI_ADNI_MEM_abc[,c("B2","CI",'sig','Atrophy')]),
               data.frame(EMCI_ADNI_MEM_a2c[,c("B3","CI",'sig')]))

fig3EFLMCI = cbind(data.frame(LMCI_ADNI_EF_abc[,c("B2","CI",'sig','Atrophy')]),
                   data.frame(LMCI_ADNI_EF_a2c[,c("B3","CI",'sig')]))
fig3MEMLMCI = cbind(data.frame(LMCI_ADNI_MEM_abc[,c("B2","CI",'sig','Atrophy')]),
                    data.frame(LMCI_ADNI_MEM_a2c[,c("B3","CI",'sig')]))

fig3EFTOTAL = cbind(data.frame(TOTAL_ADNI_EF_abc[,c("B2","CI",'sig','Atrophy')]),
                   data.frame(TOTAL_ADNI_EF_a2c[,c("B3","CI",'sig')]))
fig3MEMTOTAL = cbind(data.frame(TOTAL_ADNI_MEM_abc[,c("B2","CI",'sig','Atrophy')]),
                    data.frame(TOTAL_ADNI_MEM_a2c[,c("B3","CI",'sig')]))

figure3Func = function(fig3df,caption){
  fig3df$Atrophy = as.character(fig3df$Atrophy)
  fig3df$Atrophy[which(fig3df$Atrophy=='SLF_PT')]='SLF.PT'
  colnames(fig3df)[which(colnames(fig3df)=='Atrophy')]='ROI'
  col1 = which(fig3df[,3]==1&fig3df[,7]==0)
  col2 = which(fig3df[,3]==0&fig3df[,7]==1)
  col3 = which(fig3df[,3]==1&fig3df[,7]==1)
  col4 = which(fig3df[,3]==0&fig3df[,7]==0)


    fig3df$B3 = ifelse(fig3df[,3]==1&fig3df[,7]==0, paste0("\\cellcolor{white}{", fig3df$B3, "}"), fig3df$B3)
  fig3df[,6] = ifelse(fig3df[,3]==1&fig3df[,7]==0, paste0("\\cellcolor{white}{", fig3df[,6], "}"), fig3df[,6])
  
  fig3df$B2 = ifelse(fig3df[,3]==0&fig3df[,7]==1, paste0("\\cellcolor{white}{", fig3df$B2, "}"), fig3df$B2)
  fig3df[,2] = ifelse(fig3df[,3]==0&fig3df[,7]==1, paste0("\\cellcolor{white}{", fig3df[,2], "}"), fig3df[,2])
  
  colnames(fig3df)[1:3] = paste0("\\cellcolor{red!30}{", colnames(fig3df)[1:3], "}")
  colnames(fig3df)[5:7] = paste0("\\cellcolor{blue!30}{", colnames(fig3df)[5:7], "}")
  print(xtable(fig3df[,-c(3,7)], digits = 4,caption = caption),caption.placement = 'top',
        include.rownames = FALSE,
        booktabs = TRUE,
        add.to.row=list(
          pos=list(as.list(c(col1,col2,col3)-1))[[1]],
          command=c(rep("\\rowcolor{red!30}",ifelse(length(col1)==0,0,
                      length(seq(from=1,to=length(col1),by=1)))),
                    rep("\\rowcolor{blue!30}",ifelse(length(col2)==0,0,
                        length(seq(from=1,to=length(col2),by=1)))),
                    rep("\\rowcolor{green!30}",ifelse(length(col3)==0,0,
                        length(seq(from=1,to=length(col3),by=1)))))),
        sanitize.text.function = function(x){x},table.placement = 'H'
  )
}


## B1 tables
B1EFCN = data.frame(NL_ADNI_EF_ab2[,c('Atrophy',"B1","CI",'sig')])
B1EFEMCI = data.frame(EMCI_ADNI_EF_ab2[,c('Atrophy',"B1","CI",'sig')])
B1EFLMCI = data.frame(LMCI_ADNI_EF_ab2[,c('Atrophy',"B1","CI",'sig')])
B1EFTOTAL = data.frame(TOTAL_ADNI_EF_ab2[,c('Atrophy',"B1","CI",'sig')])

B1MEMCN = data.frame(NL_ADNI_MEM_ab2[,c('Atrophy',"B1","CI",'sig')])
B1MEMEMCI = data.frame(EMCI_ADNI_MEM_ab2[,c('Atrophy',"B1","CI",'sig')])
B1MEMLMCI = data.frame(LMCI_ADNI_MEM_ab2[,c('Atrophy',"B1","CI",'sig')])
B1MEMTOTAL = data.frame(TOTAL_ADNI_MEM_ab2[,c('Atrophy',"B1","CI",'sig')])

figure4Func = function(fig3df,caption){
  fig3df$Atrophy = as.character(fig3df$Atrophy)
  fig3df$Atrophy[which(fig3df$Atrophy=='SLF_PT')]='SLF.PT'
  colnames(fig3df)[which(colnames(fig3df)=='Atrophy')]='ROI'
  col1 = which(fig3df$sig==1)
  print(xtable(fig3df, digits = 4,caption = caption),caption.placement = 'top',
        include.rownames = FALSE,
        booktabs = TRUE,
        add.to.row=list(
          pos=list(as.list(c(col1)-1))[[1]],
          command=c(rep("\\rowcolor{orange!30}",
                        ifelse(length(col1)==0,0,
                               length(seq(from=1,to=length(col1),by=1)))))),
        sanitize.text.function = function(x){x},table.placement = 'H'
  )
}

## ab tables
abEFCN = data.frame(NL_ADNI_EF_ab[,c('Atrophy',"B3","CI",'sig')])
abEFEMCI = data.frame(EMCI_ADNI_EF_ab[,c('Atrophy',"B3","CI",'sig')])
abEFLMCI = data.frame(LMCI_ADNI_EF_ab[,c('Atrophy',"B3","CI",'sig')])
abEFTOTAL = data.frame(TOTAL_ADNI_EF_ab[,c('Atrophy',"B3","CI",'sig')])

abMEMCN = data.frame(NL_ADNI_MEM_ab[,c('Atrophy',"B3","CI",'sig')])
abMEMEMCI = data.frame(EMCI_ADNI_MEM_ab[,c('Atrophy',"B3","CI",'sig')])
abMEMLMCI = data.frame(LMCI_ADNI_MEM_ab[,c('Atrophy',"B3","CI",'sig')])
abMEMTOTAL = data.frame(TOTAL_ADNI_MEM_ab[,c('Atrophy',"B3","CI",'sig')])


### Demographic Tables

demData = readRDS("~/Documents/Davis 2016-2017/ADC_Evan/Neuroimaging-Simulation-Study/R_SERVER2/fdata1.RDA")

# The data set is sorted by subject RID.
first <- which(!duplicated(demData$RID))
last <- c(first[-1] -1, nrow(demData))

demFirst = demData[first,]
demLast = demData[last,]
demLast$MEMbl = demFirst$ADNI_MEM
demLast$EFbl = demFirst$ADNI_EF
demLast$MEM2yr = (demLast$ADNI_MEM - demLast$MEMbl)/as.numeric(demLast$time2)*2
demLast$EF2yr = (demLast$ADNI_EF - demLast$EFbl)/as.numeric(demLast$time2)*2

library(plyr)

ddply(demLast,.(DX),summarise,
      N = length(DX),
      Age = mean(AGE,na.rm =TRUE),
      "Age sd" = sd(AGE,na.rm = TRUE),
      Educ = mean(PTEDUCAT,na.rm =TRUE),
      "Educ sd" = sd(PTEDUCAT,na.rm = TRUE),
      "CSF AB" = mean(ABETA_bl,na.rm = TRUE),
      "CSF AB sd" = sd(ABETA_bl,na.rm = TRUE),
      "CSF Tau" = mean(TAU_bl,na.rm = TRUE),
      "CSF Tau sd" = sd(TAU_bl,na.rm = TRUE),
      "MEM Baseline" = mean(MEMbl,na.rm = TRUE),
      "MEM Baseline sd" = sd(MEMbl,na.rm = TRUE), 
      "EF Baseline" = mean(EFbl,na.rm = TRUE),
      "EF Baseline sd" = sd(EFbl,na.rm = TRUE),
      "MEM 2yr Change" = mean(MEM2yr,na.rm = TRUE),
      "MEM 2yr Change sd" = sd(MEM2yr,na.rm = TRUE), 
      "EF 2yr Change" = mean(EF2yr,na.rm = TRUE),
      "EF 2yr Change sd" = sd(EF2yr,na.rm = TRUE),
      "MTL 2yr Change" = mean(Atrophy_MediatTemporal,na.rm = TRUE),
      "MTL 2yr Change sd" = sd(Atrophy_MediatTemporal,na.rm = TRUE),
      "LTR 2yr Change" = mean(Atrophy_LateralTemporal,na.rm = TRUE),
      "LTR 2yr Change sd" = sd(Atrophy_LateralTemporal,na.rm = TRUE))
      

## FINAL CI graphics
library(dplyr)
fig_CN  = fig2CN[,-1]
fig_CN = fig_CN[which(fig_CN$Beta!=" "),]

fig_tidy = function(df){
df_tidy = select(tbl_df(df),Variable:sig)
df_tidy = filter(df_tidy,Beta >" ")

Figure_df  = df_tidy %>% separate(CI,c("lowerCI","upperCI"), ", ") %>% mutate(Beta = as.numeric(Beta))
Figure_df$upperCI = as.numeric(unlist(lapply(Figure_df$upperCI,function(x){gsub("\\]","",x)})))
Figure_df$lowerCI = as.numeric(unlist(lapply(Figure_df$lowerCI,function(x){gsub("\\[","",x)})))
Figure_df$sig = ifelse(Figure_df$sig==0,'Not Significant','Significant')
Figure_df = Figure_df[order(Figure_df$Beta),] %>% filter(Variable!='MediatTemporal')
Figure_df$Variable = factor(Figure_df$Variable, levels = Figure_df$Variable)
df1 =  filter(Figure_df,Variable%ni%c('EF','MEM','Tau'))

list(df1 = df1)
}
B2_EMCIdf_mem = with(fig3MEMEMCI,data.frame(Variable = Atrophy,Beta = B2,CI = CI, sig = sig))
B3_EMCIdf_mem = with(fig3MEMEMCI,data.frame(Variable = Atrophy,Beta = B3,CI = fig3MEMEMCI[,ncol(fig3MEMEMCI)-1], sig = fig3MEMEMCI[,ncol(fig3MEMEMCI)]))

B2_LMCIdf_mem = with(fig3MEMLMCI,data.frame(Variable = Atrophy,Beta = B2,CI = CI, sig = sig))
B3_LMCIdf_mem = with(fig3MEMLMCI,data.frame(Variable = Atrophy,Beta = B3,CI = fig3MEMLMCI[,ncol(fig3MEMEMCI)-1], sig = fig3MEMLMCI[,ncol(fig3MEMLMCI)]))

B2_LMCIdf_ef = with(fig3EFLMCI,data.frame(Variable = Atrophy,Beta = B2,CI = CI, sig = sig))
B3_LMCIdf_ef = with(fig3EFLMCI,data.frame(Variable = Atrophy,Beta = B3,CI = fig3EFLMCI[,ncol(fig3EFEMCI)-1], sig = fig3EFLMCI[,ncol(fig3EFLMCI)]))


CNfig = fig_tidy(fig2CN)

B2_EMCIfig_mem = fig_tidy(B2_EMCIdf_mem)
B3_EMCIfig_mem = fig_tidy(B3_EMCIdf_mem)
B2_EMCIfig_mem$df1$Variable = factor(B2_EMCIfig_mem$df1$Variable,level = levels(B3_EMCIfig_mem$df1$Variable))

EMCI = rbind(B2_EMCIfig_mem$df1,B3_EMCIfig_mem$df1)
EMCI$grp = c(rep('b2',10),rep('b3',10))


B2_LMCIfig_mem = fig_tidy(B2_LMCIdf_mem)
B3_LMCIfig_mem = fig_tidy(B3_LMCIdf_mem)
B2_LMCIfig_mem$df1$Variable = factor(B2_LMCIfig_mem$df1$Variable,level = levels(B3_LMCIfig_mem$df1$Variable))


B2_LMCIfig_ef = fig_tidy(B2_LMCIdf_ef)
B3_LMCIfig_ef = fig_tidy(B3_LMCIdf_ef)
B2_LMCIfig_ef$df1$Variable = factor(B2_LMCIfig_ef$df1$Variable,level = levels(B3_LMCIfig_ef$df1$Variable))


LMCI_mem = rbind(B2_LMCIfig_mem$df1,B3_LMCIfig_mem$df1)
LMCI_mem$grp = c(rep('b2',10),rep('b3',10))

LMCI_ef = rbind(B2_LMCIfig_ef$df1,B3_LMCIfig_ef$df1)
LMCI_ef$grp = c(rep('b2',10),rep('b3',10))

library(ReporteRs)
library(ggplot2)
library(gtable)
cn = ylim(-.0001,.00024)
emci = ylim(-0.00022,0.0007)
lmci = ylim(-0.0009,0.002)
lmci2 = ylim(-0.001,0.002)

p = ggplot(EMCI, 
           aes(x = Variable, y = Beta, color = sig)) + 
  geom_point() + 
  facet_grid(~grp,scales = 'free') + 
  coord_flip() + theme_classic()+
  theme(
    legend.position = 'bottom',
    axis.text.x = element_text(size = 24,face = 'bold',colour = 'black'), 
    axis.text.y = element_text(size = 32,face = 'bold',color = 'black'), 
    title = element_text(size = 16, face = 'bold',color = 'black'))+
geom_errorbar(aes(ymin = lowerCI, ymax=upperCI), width=.3, position=position_dodge(0.05),size = 1.1) +
  geom_hline(yintercept = 0,color = 'black',linetype = 2,size = 1.3) + 
scale_colour_manual(values = c('black','red')) +
  labs(color = "")+ylab('')+xlab('') +emci +theme(legend.position = 'bottom')+
ylim(-.0001,.0005)+scale_y_continuous(breaks = c(0.0000,0.0003,0.0005), labels = c('0','3e-4','5e-4'))+
  theme(panel.spacing = unit(18, "lines"))

p
p1 = p+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
g <- ggplotGrob(p)
g1 = ggplotGrob(p1)
axis <- gtable_filter(g, "axis-l")[["grobs"]][[1]][["children"]][["axis"]][,1]
g[["grobs"]][[4]][["children"]][["axis"]] <- NULL

# build plot & add axis to LHS of left facet
panels <- subset(g$layout, name %in% grep('panel',g$layout$name,value = TRUE))
g2 <- gtable_add_grob(g1, grobs=axis, t = unique(panels$t), l=tail(panels$l, -1)-1)

# Locate the tops of the plot panels
# Locate the tops of the plot panels
panels <- grep("panel", g2$layout$name)
top <- unique(g2$layout$t[panels])

g2 = g2[-(top-1), ]
grid.newpage()

png('gEMCI.png',width = 1200, height = 500)

grid.draw(g2)

dev.off()



### LMCI MEM

p = ggplot(LMCI_mem, 
           aes(x = Variable, y = Beta, color = sig)) + 
  geom_point() + 
  facet_grid(~grp,scales = 'free') + 
  coord_flip() + theme_classic()+
  theme(
    legend.position = 'bottom',
    axis.text.x = element_text(size = 24,face = 'bold',colour = 'black'), 
    axis.text.y = element_text(size = 32,face = 'bold',color = 'black'), 
    title = element_text(size = 16, face = 'bold',color = 'black'))+
  geom_errorbar(aes(ymin = lowerCI, ymax=upperCI), width=.3, position=position_dodge(0.05),size = 1.1) +
  geom_hline(yintercept = 0,color = 'black',linetype = 2,size = 1.3) + 
  scale_colour_manual(values = c('black','red')) +
  labs(color = "")+ylab('')+xlab('') +lmci +theme(legend.position = 'bottom')+theme(panel.spacing = unit(18, "lines"))+
  scale_y_continuous(breaks = c(0.0000,0.0005,0.001), labels = c('0','5e-4','1e-3'))

p
p1 = p+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
g <- ggplotGrob(p)
g1 = ggplotGrob(p1)
axis <- gtable_filter(g, "axis-l")[["grobs"]][[1]][["children"]][["axis"]][,1]
g[["grobs"]][[4]][["children"]][["axis"]] <- NULL

# build plot & add axis to LHS of left facet
panels <- subset(g$layout, name %in% grep('panel',g$layout$name,value = TRUE))
g2 <- gtable_add_grob(g1, grobs=axis, t = unique(panels$t), l=tail(panels$l, -1)-1)

# Locate the tops of the plot panels
# Locate the tops of the plot panels
panels <- grep("panel", g2$layout$name)
top <- unique(g2$layout$t[panels])

g2 = g2[-(top-1), ]
grid.newpage()

png('gLMCImem.png',width = 1200, height = 500)

grid.draw(g2)

dev.off()




### LMCI EF

p = ggplot(LMCI_ef, 
           aes(x = Variable, y = Beta, color = sig)) + 
  geom_point() + 
  facet_grid(~grp,scales = 'free') + 
  coord_flip() + theme_classic()+
  theme(
    legend.position = 'bottom',
    axis.text.x = element_text(size = 24,face = 'bold',colour = 'black'), 
    axis.text.y = element_text(size = 30,face = 'bold',color = 'black'), 
    title = element_text(size = 16, face = 'bold',color = 'black'))+
  geom_errorbar(aes(ymin = lowerCI, ymax=upperCI), width=.3, position=position_dodge(0.05),size = 1.1) +
  geom_hline(yintercept = 0,color = 'black',linetype = 2,size = 1.3) + 
  scale_colour_manual(values = c('black','red')) +
  labs(color = "")+ylab('')+xlab('') +lmci2 +theme(legend.position = 'bottom')+theme(panel.spacing = unit(18, "lines"))+
  scale_y_continuous(breaks = c(0.0000,0.001), labels = c('0','1e-3'))

p
p1 = p+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
g <- ggplotGrob(p)
g1 = ggplotGrob(p1)
axis <- gtable_filter(g, "axis-l")[["grobs"]][[1]][["children"]][["axis"]][,1]
g[["grobs"]][[4]][["children"]][["axis"]] <- NULL

# build plot & add axis to LHS of left facet
panels <- subset(g$layout, name %in% grep('panel',g$layout$name,value = TRUE))
g2 <- gtable_add_grob(g1, grobs=axis, t = unique(panels$t), l=tail(panels$l, -1)-1)

# Locate the tops of the plot panels
# Locate the tops of the plot panels
panels <- grep("panel", g2$layout$name)
top <- unique(g2$layout$t[panels])

g2 = g2[-(top-1), ]
grid.newpage()

png('gLMCIef.png',width = 1200, height = 500)

grid.draw(g2)

dev.off()








plot1CI = function(df,name,ylimit){

  themePlots = theme(legend.position = 'bottom',axis.text.x = element_text(size = 12,face = 'bold',colour = 'black'), 
                     axis.text.y = element_text(size = 12,face = 'bold',color = 'black'), 
                     title = element_text(size = 16, face = 'bold',color = 'black')) + theme_classic()
  
  plot1 = ggplot(data = df$df1, aes(x = Variable, y = Beta,color = sig))+geom_point() + 
  geom_errorbar(aes(ymin = lowerCI, ymax=upperCI), width=.3, position=position_dodge(0.05),size = 1) + 
  coord_flip() + 
  geom_hline(yintercept = 0,color = 'black',linetype = 2,size = 1.3) + 
  scale_colour_manual(values = c('black','red')) +
  labs(color = "")+ylab('')+xlab('') +ylimit + themePlots +theme(legend.position = 'bottom')

doc = docx()

doc = addPlot(doc = doc, fun = print, x = plot1)
writeDoc(doc, file = name, height = 2, width = 2)
list(plot1)
}

cn = ylim(-.0001,.00024)
emci = ylim(-0.00022,0.0007)
lmci = ylim(-0.0009,0.002)
lmci2 = ylim(-0.001,0.002)

df.cn = list(CNfig)
df.cn.list = list('cn.docx')
df.list = list(B2_EMCIfig_mem, B3_EMCIfig_mem)
df.listLMCI = list(B2_LMCIfig_ef,B3_LMCIfig_ef)
df.listLMCI2 = list(B2_LMCIfig_mem,B3_LMCIfig_mem)
file.list = list('EMCI1.docx','EMCI2.docx')
file.listLMCI = list('l1.docx','l2.docx')
file.listLMCI2 = list('l3.docx','l4.docx')

plotCN =  mapply(function(x,y){plot1CI(x,y,cn)}, x = df.cn, y = df.cn.list)
plotR = mapply(function(x,y){plot1CI(x,y,emci)}, x = df.list, y = file.list)
plotR1 = mapply(function(x,y){plot1CI(x,y,lmci)}, x = df.listLMCI, y = file.listLMCI)
plotR2 = mapply(function(x,y){plot1CI(x,y,lmci2)}, x = df.listLMCI2, y = file.listLMCI2)


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(plotR[[1]])

png('plotCN.png',width = 800,height = 800)
plot = grid.arrange(arrangeGrob(plotCN[[1]] + theme(legend.position="none"),
                                    nrow=1),
                        mylegend,nrow = 2,heights=c(25, .5))
dev.off()



png('plotEMCI.png',width = 800,height = 800)
plotEMCI = grid.arrange(arrangeGrob(plotR[[1]] + theme(legend.position="none") +
                                    theme(axis.title.y=element_blank(),
                                          axis.text.y=element_blank(),
                                          axis.ticks.y=element_blank())+ylim(-.0001,.0005)+scale_y_continuous(breaks = c(0.0000,0.00025,0.00050)),
                         plotR[[2]] + theme(legend.position="none")+scale_y_continuous(breaks = c(0.0000,0.00025,0.00050)),
                         nrow=1),
             mylegend,nrow = 2,heights=c(1, 1))
dev.off()

png('plotLMCI_EF.png',width = 800,height = 800)
plotEMCI = grid.arrange(arrangeGrob(plotR1[[1]] + theme(legend.position="none")+ggtitle(expression(beta[2]))+theme(title  = element_text(size = 16,face = 'bold'))+
                                      theme(axis.title.x=element_blank(),
                                            axis.text.x=element_blank(),
                                            axis.ticks.x=element_blank()),
                                    plotR1[[2]] + theme(legend.position="none")+ggtitle(expression(beta[3]))+theme(title  = element_text(size = 16,face = 'bold')),
                                    nrow=1),
                        mylegend,nrow = 2,heights=c(1, 1))
dev.off()
png('plotLMCI_MEM.png',width = 800,height = 800)
plotEMCI = grid.arrange(arrangeGrob(plotR2[[1]] + theme(legend.position="none")+ggtitle(expression(beta[2]))+theme(title  = element_text(size = 16,face = 'bold')),
                                    plotR2[[2]] + theme(legend.position="none")+ggtitle(expression(beta[3]))+theme(title  = element_text(size = 16,face = 'bold')),
                                    nrow=1),
                        mylegend,nrow = 2,heights=c(1, 1))
dev.off()