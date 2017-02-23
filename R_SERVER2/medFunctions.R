################################
### Mediation Model Functions ##
################################

################################
## medBoot(B,data,dBoot) #######
################################

# B = numeric, # of bootstrap samples
# data = data.frame, final data set, fdata
# dBoot = bootstrap samples (by RID)
# returns a data.frame with bootstrap sample data, input for medmod()

medBoot = function(B,data,dBoot){
  dfBoot = dBoot[,c(2,(B+2))]
  df = merge(data,dfBoot,by.x = 'RID',by.y = colnames(dfBoot)[2],all.x = FALSE)
  df =df[order(df$RID,df$ID,df$VISCODE2),]
  return(df)
}

######################################
## medmod(cogVar,atVar,dxgroup,df) ###
######################################

# cogVar: character, cognitive variable, options = 'ADNI_MEM' or 'ADNI_EF'
# atVar: character, atrophy ROI, options = see atrophy list.
# dxgroup: character, options = 'NL','EMCI','LMCI'
# df: data.frame, using medBoot function to generate bootstrapped data set.
# returns the vector resB = c(a,a2,b,b2,c,d,a*b,b*c,a*b*c,a*b2,a2*c)
medmod = function(cogVar,atVar,dxgroup,df){
  df$cog = df[,cogVar]
  df$atrophy = df[,atVar]
  subdf = subset(df,DX %in% dxgroup)
  
  ######################################################################################  
  ######################################################################################
  
  ## a: ABeta ---> Tau
  
  fmA = lm(taust~abetast+Age75+Educ16+Sex,data = subset(subdf,VISCODE2=='bl'))
  
  ## d: ABeta ---> Cog
  
  lmeD = if(class(try(
             lme(cog ~ time2 * (abetast + Age75 + Sex + Educ16),
             random = ~(1+time2)|ID,
             data = subdf,
             na.action = na.omit),
             silent = TRUE))=='try-error'){
    ### if random slope won't converge, just run random intercept.
             lme(cog ~ time2 * (abetast + Age75 + Sex + Educ16),
             random = ~(1)|ID,
             data = subdf,
             na.action = na.omit)}else{
    ### if convergences, take full model           
             lme(cog ~ time2 * (abetast +
                         Age75 + Sex + Educ16),
        random = ~(1+time2)|ID,
        data = subdf,
        na.action = na.omit
    )
    
  }
  
  ######################################################################################  
  ######################################################################################
  ## Run D, B2, and C simultaneously estimating Cog with random effects. 
  
  ### Atrophy Models
  
  ## a2: ABeta ---> Atrophy
  ## b:  Tau ---> Atrophy
  
  fmA2 = lm(atrophy~abetast+taust+Age75+Educ16+Sex,data = subset(subdf,VISCODE2=='bl'))
  
  ## b2: Tau ---> Cog
  ## c: Atrophy ---> Cog
  
  lmeB2 = if(class(try(
              lme(cog ~ time2 * (abetast + taust + atrophy +Age75 + Sex + Educ16),
              random = ~(1+time2)|ID,
              data = subdf,
              na.action = na.omit),
              silent = TRUE))=='try-error'){
    ### if random slope won't converge, just run random intercept.
              lme(cog ~ time2 * (abetast + taust + atrophy +Age75 + Sex + Educ16),
              random = ~(1)|ID,
              data = subdf,
              na.action = na.omit)}else{
    ### if convergences, take full model           
                lme(cog ~ time2 * (abetast + taust + atrophy +Age75 + Sex + Educ16),
                random = ~(1+time2)|ID,
                data = subdf,
                na.action = na.omit
                )  
              }
  
  ### Extract Coefficients 
  
  a = as.numeric(coef(fmA)['abetast'])
  d = as.numeric(fixed.effects(lmeD)['time2:abetast'])
  a2 = as.numeric(coef(fmA2)['abetast'])
  b = as.numeric(coef(fmA2)['taust'])
  b2 = as.numeric(fixed.effects(lmeB2)['time2:taust'])
  c = as.numeric(fixed.effects(lmeB2)['time2:atrophy'] )
  dP = as.numeric(fixed.effects(lmeB2)['time2:abetast'] )
  
  ######################################################################################  
  ######################################################################################
  resB = data.frame(a = a,a2 = a2,b = b,b2 = b2,c = c,d = d,dP = dP,
                    ab = a*b,bc = b*c,abc = a*b*c,ab2 = a*b2,a2c = a2*c)
  resB
}
