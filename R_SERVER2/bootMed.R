

################################
### BOOTSTRAP Function #########
################################

bootMed = function(B,atVar = atrophy,cogVar = cog,dxVar = dx,dBoot = dataB){
resultsB = 
  lapply(atVar,function(x){
   lapply(cogVar,function(y){
    lapply(1:length(dxVar),function(z){
        rData = NULL
        count=1
        while(count<=B){
        bList = data.frame(B.Run = count, Atrophy = x, Cog = y, DX = names(dx[z]))
        rList = medmod(cogVar = y,atVar = x,dxgroup = unlist(dx[z]),df = data.frame(dBoot[count]))
        fList = cbind(bList,rList)
        rData = rbind(rData,fList)
        count = count+1
      }
    rData
})})})

return(list(ResultsBoot = resultsB))
}



