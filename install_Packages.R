# step 1: install the correct Rtools for R complier 
# step 2: modify the "pack_root" to your directory, or send your Sys.info()[['nodename']] to us 
# step 3: uncomment the install part and install from source by the following codes

if(Sys.info()[['nodename']]=='US-WASH-23C9KQ1'){
  pack_root<<-"C:/MuniGit"
}else if(Sys.info()[['nodename']]=='AU-EQUITYPC'){
  pack_root<<-"C:/Source"
}else if(Sys.info()[['nodename']] %in%  c('ERMUNIQC65',"ERMUNIQC66",'ERMUNIQC85',"ERMUNIQC86")){
  pack_root <<- "D:"
}
setwd(pack_root)
set.seed(290875)

lib.path = .libPaths()[1]  # path where to install the packages e.g "C:/Program Files/R/library" 

#c:\program Files\R\R-3.0.1\bin\x64\R CMD Install -l "C:\Program Files\R\library" c:\..\party
install.packages('sandwich', lib=lib.path)
install.packages('strucchange', lib=lib.path)
install.packages('coin', lib=lib.path)
install.packages( paste( pack_root, '/PARTY/modeltools_0.2-21/modeltools', sep = ''),verbose=T,
                  repos = NULL, lib=lib.path, type = 'source')
install.packages( paste( pack_root, '/PARTY/party_1.0-10/party', sep = ''),verbose=T,
                  repos = NULL, lib=lib.path, type = 'source')
install.packages( paste( pack_root, '/PARTY/mobForest_1.2/mobForest', sep = ''),verbose=T,
                  repos = NULL, lib=lib.path, type = 'source')

library(mobForest)
# load test data
# test for mobForest
data("BostonHousing", package = "mlbench")
ntree = 2
objfun = correl.obj
# here the correl.obj is embeded in the modeltools, but you may supply your own objfun
# also the correl.obj is calcualting the negative correlation ( for minimizing)
objfun_method= 'min' # 'min' or 'sum' which defines the method to make the splitting 

caltime = system.time({
  rfout <-
    mobForestAnalysis(formula = as.formula(medv ~ lstat),
                      partitionVariables = c("rad", "tax", "crim"),
                      mobForest.controls = 
                        mobForest_control(ntree = ntree, mtry = 1, replace = T,
                                          alpha = 0.05, bonferroni = TRUE,
                                          objfun = objfun,
                                          
                                          minsplit = 10
                                          ,objfun_method= objfun_method
                        ),
                      data = BostonHousing, 
                      processors = min ( ntree, 6), 
                      model = linearModel)
})

print(caltime/ntree)
# check out the output
rfout@mf.trees.prop # show the objective functions, AdjR2 for each node each tree
rfout@mf.trees # save the trees in the forest 

tree_id = '1'
test.mob = rfout@mf.trees[[tree_id]]
plot(test.mob )

# outsample prediction  
pred = getPredictedValues(object=rfout, newdata=T, 
                          newTestData=BostonHousing[1,])
pred$predMat # average values accross all trees 
pred$response # list of responses from all trees 
pred$node # predicted leaf node from all trees 
pred$predConditions # paths to each tree for each variable 

# prediction output eg. 
# > pred$predMat
# PredMean PredStdev  Residual AvgNodeAdjR2 AvgNodeObjFunValue
# [1,] 29.20678 0.3749711 -5.206776     0.575237         -0.8429686
# > pred$response
# [,1]     [,2]
# [1,] 28.94163 29.47192
# > pred$predConditions
# [[1]]
# tax_LessThan tax_MoreThan nodeID objFunValue     AdjR2 tree_id
# 1           NA          243      5  -0.8410226 0.5833734       1
# 2           NA           NA      1  -0.8449145 0.5671005       2
# 
# > pred$predMat # average values accross all trees 
# PredMean PredStdev  Residual AvgNodeAdjR2 AvgNodeObjFunValue
# [1,] 29.20678 0.3749711 -5.206776     0.575237         -0.8429686
# > pred$response # list of responses from all trees 
# [,1]     [,2]
# [1,] 28.94163 29.47192
# > pred$node # predicted leaf node from all trees 
# [,1] [,2]
# [1,]    5    1
# > pred$predConditions # paths to each tree for each variable
# [[1]]
# tax_LessThan tax_MoreThan nodeID objFunValue     AdjR2 tree_id
# 1           NA          243      5  -0.8410226 0.5833734       1
# 2           NA           NA      1  -0.8449145 0.5671005       2
