# step 1: install the correct Rtools for R complier 
# step 2: modify the "pack_root" to your directory, or send your Sys.info()[['nodename']] to us 
# step 3: uncomment the install part and install from source by the following codes

if(Sys.info()[['nodename']]=='US-WASH-23C9KQ1'){
  pack_root<<-"C:/MuniGit"
}else if(Sys.info()[['nodename']]=='AU-EQUITYPC'){
  pack_root<<-"C:/Source"
}else if(Sys.info()[['nodename']] %in%  c('ERMUNIQC65',"ERMUNIQC66",'ERMUNIQC85',"ERMUNIQC86")){
  pack_root <<- "D:/models"
}
setwd(pack_root)
set.seed(290875)

# 
# install.packages( paste( pack_root, '/PARTY/party_1.0-10/party', sep = ''),verbose=T,
#                   repos = NULL, type = 'source'
#                   )
# install.packages( paste( pack_root, '/PARTY/mobForest_1.2/mobForest', sep = ''),verbose=T,
#                   repos = NULL, type = 'source')
# install.packages( paste( pack_root, '/PARTY/modeltools_0.2-21/modeltools', sep = ''),verbose=T,
#                   repos = NULL, type = 'source')
library(modeltools)
library(party)
library(mobForest)

data("BostonHousing", package = "mlbench")
BostonHousing <- BostonHousing[,c("rad", "tax", "crim", "medv", "lstat")]
rfout <- mobForestAnalysis(as.formula(medv ~ lstat),
                           c("rad", "tax", "crim"), mobForest.controls =
                             mobForest_control(ntree = 3, mtry = 2, replace = TRUE,
                                               alpha = 0.05, bonferroni = TRUE,
                                               
                                               objfun = correl.obj,
                                               # here the correl.obj is embeded in the modeltools, but you may supply your own objfun
                                               # also the correl.obj is calcualting the negative correlation ( for minimizing)
                                               minsplit = 25),
                           data = BostonHousing, processors = 3, model = linearModel)

pred = getPredictedValues(object=rfout, newdata=T, 
                          newTestData=BostonHousing[1, ])
#
# in the pred result, you may check statistics for each tree, each single node, and some average stat accross the trees
# documentation is not supplied, please let me know if there is any question. 

# # eg.
# > pred
# response
# [,1]     [,2]     [,3]
# [1,] 24.95467 29.74449 29.03649
# 
# $node
# [,1] [,2] [,3]
# [1,]   12    7    1
# 
# $nodeAdjR2
# [,1]      [,2]     [,3]
# [1,] 0.3608695 0.5866803 0.519687
# 
# $nodeObjFunValue
# [,1]       [,2]       [,3]
# [1,] -0.6809822 -0.7729968 -0.8159765
# 
# $predMat
# PredMean PredStdev  Residual AvgNodeAdjR2 AvgNodeObjFunValue
# [1,] 27.91188  2.585372 -3.911881    0.4890789         -0.7566518
# 
# $predConditions
# $predConditions[[1]]
# tax_LessThan tax_MoreThan nodeID objFunValue     AdjR2 crim_LessThan crim_MoreThan
# 1          305          277     12  -0.6809822 0.3608695            NA            NA
# 2           NA          242      7  -0.7729968 0.5866803            NA            NA
# 3           NA           NA      1  -0.8159765 0.5196870            NA            NA
# tree_id rad_LessThan rad_MoreThan
# 1       1           NA           NA
# 2       2            2           NA
# 3       3           NA           NA
