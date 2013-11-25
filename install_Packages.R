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
install.packages( paste( pack_root, '/PARTY/party_1.0-10/party', sep = ''),verbose=T,
                  repos = NULL, type = 'source')
install.packages( paste( pack_root, '/PARTY/mobForest_1.2/mobForest', sep = ''),verbose=T,
                  repos = NULL, type = 'source')
install.packages( paste( pack_root, '/PARTY/modeltools_0.2-21/modeltools', sep = ''),verbose=T,
                  repos = NULL, type = 'source')
library(ggplot2)
library(modeltools)
library(party)
library(mobForest)

# load test data
data("BostonHousing", package = "mlbench")

# test for mob 
BostonHousing$lstat <- log(BostonHousing$lstat)
BostonHousing$rm <- BostonHousing$rm^2
BostonHousing$chas <- factor(BostonHousing$chas, levels = 0:1, labels = c("no", "yes"))
BostonHousing$rad <- factor(BostonHousing$rad, ordered = TRUE)
objfun_method.list = c( 'min', 'sum')
correl.obj.list = list(correl = correl.obj
#                        , correl.obj.x_y = correl.obj.x_y
#                        , deviance = deviance
                       )
fmBH.list = list()
all.obj = data.frame()
for ( objfun_method in objfun_method.list){
  for ( obj in names(correl.obj.list )){
    
    fmBH <- mob(medv ~ lstat + rm | zn + indus + chas + nox + age + dis + rad + tax + crim + b + ptratio,
                control = mob_control(minsplit = 5, 
                                      objfun = correl.obj.list[[obj]],
                                      objfun_method = objfun_method,
                                      verbose = F),
                data = BostonHousing, model = linearModel)
    type = paste(objfun_method, obj, sep = '_' )
    fmBH.list[[type]] = fmBH
    temp.condition = Condition.SplittingNode(fmBH@tree, objfun=correl.obj.list[[obj]])
    temp.condition$type = type
    all.obj = plyr::rbind.fill(temp.condition, all.obj)
  }
    
}
undebug(predict)
BostonHousing$node_mob = predict(object= fmBH, newdata=BostonHousing, type = 'node')
sort(unique(BostonHousing$node_mob))
sort(party:::terminal_nodeIDs(fmBH@tree))
debug(predict)
ggplot(data=all.obj, 
       mapping=aes(x = nodeID, 
                y = -objFunValue,
                col =type )) + 
  geom_line()

# test for mobForest
data("BostonHousing", package = "mlbench")
rfout <-
  mobForestAnalysis(formula = as.formula(medv ~ lstat),
                    partitionVariables = c("rad", "tax", "crim"),
                    mobForest.controls = 
                      mobForest_control(ntree = 3, mtry = 2, replace = TRUE,
                                        alpha = 0.05, bonferroni = TRUE,
                                        objfun = deviance,
                                        
                                        # here the correl.obj is embeded in the modeltools, but you may supply your own objfun
                                       # also the correl.obj is calcualting the negative correlation ( for minimizing)
                                        minsplit = 25
#                                         ,objfun_method= 'min'
                                        ),
                    data = BostonHousing, 
                    processors = 3, 
                    model = linearModel)
test.nodeID = 9
tree_id = '2'
test.mob = rfout@mf.trees[[tree_id]]
plot(test.mob )
# test.weights = test.mob@tree$left$weights
# # str(test.mob@tree$left$model$ModelEnv@get("designMatrix"))
# lm.test1 = lm(medv ~  lstat, data = BostonHousing,weights= test.weights)
# lm.test2 = lm(medv ~ lstat, data = BostonHousing[which(test.weights==1),])
# test$model$weights == lm.test1$weights
# str(test$model$ModelEnv@get("designMatrix"))
# z <- lm.wfit(object@get("designMatrix"),
#              object@get("responseMatrix"), weights, ...)


Condition.SplittingNode(test.mob@tree, objfun=correl.obj)
# test.model = nodes(test.mob, test.nodeID)[[1]]
# length(which(test.model$weights>0))
BostonHousing$node = 
  predict(object= test.mob,
          newdata = BostonHousing,#[493,c('lstat',"rad", "tax", "crim")], 
          type = 'node')
head(BostonHousing)
sort(unique(BostonHousing$node))
# debug(predict)
# lm.model = lm( medv ~ lstat,
#                data =BostonHousing[which(BostonHousing$node == test.nodeID),] 
#                )

pred = getPredictedValues(object=rfout, newdata=T, 
                          newTestData=BostonHousing)
rfout@mf.trees.prop[[tree_id]]
pred$node[1,as.integer(tree_id)]
BostonHousing[1, ]
BostonHousing$node[1]
sort(unique(pred$node[,as.integer(tree_id)]))
rfout@mf.trees.prop[[1]]
#
# in the pred result, you may check statistics for each tree, each single node, and some average stat accross the trees
# documentation is not supplied, please let me know if there is any question. 

# eg.
# > pred
# response
# [,1]     [,2]     [,3]
# [1,] 24.95467 29.74449 29.03649
# $node
# [,1] [,2] [,3]
# [1,]   12    7    1
# $nodeAdjR2
# [,1]      [,2]     [,3]
# [1,] 0.3608695 0.5866803 0.519687
# $nodeObjFunValue
# [,1]       [,2]       [,3]
# [1,] -0.6809822 -0.7729968 -0.8159765
# $predMat
# PredMean PredStdev  Residual AvgNodeAdjR2 AvgNodeObjFunValue
# [1,] 27.91188  2.585372 -3.911881    0.4890789         -0.7566518
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
