# step 1: install the correct Rtools for R complier 
# step 2: modify the "pack_root" to your directory, or send your Sys.info()[['nodename']] to us 
# step 3: install from source by the following 

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
library(party)
library(mobForest)

data("BostonHousing", package = "mlbench")
BostonHousing <- BostonHousing[,c("rad", "tax", "crim", "medv", "lstat")]
rfout <- mobForestAnalysis(as.formula(medv ~ lstat),
                           c("rad", "tax", "crim"), mobForest.controls =
                             mobForest_control(ntree = 3, mtry = 2, replace = TRUE,
                                               alpha = 0.05, bonferroni = TRUE,
                                               minsplit = 25),
                           data = BostonHousing, processors = 3, model = linearModel)
pred = getPredictedValues(object=rfout, newdata=T, 
                          newTestData=BostonHousing[1, ])

#
# in the pred result, you may check statistics for each tree, each single node, and some average stat accross the trees
# documentation is not supplied, please let me know if there is any question. 

