# Version 4 was made because parseFormula() had a bug that needed to be rectified. 
# Version 5 for adding the generalized linear model functionality

get_condition = function (condition, x ){
  # x should be a psplit class
  # condition is a data frame 
  left_condition = condition 
  right_condition = condition 
  if( class(x)=="orderedSplit"){
    max.name = paste(x$variableName, 'LessThan', sep = '_')
    min.name = paste(x$variableName, 'MoreThan', sep = '_')
    if ( max.name%in% names( condition)){
      left_condition[[max.name]] = min(left_condition[[max.name]] , x$splitpoint, na.rm=T)
      right_condition[[min.name]] = max(right_condition[[min.name]] , x$splitpoint, na.rm=T)
    }else{
      left_condition[[max.name]] = x$splitpoint
      right_condition[[min.name]] = x$splitpoint
      left_condition[[min.name]] = NA
      right_condition[[max.name]] = NA
    }
  }
  
  if( class(x) == "nominalSplit"){
    levels <- attr(x$splitpoint, "levels")
    tab <- x$table
    left_condition[[x$variableName]] <- paste(c(levels[as.logical(x$splitpoint) & (tab > 0)]), collapse=";" )
    right_condition[[ x$variableName]] <-paste(c(levels[!as.logical(x$splitpoint) & (tab > 0)]), collapse=";" )
  }
  return ( list('left'=left_condition,
                'right'= right_condition
  ))
}

Condition.SplittingNode <- function(x, condition =NULL, objfun = NULL, ...) {
  # get condition path and objective function value for splitting node
  if (x$terminal){
    # this should be a very rare case, because only the node 1 == terminal node
    objFunValue = 0 
    if ( !is.null( objfun)){
      objFunValue = objfun(x$model) 
    }
    return (data.frame(nodeID = x$nodeID, 
                       objFunValue = objFunValue, 
                       AdjR2 = summary(x$model)$adj.r.squared
                       ))
  }
  
  condition.result = data.frame()
  if(is.null ( condition)){
    condition = list()
  }
  condition = get_condition(condition,x$psplit )
  left_condition = condition[['left']]
  right_condition = condition[['right']]
  
  condition.result.left = data.frame()
  condition.result.right = data.frame()
  if ( x$left$terminal){
    objFunValue = 0 
    if ( !is.null( objfun)){
      objFunValue = objfun(x$left$model) 
    }
    condition.result.left = data.frame(t(sapply(left_condition,c)))
    condition.result.left$nodeID = x$left$nodeID
    condition.result.left$objFunValue = objFunValue
    condition.result.left$AdjR2 = summary(x$left$model)$adj.r.squared
    
  }else{
    # pass to left 
    #     left_condition = rbind( condition, left_condition)
    
    condition.result = 
      plyr:::rbind.fill( condition.result, 
                  Condition.SplittingNode (x$left, condition = left_condition, objfun=objfun ) 
      )
  }
  
  if( x$right$terminal){
    
    objFunValue = 0 
    if ( !is.null( objfun)){
      objFunValue = objfun(x$right$model)
    }
    
    condition.result.right = data.frame(t(sapply(right_condition,c)))
    condition.result.right$nodeID = x$right$nodeID
    condition.result.right$objFunValue = objFunValue
    condition.result.right$AdjR2 = summary(x$right$model)$adj.r.squared
    
    #     condition.result.right = data.frame(nodeID = x$right$nodeID, 
    #                                         path = right_condition, 
    #                                         objFunValue = objFunValue)
    condition.result = plyr:::rbind.fill(condition.result,condition.result.right, condition.result.left  )
  }else{
    # pass condition to right 
    
    condition.result = 
      plyr:::rbind.fill( condition.result,
                  condition.result.left, 
                  Condition.SplittingNode (x$right, condition = right_condition, objfun=objfun ) 
      )
  }
  return ( condition.result)
}

mob_RF_Tree <- function (mainModel, partitionVars, mtry, weights, data = list(), na.action = na.omit, model = glinearModel, control = mob_control(), ...)
{	
  # sampling the variables to split
	formula = formula(paste(mainModel, paste(sample(partitionVars)[1:mtry], collapse=" + "), sep=" | "))	# add	
	if (inherits(formula, "formula")) {
        mobpp <- function(formula, data, model) {
            ff <- attr(modeltools:::ParseFormula(formula), "formula")
			#ff <- attr(ParseFormula(formula), "formula")
            ff$input[[3]] <- ff$input[[2]]
            ff$input[[2]] <- ff$response[[2]]
            dpp(model, as.formula(ff$input), other = list(part = as.formula(ff$blocks)), data = data, na.action = na.action)
        }
        formula <- mobpp(formula, data, model)
    }
    if (missing(weights))
        weights <- rep(1, dimension(formula, "part")[1])
    fm <- fit(model, formula, ...)
    where <- integer(length(weights))
    mob_fit <- function(obj, mf, weights, control) {
        obj <- reweight(obj, weights)
        if (inherits(obj, "try-error")) {
            node <- list(nodeID = NULL, weights = weights, criterion = list(statistic = 0, 
                criterion = 0, maxcriterion = 0), terminal = TRUE, 
                psplit = NULL, ssplits = NULL, prediction = 0, 
                left = NULL, right = NULL, sumweights = as.double(sum(weights)), 
                model = obj)
            class(node) <- "TerminalNodeModel"
            node$nodeID <- as.integer(nodeid)
            where[weights > 0] <<- as.integer(nodeid)
            nodeid <<- nodeid + 1
            return(node)
        }
        # XXX here is a problem: if the mob_fit_setupnode returns a False or Inf, then, 
        thisnode <- party:::mob_fit_setupnode(obj, mf, weights, control)
        # XXX fixed a bug here to handle special case of terminal node
        if (identical(FALSE, thisnode)) {
          node <- list(nodeID = NULL, weights = weights,
                       criterion = list(statistic = 0, 
                                        criterion = 0, 
                                        maxcriterion = 0),
                       terminal = TRUE, psplit = NULL, ssplits = NULL,
                       prediction = 0, left = NULL, right = NULL, 
                       sumweights = as.double(sum(weights)), 
                       model = obj)
          class(node) <- "TerminalModelNode"  
          node$nodeID <- as.integer(nodeid)
          where[weights > 0] <<- as.integer(nodeid)
          nodeid <<- nodeid + 1
          return(node)
        }
        
        
        
        thisnode$nodeID <- as.integer(nodeid)
        where[weights > 0] <<- as.integer(nodeid)
        nodeid <<- nodeid + 1
        thisnode$model <- obj
        if (!thisnode$terminal) {
            childweights <- party:::mob_fit_childweights(thisnode, mf,weights)
            if (any(sapply(childweights, sum) == 0)) {
                thisnode$terminal <- TRUE
                class(thisnode) <- "TerminalModelNode"
                return(thisnode)
            }
			mf = formula(paste(mainModel, paste(sample(partitionVars)[1:mtry], collapse=" + "), sep=" | "))	# add			
			mf <- mobpp(mf, data, model) # add
			thisnode$left <- mob_fit(obj, mf, weights = childweights$left, control)			
            thisnode$right <- mob_fit(obj, mf, weights = childweights$right, control)
        }
        return(thisnode)
    }
    nodeid <- 1
    tr <- mob_fit(fm, formula, weights = weights, control = control)
    y <- formula@get("response")
    yy <- new("VariableFrame", nrow(y), ncol(y))
    yy@variables <- formula@get("response")
    rval <- new("mob", tree = tr, responses = yy, data = formula, 
        where = where)
    return(rval)
}

computeR2 <- function(response, predictions)
{
	mean.predictions = apply(predictions, 1, mean, na.rm=T)
	sse = sum((response[,1] - mean.predictions)**2, na.rm=T)
	ssto = sum((response[,1] - colMeans(response, na.rm=T))**2, na.rm=T)
	R2 = 1 - (sse/ssto)	
	if(R2 < 0) R2 = 0
	return(R2)
}

computeAcc <- function(response, predictions, prob.cutoff)
{
	#pred.class = rep(0, length(response[,1]))
	#prop.1 = apply(predictions, 1, function(x) mean(x == 1))
	#pred.class[which(prop.1 > 0.5)] = 1
	#acc = mean(pred.class == response[,1])
	#return(acc)
	pred.class = rep(0, length(response[,1]))
	levels(response[,1]) <- list("0"=levels(response[,1])[1], "1"=levels(response[,1])[2])
	meanPred = apply(predictions, 1, mean, na.rm=T)
	pred.class[which(meanPred > prob.cutoff)] = 1
	acc = mean(pred.class == response[,1])
	return(acc)
}

computeMSE <- function(response, predictions)
{
  mean.predictions = apply(predictions, 1, mean, na.rm=T)
  sse = sum((response[,1] - mean.predictions)**2, na.rm=T)
  n = length(which(!is.na((response[,1] - colMeans(response, na.rm=T)))))
  #ssto = sum((response[,1] - colMeans(response, na.rm=T))**2, na.rm=T)
  #R2 = 1 - (sse/ssto)	
  #if(R2 < 0) R2 = 0
  return(sse/n)
}


treePredictions <- function(j, data, tree)
{
  # prediction for one single tree, given the newdata and tree 
  # XXX this is a good example to extract the terminal node for one single data 
  # and the 'predict_response' method is defined in the modeltools package
	tr = tree
	dat = data
	while(tr$terminal == FALSE)
	{
		newvar = tr$psplit$variableName
		if(class(tr$psplit) == "nominalSplit")
		{
			leftSplitVars = levels(tr$psplit$splitpoint)[as.logical(tr$psplit$splitpoint)]
			if (!is.na(match(as.character(dat[j,newvar]), leftSplitVars)))
			{	
				tr = tr$left
			} else {
				tr = tr$right
			}				
		} else {	
      # XXX should change all the ordered variable split into numerical
			if (as.numeric(dat[j,newvar]) <= as.numeric(tr$psplit$splitpoint)) 
			{	
				tr = tr$left
			} else {
				tr = tr$right
			}			
		}			
	}
# 	return(as.numeric(tr$model$predict_response(dat[j,])))	
	return(list ( response = as.numeric(tr$model$predict_response(dat[j,])),
                node = tr$nodeID
#                 , AdjR2 = summary(tr$model)$adj.r.squared
                ))	
} 

# bootstrap over a list of tree numbers and generate tree for each of them
bootstrap <- function(i, data, mainModel, partitionVars, mtry, newTestData, mob.controls, 
                      fraction, replace, model, family, prob.cutoff)
{    
  if(is.null(prob.cutoff)) prob.cutoff = 0.5
  # sample the data sets 
  data.subinds = sample(nrow(data), replace = replace)[1:round(fraction*nrow(data))]	
  data.sub = data[data.subinds,]
  fmBH = NULL
  if(model@name == "linear regression model")
    fmBH <- mob_RF_Tree(mainModel = mainModel, 
                        partitionVars = partitionVars, 
                        mtry = mtry, 
                        control = mob.controls,
                        data = data.sub, 
                        model = model)
  if(model@name == "generalized linear regression model") # only difference is the family in difference GLM classes 
    fmBH <- mob_RF_Tree(mainModel = mainModel, 
                        partitionVars = partitionVars,
                        mtry = mtry,
                        control = mob.controls,
                        data = data.sub, 
                        model = model, 
                        family = family)
  
  oob.inds = setdiff(1:nrow(data), data.subinds)
  oob.sub = data[oob.inds,]
  
#   pred = sapply(1:nrow(data), treePredictions, data = data, tree = fmBH@tree)
  pred = sapply(1:nrow(data), 
                function ( j){
                  temp.pred = treePredictions(j, data = data, tree = fmBH@tree)
                    return ( temp.pred[['response']])
                }
                )
  
  # get the response variable from the model formula, instead of quoting the column names
  obs.outcome <- ModelEnvFormula(as.formula(mainModel), data = data)@get("response")
  
  
  ret = NULL
  if ((is.null(fmBH@tree$model$family)) ||
        (fmBH@tree$model$family$link == "log" && 
           fmBH@tree$model$family$family == "poisson")
  )
  {
    gen.rsq = computeR2(obs.outcome, matrix(pred, ncol=1))
    mse.gen = sum((obs.outcome - pred)**2, na.rm=T)/nrow(data)
    oob.rsq = computeR2(matrix(obs.outcome[oob.inds,], ncol=1), matrix(pred[oob.inds], ncol=1))
    #mse.oob = sum((obs.outcome[oob.inds,] - pred[oob.inds])**2, na.rm=T)/length(oob.inds)
    mse.oob = computeMSE(matrix(obs.outcome[oob.inds,],ncol=1), matrix(pred[oob.inds], ncol=1))
    oob.rsq.perm = rep(0, length(partitionVars))
    oob.mse.perm = rep(0, length(partitionVars))
    varImp = rep(0, length(partitionVars))
    for(p in 1:length(partitionVars))
    {		#XXX permutation over the partitioning variables to get the variable importance in each single tree
      oob.perm = oob.sub
      oob.perm[,partitionVars[p]] = sample(oob.sub[,partitionVars[p]])
#       oob.pred.perm = sapply(1:nrow(oob.perm), 
#                              treePredictions, 
#                              data = oob.perm, tree = fmBH@tree)			
      oob.pred.perm = sapply(1:nrow(oob.perm), 
                    function ( j){
                      temp.pred = treePredictions(j, data = oob.perm, tree = fmBH@tree)
                      return ( temp.pred[['response']])
                    }
      )
      
      oob.rsq.perm[p] = computeR2(matrix(obs.outcome[oob.inds,],ncol=1), matrix(oob.pred.perm, ncol=1))
      oob.mse.perm[p] = computeMSE(matrix(obs.outcome[oob.inds,],ncol=1), matrix(oob.pred.perm, ncol=1))
      #oob.rsq.perm[p] = 1 - (sse.perm/ssto.oob)		
    }
    
    pred.new = c()
    newdat.rsq = c()
    if(nrow(newTestData) > 0) 
    {
#       pred.new = sapply(1:nrow(newTestData), 
#                         treePredictions, 
#                         data = newTestData, tree = fmBH@tree)
      pred.new = sapply(1:nrow(newTestData), 
                             function ( j){
                               temp.pred = treePredictions(j, data = newTestData, tree = fmBH@tree)
                               return ( temp.pred[['response']])
                             }
      )
      
      obs.newdat <- ModelEnvFormula(as.formula(paste(mainModel, partitionVars, sep=" | ")),
                                    data = newTestData)@get("response")
      newdat.rsq = computeR2(obs.newdat, matrix(pred.new, ncol=1))
    }
    ret <- list(oob.inds, oob.rsq, pred, 
                (oob.mse.perm - mse.oob),
                mse.oob, gen.rsq, mse.gen, pred.new, 
                newdat.rsq, fmBH@tree)	
    names(ret) <- c("oob.inds", "oob.R2", "pred",
                    "rawVarImp",
                    "mse.oob", "gen.R2", "mse.gen", "pred.new", "newdat.R2", 
                    'mf.single.tree')	
  }
  
  if (model@name == "generalized linear regression model")
  {		
    if (fmBH@tree$model$family$link == "logit")
    {	
      levels(obs.outcome[,1]) <- list("0"=levels(obs.outcome[,1])[1], "1"=levels(obs.outcome[,1])[2])
      pred.class = rep(0, length(pred))
      pred.class[which(pred > prob.cutoff)] = 1		
      gen.acc = length(which(pred.class == obs.outcome[,1]))
      oob.acc = length(which(pred.class[oob.inds] == obs.outcome[oob.inds,1]))
      oob.acc.perm = rep(0, length(partitionVars))
      varImp = rep(0, length(partitionVars))
      for(p in 1:length(partitionVars))
      {		
        oob.perm = oob.sub
        oob.perm[,partitionVars[p]] = sample(oob.sub[,partitionVars[p]])
#         oob.pred.perm = sapply(1:nrow(oob.perm), treePredictions, data = oob.perm, tree = fmBH@tree)
        oob.pred.perm = sapply(1:nrow(oob.perm), 
                               function ( j){
                                 temp.pred = treePredictions(j, data = oob.perm, tree = fmBH@tree)
                                 return ( temp.pred[['response']])
                               }
        )
        
        
        pred.class.perm = rep(0, length(oob.pred.perm))
        pred.class.perm[which(oob.pred.perm > prob.cutoff)] = 1
        #oob.acc.perm[p] = computeR2(obs.outcome[oob.inds,], matrix(oob.pred.perm, ncol=1))
        oob.acc.perm[p] = length(which(pred.class.perm == obs.outcome[oob.inds,1]))
      }
      
      pred.new = c()
      newdat.acc = c()
      if(nrow(newTestData) > 0)
      {
#         pred.new = sapply(1:nrow(newTestData), 
#                           treePredictions, 
#                           data = newTestData,
#                           tree = fmBH@tree)
        pred.new = sapply(1:nrow(newTestData), 
                          function ( j){
                            temp.pred = treePredictions(j, data = newTestData, tree = fmBH@tree)
                            return ( temp.pred[['response']])
                          }
        )
        
        
        predNew.class = rep(0, length(pred.new))
        predNew.class[which(pred > prob.cutoff)] = 1
        obs.newdat <- ModelEnvFormula(as.formula(paste(mainModel, partitionVars, sep=" | ")), data = newTestData)@get("response")
        levels(obs.newdat[,1]) <- list("0"=levels(obs.newdat[,1])[1], "1"=levels(obs.newdat[,1])[2])
        newdat.acc = length(which(predNew.class == obs.newdat[,1]))/nrow(newTestData)
      }
      
      ret <- list(oob.inds, 
                  (oob.acc/length(oob.inds)),
                  pred, 
                  (oob.acc - oob.acc.perm)/length(oob.inds),
                  (gen.acc/nrow(data)),
                  pred.new, 
                  newdat.acc,
                  fmBH@tree
      )
      
      names(ret) <- c("oob.inds", "oob.acc", "pred", "rawVarImp",
                      "gen.acc", "pred.new", "newdat.acc", 
                      'mf.single.tree')
      
      # XXX here is the place what the tree returns 
    }
  }
  return(ret)
}


getmobForestObject.LM <- function(object, mainModel, 
                                  partitionVars, data, newTestData, ntree, fam,objfun)
{		
	B = ntree
	pp.out <- object
	varImpMatrix = matrix(0, nrow=length(partitionVars), ncol=B)
	rownames(varImpMatrix) = partitionVars	
	oob.R2 = c()
	oob.mse = c()
	general.R2 = c()
	general.mse = c()
	general.predictions = matrix(NA, ncol=B, nrow= length(pp.out[[1]]$pred))
	oob.predictions = matrix(NA, ncol=B, nrow= length(pp.out[[1]]$pred))
	mf.trees = list()
	mf.trees.prop = list()
	for(i in 1:B)
	{
		oob.R2[i] = pp.out[[i]]$oob.R2
		oob.mse[i] = pp.out[[i]]$mse.oob
		general.R2[i] = pp.out[[i]]$gen.R2
		general.mse[i] = pp.out[[i]]$mse.gen
		general.predictions[,i] = pp.out[[i]]$pred
		oob.predictions[pp.out[[i]][[1]],i] = pp.out[[i]]$pred[pp.out[[i]]$oob.inds]
		varImpMatrix[,i] = pp.out[[i]]$rawVarImp	
		mf.trees[[i]] = pp.out[[i]]$mf.single.tree
		temp.prop = Condition.SplittingNode(pp.out[[i]]$mf.single.tree, 
		                                    condition = NULL, 
		                                    objfun=objfun ) # save the tree's property here: objective function, and path, etc
		row.names(temp.prop) = temp.prop$nodeID
		mf.trees.prop[[i]] =  temp.prop
		
	}
	names(mf.trees) = c(1:B)
  names(mf.trees.prop) = c(1:B)
  
	obs.outcome <- ModelEnvFormula(as.formula(mainModel), data = data)@get("response")
	oobRes = obs.outcome[,1] - apply(oob.predictions, 1, mean, na.rm=T)
	oobPred <- prediction_output(apply(oob.predictions, 1, mean, na.rm=T), apply(oob.predictions, 1, sd, na.rm=T), oobRes, oob.R2, oob.mse, computeR2(obs.outcome, oob.predictions), "OOB")	
	genRes = obs.outcome[,1] - apply(general.predictions, 1, mean, na.rm=T)	
	generalPred <- prediction_output(apply(general.predictions, 1, mean, na.rm=T), apply(general.predictions, 1, sd, na.rm=T), genRes, general.R2, general.mse, computeR2(obs.outcome, general.predictions), "General")		
	
	newdataPred = prediction_output()
	newdata = newTestData
	newdata.obs = data.frame(matrix(0,0,0))
	if(nrow(newdata) > 0) 
	{
		pred.newdata = lapply(1:B, function(x) pp.out[[x]]$pred.new)
		newdata.predictions = matrix(unlist(pred.newdata), ncol=B)
		newdata.obs <- ModelEnvFormula(as.formula(paste(mainModel, paste(partitionVars, collapse=" + "), sep=" | ")), data = newdata)@get("response")		
		newdata.R2 = unlist(lapply(1:B, function(x) pp.out[[x]]$newdat.R2))
		newdatRes = newdata.obs[,1] - apply(newdata.predictions, 1, mean, na.rm=T)
		#newdataPred <- prediction_output(apply(newdata.predictions, 1, function(x) sum(oob.R2*x, na.rm=T)/sum(oob.R2, na.rm=T)), newdatRes, newdata.R2, numeric(), computeR2(newdata.obs, newdata.predictions), "Newdata")		
		newdataPred <- prediction_output(predMean = apply(newdata.predictions, 1, mean, na.rm=T), predSd = apply(newdata.predictions, 1, sd, na.rm=T), residual = newdatRes, R2 = newdata.R2,  overallR2 = computeR2(newdata.obs, newdata.predictions), predType = "Newdata")  			
	}
	varImpObj <- varimp_output(varImpMatrix)
	mfout <- mobForest_output(oobPred, generalPred, newdataPred, varImpObj, 
                            paste(mainModel, paste(partitionVars, collapse=" + "), sep=" | "),
                            fam = "", train.response = obs.outcome, 
                            new.response = newdata.obs, 
	                          mf.trees = mf.trees, 
	                          mf.trees.prop = mf.trees.prop)
	return(mfout)
}

getmobForestObject.GLM <- function(object, mainModel, partitionVars, data, newTestData, ntree, fam, prob.cutoff,objfun)
{		
	if(is.null(prob.cutoff)) prob.cutoff = 0.5
      
  B = ntree
	pp.out <- object
	varImpMatrix = matrix(0, nrow=length(partitionVars), ncol=B)
	rownames(varImpMatrix) = partitionVars	
	oob.acc = c()
	general.acc = c()
	general.predictions = matrix(NA, ncol=B, nrow= length(pp.out[[1]]$pred))
	oob.predictions = matrix(NA, ncol=B, nrow= length(pp.out[[1]]$pred))
	mf.trees = list()
	mf.trees.prop = list()
	for(i in 1:B)
	{
		oob.acc[i] = pp.out[[i]]$oob.acc
		general.acc[i] = pp.out[[i]]$gen.acc
		general.predictions[,i] = pp.out[[i]]$pred
		oob.predictions[pp.out[[i]]$oob.inds,i] = pp.out[[i]]$pred[pp.out[[i]]$oob.inds]
		varImpMatrix[,i] = pp.out[[i]]$rawVarImp	
		mf.trees[[i]] = pp.out[[i]]$mf.single.tree
    temp.prop = Condition.SplittingNode(pp.out[[i]]$mf.single.tree, 
                                        condition = NULL, 
                                        objfun=objfun ) # save the tree's property here: objective function, and path, etc
    row.names(temp.prop) = temp.prop$nodeID
    
		mf.trees.prop[[i]] =  temp.prop
		
    
	}
	names(mf.trees ) = c(1:B)
  names(mf.trees.prop) = c(1:B)
  
	obs.outcome <- ModelEnvFormula(as.formula(mainModel), data = data)@get("response")
	#levels(obs.outcome[,1]) <- list("0"=levels(obs.outcome[,1])[1], "1"=levels(obs.outcome[,1])[2])
	
	#OOBpred.class = rep(0, length(obs.outcome[,1]))
	#oobtemp = as.numeric(oob.predictions)
	#oobtemp[which(oobtemp > 0.5)] = 1
	#oobtemp[which(oobtemp <= 0.5)] = 0
	#prop.1 = apply(matrix(oobtemp, ncol=B), 1, function(x) mean(x == 1, na.rm=T))
	#OOBpred.class[which(prop.1 > 0.5)] = 1	
	
	#Gpred.class = rep(0, length(obs.outcome[,1]))
	#gtemp = as.numeric(general.predictions)
	#gtemp[which(gtemp > 0.5)] = 1
	#gtemp[which(gtemp <= 0.5)] = 0
	#prop.1 = apply(matrix(gtemp, ncol=B), 1, function(x) mean(x == 1))
	#Gpred.class[which(prop.1 > 0.5)] = 1
	
	#oobPred <- prediction_output(predMean = OOBpred.class, predSd = rep(NA, length(OOBpred.class)) , residual = rep(NA, length(OOBpred.class)), R2 = oob.acc, overallR2 = computeAcc(obs.outcome, oob.predictions), predType = "OOB")
	#generalPred <- prediction_output(predMean = Gpred.class, predSd = rep(NA, length(Gpred.class)), residual = rep(NA, length(Gpred.class)), R2 = general.acc,  overallR2 = computeAcc(obs.outcome, general.predictions), predType = "General")
	oobRes = rep(NA, nrow(obs.outcome))
	oobPred <- prediction_output(predMean = apply(oob.predictions, 1, mean, na.rm=T), predSd = apply(oob.predictions, 1, sd, na.rm=T) , residual = oobRes, R2 = oob.acc, overallR2 = computeAcc(obs.outcome, oob.predictions, prob.cutoff), predType = "OOB")	
	genRes = rep(NA, nrow(obs.outcome))
	generalPred <- prediction_output(predMean = apply(general.predictions, 1, mean, na.rm=T), predSd = apply(general.predictions, 1, sd, na.rm=T) , residual = genRes, R2 = general.acc, overallR2 = computeAcc(obs.outcome, general.predictions, prob.cutoff), predType = "General")
		
	newdataPred = prediction_output()
	newdata = newTestData
	newdata.obs = data.frame(matrix(0,0,0))
	if(nrow(newdata) > 0) 
	{
		pred.newdata = lapply(1:B, function(x) pp.out[[x]]$pred.new)
		newdata.predictions = matrix(unlist(pred.newdata), ncol=B)
		newdata.obs <- ModelEnvFormula(as.formula(paste(mainModel, paste(partitionVars, collapse=" + "), sep=" | ")), data = newdata)@get("response")
		#levels(newdata.obs[,1]) <- list("0"=levels(newdata.obs[,1])[1], "1"=levels(newdata.obs[,1])[2])
		newdata.acc = unlist(lapply(1:B, function(x) pp.out[[x]]$newdat.acc))
		
		#NewDatapred.class = rep(0, length(newdata.obs[,1]))
		#prop.1 = apply(newdata.predictions, 1, function(x) mean(x == 1))
		#NewDatapred.class[which(prop.1 > 0.5)] = 1
		
		#newdataPred <- prediction_output(predMean = NewDatapred.class, predSd = rep(NA, length(NewDatapred.class)), residual = (1 - apply(newdata.predictions,1, mean, na.rm=T)), R2 = newdata.acc,  overallR2 = computeAcc(newdata.obs, newdata.predictions), predType = "Newdata")
		newdatRes = rep(NA, nrow(newdata.obs))
		newdataPred <- prediction_output(predMean = apply(newdata.predictions, 1, mean, na.rm=T), predSd = apply(newdata.predictions, 1, sd, na.rm=T), residual = newdatRes, R2 = newdata.acc,  overallR2 = computeAcc(newdata.obs, newdata.predictions,prob.cutoff), predType = "Newdata")		
	}
	varImpObj <- varimp_output(varImpMatrix)
	mfout <- mobForest_output(oobPred, generalPred, newdataPred, varImpObj,
                            paste(mainModel, paste(partitionVars, collapse=" + "), sep=" | "), 
                            fam = fam, train.response = obs.outcome, 
                            new.response = newdata.obs, 
	                          mf.trees = mf.trees, 
	                          mf.trees.prop = mf.trees.prop
                            )
	return(mfout)
}

stringFormula <- function(formula)
{
	fterms = as.character(formula(terms(formula)))
	outc = fterms[2]
	mod.part = fterms[3]
	mod = paste(outc, "~", mod.part)	
	return(mod)
}

mobForestAnalysis <- function(formula, partitionVariables, 
                              data, 
                              mobForest.controls = mobForest_control(),
                              newTestData = as.data.frame(matrix(0,0,0)), 
                              processors = 1, model = linearModel, 
                              family = NULL, prob.cutoff = NULL)
{
	#library(party)
	mod <- stringFormula(formula)
	partitionVars = partitionVariables
	B = mobForest.controls@ntree
	mtry = mobForest.controls@mtry	
	objfun = mobForest.controls@mob.control$objfun
  
	if(mtry == 0) mtry = round(length(partitionVars)/3)	
	fraction = mobForest.controls@fraction
	if (mobForest.controls@replace == TRUE) fraction = 1		
	#library(parallel)
	cl <- makeCluster(getOption("cl.cores", processors))
	clusterEvalQ(cl, library(party))
	clusterExport(cl, c("mob_RF_Tree", "treePredictions", "computeR2", "computeAcc",
                      "computeMSE"))
  # XXX fmBH is the tree object for each tree id, and added new. 
  ## cluster apply 
	pp.out <- clusterApply(cl, 1:B, 
                         bootstrap, # bootstrap function for randomization and for each single Tree
                         data = data, 
                         mainModel = mod, 
                         partitionVars = partitionVars,
                         mtry = mtry, 
                         newTestData = newTestData, 
                         mob.controls = mobForest.controls@mob.control, 
                         fraction = fraction, 
                         replace = mobForest.controls@replace, 
                         model = model, 
                         family = family, 
                         prob.cutoff = prob.cutoff)
	stopCluster(cl)
	obs.outcome <- ModelEnvFormula(as.formula(paste(mod, partitionVars, sep=" | ")), data = data)@get("response")
	mfObj = NULL
	if (model@name == "generalized linear regression model") {
		if (family$family == "binomial") {
			mfObj <- getmobForestObject.GLM(pp.out, mainModel = mod, 
                                      partitionVars = partitionVars, 
                                      data = data, 
                                      newTestData = newTestData, 
                                      ntree = B, fam = "binomial",
                                      prob.cutoff = prob.cutoff,
			                                objfun = objfun)
		}
		if (family$family == "poisson") {
			mfObj <- getmobForestObject.LM(pp.out, 
                                     mainModel = mod, 
                                     partitionVars = partitionVars, 
                                     data = data, 
                                     newTestData = newTestData, 
                                     ntree = B, 
                                     fam = "poisson", 
			                               objfun = objfun)
		}
	}
	if (model@name == "linear regression model") {
		mfObj <- getmobForestObject.LM(pp.out, 
                                   mainModel = mod, 
                                   partitionVars = partitionVars, 
                                   data = data, 
                                   newTestData = newTestData, 
                                   ntree = B, 
                                   fam = "", 
		                               objfun = objfun)
	}
	return(mfObj)	
}
