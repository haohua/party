# print -> print.mob -> print.SplittingNode -> print.nominal / print.ordinal
# $Id: Print.R 282 2006-10-09 14:23:30Z hothorn $

prettysplit <- function(x, inames = NULL, ilevels = NULL) {
    if (length(x) == 4)
        names(x) <- c("variableID", "ordered", "splitpoint", "splitstatistics")
    if (length(x) == 5)
        names(x) <- c("variableID", "ordered", "splitpoint", "splitstatistics",
                      "toleft")
    if (length(x) == 6)
        names(x) <- c("variableID", "ordered", "splitpoint", "splitstatistics",
                      "toleft", "table")
    if (x$ordered) {
        class(x) <- "orderedSplit"
    } else {
        class(x) <- "nominalSplit"
    }
    if (!is.null(ilevels)) {
        if (!is.null(ilevels[x[["variableID"]]]))
            attr(x$splitpoint, "levels") <- ilevels[[x[["variableID"]]]]
    }
    if (!is.null(inames)) x$variableName <- inames[x[["variableID"]]]
    return(x)
}

prettytree <- function(x, inames = NULL, ilevels = NULL) {
    names(x) <- c("nodeID", "weights", "criterion", "terminal",
                  "psplit", "ssplits", "prediction", "left", "right")
    if (is.null(inames) && extends(class(x), "BinaryTree"))
        inames <- names(x@data@get("input"))
    names(x$criterion) <- c("statistic", "criterion", "maxcriterion")
    names(x$criterion$criterion) <- inames
    names(x$criterion$statistic) <- inames

    if (x$terminal) {
        class(x) <- "TerminalNode"
        return(x)
    }

    x$psplit <- prettysplit(x$psplit, inames = inames, ilevels = ilevels)
    if (length(x$ssplit) > 0)
        x$ssplit <- lapply(x$ssplit, prettysplit, inames = inames, 
                           ilevels = ilevels)

    class(x) <- "SplittingNode"
    x$left <- prettytree(x$left, inames = inames, ilevels = ilevels)   
    x$right <- prettytree(x$right, inames = inames, ilevels = ilevels)    
    return(x)
}
 
print.TerminalNode <- function(x, n = 1, ...) {
    cat(paste(paste(rep(" ", n - 1), collapse = ""), x$nodeID, ")* ", 
                    sep = "", collapse = ""),
        "weights =", sum(x$weights), "\n")
}
 
print.SplittingNode <- function(x, n = 1, ...) {
    # print the splitting node 
    cat(paste(paste(rep(" ", n - 1), collapse = ""), x$nodeID, ") ", sep=""))
    print(x$psplit, left = TRUE) # use print.nominalSplit 
    cat(paste("; criterion = ", round(x$criterion$maxcriterion, 3), 
              ", statistic = ", round(max(x$criterion$statistic), 3), "\n", 
              collapse = "", sep = ""))
    print(x$left, n + 2) # after split once, the code is increased by 2 
    cat(paste(paste(rep(" ", n - 1), collapse = ""), x$nodeID, ") ", sep=""))
    print(x$psplit, left = FALSE)
    cat("\n")
    print(x$right, n + 2)
}
print.orderedSplit <- function(x, left = TRUE, ...) {
  if (!is.null(attr(x$splitpoint, "levels"))) {
    sp <- attr(x$splitpoint, "levels")[x$splitpoint]
  } else {
    sp <- x$splitpoint
  }
  if (!is.null(x$toleft)) left <- as.logical(x$toleft) == left
  if (left) {
    cat(x$variableName, "<=", sp)
  } else {
    cat(x$variableName, ">", sp)
  }
}


print.nominalSplit <- function(x, left = TRUE, ...) {
  
  levels <- attr(x$splitpoint, "levels")
  
  ### is > 0 for levels available in this node
  tab <- x$table
  
  if (left) {
    lev <- levels[as.logical(x$splitpoint) & (tab > 0)]
  } else {
    lev <- levels[!as.logical(x$splitpoint) & (tab > 0)]
  }
  
  txt <- paste("{", paste(lev, collapse = ", "), "}", collapse = "", sep = "")
  cat(x$variableName, "==", txt)
}
Conditioon.mob = function(x, ...){
  if ( class( x@tree) == 'TerminalModelNode'){
    Condition.TerminalNode(x@tree)
  }else if ( class( x@tree) == "SplittingNode"){
    Condition.SplittingNode(x@tree)
  }
}
Condition.TerminalNode <- function(x, n = 1, ...) {
  cat(paste(paste(rep(" ", n - 1), collapse = ""), x$nodeID, ")* ", 
            sep = "", collapse = ""),
      "weights =", sum(x$weights), "\n")
}


try(require(plyr))
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
  # get conditions for splitting node 
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
    
  }else{
    # pass to left 
#     left_condition = rbind( condition, left_condition)
    
    condition.result = 
      rbind.fill( condition.result, 
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
#     condition.result.right = data.frame(nodeID = x$right$nodeID, 
#                                         path = right_condition, 
#                                         objFunValue = objFunValue)
    condition.result = rbind.fill(condition.result,condition.result.right, condition.result.left  )
  }else{
    # pass condition to right 
    
    condition.result = 
      rbind.fill( condition.result,
             condition.result.left, 
             Condition.SplittingNode (x$right, condition = right_condition, objfun=objfun ) 
      )
  }
  return ( condition.result)
}


print.BinaryTreePartition <- function(x, ...)
    print(x@tree)

print.BinaryTree <- function(x, ...) {
    cat("\n")
    cat("\t Conditional inference tree with", length(unique(where(x))), 
        "terminal nodes\n\n")
    y <- x@responses
    if (y@ninputs > 1) {
        cat("Responses:", paste(names(y@variables), 
                                collapse = ", "), "\n")
    }  else {
        cat("Response: ", names(y@variables), "\n")
    }
    inames <- names(x@data@get("input"))
    if (length(inames) > 1) {
        cat("Inputs: ", paste(inames, collapse = ", "), "\n")
    } else {
        cat("Input: ", inames, "\n")
    }
    cat("Number of observations: ", x@responses@nobs, "\n\n")
    print(x@tree)
}

print.RandomForest <- function(x, ...) {
    cat("\n")
    cat("\t Random Forest using Conditional Inference Trees\n")
    cat("\n")
    cat("Number of trees: ", length(x@ensemble), "\n")
    cat("\n")
    y <- x@responses
    if (y@ninputs > 1) {
        cat("Responses:", paste(names(y@variables),
                                collapse = ", "), "\n")
    }  else {
        cat("Response: ", names(y@variables), "\n")
    }
    inames <- names(x@data@get("input"))
    if (length(inames) > 1) {
        cat("Inputs: ", paste(inames, collapse = ", "), "\n")
    } else {
        cat("Input: ", inames, "\n")
    }
    cat("Number of observations: ", x@responses@nobs, "\n\n")
    invisible(x)
}

setMethod("show", "BinaryTree", function(object) print(object))
setMethod("show", "RandomForest", function(object) print(object))

