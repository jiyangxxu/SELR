#' Internal functions called by other functions.
#' @noRd
#' @name selr-internal
#'
#'
#'
#' @aliases preprocess risk target convertW convertY
#'
#'
#' @keywords internal
preprocess <- function(X,y,Y,W) {
  # # check parameter type
  # if(!(is.data.frame(X))){
  #   stop("bad input for argument 'X'")
  # }
  # if(!(typeof(y)=="list")){
  #   stop("bad input for argument 'y'")
  # }
  # if(!(is.matrix(W))){
  #   stop("bad input for argument 'W'")
  # }
  # if(!(is.matrix(Y))){
  #   stop("bad input for argument 'Y'")
  # }

  # check matrix size
  if(length(unique(unlist(y)))*(length(unique(unlist(y)))-1)!=length(W)){
    stop("incorrect size for argument 'W'")
  }
  if(length(unique(unlist(y)))*(length(unique(unlist(y)))-1)!=length(Y)){
    stop("incorrect size for argument 'Y'")
  }
  if(nrow(X)!=length(unlist(y))){
    stop("incorrect match for argument 'X' and 'y'")
  }

  # check binary coding of W,Y
  if(any(W!=0 & W!=1)){
    stop("bad input for argument 'W'")
  }
  if(any(Y!=0 & Y!=1)){
    stop("bad input for argument 'Y'")
  }

  # subset
  if(any(W==0 & Y==1)){
    stop("'Y' must be subset of 'W'")
  }

  # add index of a sample
  X[["id"]] <- rep(1:nrow(X))

  # convert y to ordinal numbers for comparisons later
  if (all(sapply(y, class) == "factor")){
    y <- as.character(y)
  }

  y_star <- y
  if (all(sapply(y, class) == "integer")) {
    # not sorted from 1
    if (y[order(unlist(y))][1]!=1) {
      y_star <- as.list(match(unlist(y), unique(sort(unlist(y)))))
    }
  } else if (all(sapply(y, class) == "character")) {
    # sort based on input order
    y_star <- as.list(match(unlist(y), unique(unlist(y))))
  }
  y_df <- as.data.frame(cbind(y_star, y=unlist(y)))


  ret <- data.frame()
  for(i in seq_len(ncol(W))){

    temp <- data.frame()

    for(j in seq_len(nrow(W))){

      index <- which(y_df$y_star == j)
      jdf <- X[index,]
      jdf[["y"]] <- y_df[index,][["y"]]
      jdf[["response"]] <- NA

      if(W[j,i]==1&Y[j,i]==1){  # W_{j,i} and Y_{j,1} both are 1, weight = 1
        jdf[["response"]] <- 1
        jdf[["weights"]] <- 1
      } else if(W[j,i]==1&Y[j,i]==0){  # W_{j,i}=1 but Y_{j,1}=0, weight = 1
        jdf[["response"]] <- 0
        jdf[["weights"]] <- 1
      } else if(W[j,i]==0&Y[j,i]==0){  # W_{j,i}=0, weight = 0
        jdf[["response"]] <- 0
        jdf[["weights"]] <- 0
      }

      temp <- rbind(temp,jdf)
    }

    temp[["levmod"]] <- i
    ret <- rbind(ret,temp)
    ret<-ret[order(ret[["id"]]), ]
    ret[["levmod"]]<-as.character(ret[["levmod"]])

  }

  ret[["y"]] <- NULL
  return(ret)
}
















# helper function for constructing the risk sets of three standard models

risk <- function(level, ratio="nonstandard",reverse=F) {
  if (level<=2) {
    stop("invalid input for 'level'")
  }

  ret <- matrix(0,level,level-1)  # initialization

  switch(ratio,
         cumulative=ret<-matrix(1,level,level-1),
         adjacent = {
           for (i in 1:(level-1)) {
             ret[i,i]   <- 1
             ret[i+1,i] <- 1
           }
         },
         continuation = {
           if (!reverse) {
             for (i in 1:(level-1)) {
               for (j in 1:(i+1)) {
                 ret[j,i] <- 1
               }
             }
           } else {
             for (i in 1:(level-1)) {
               for (j in i:level) {
                 ret[j,i] <- 1
               }
             }
           }
         },
         sequential = {
           if (!reverse) {
             for (i in 1:(level-1)) {
               for (j in i:level) {
                 ret[j,i] <- 1
               }
             }
           } else {
             for (i in 1:(level-1)) {
               for (j in 1:(i+1)) {
                 ret[j,i] <- 1
               }
             }
           }
         },
         baseline = {
           if (!reverse){
             for (i in 1:(level-1)){
               ret[i,i] <- 1
               ret[level,i] <- 1
             }
           } else {
             for (i in 1:(level-1)){
               ret[1,i] <- 1
               ret[i+1,i] <- 1
             }
           }
         },
         twolevel = {
           m <- 0
           if ((level %% 2) != 0) {
             m <- round((level+1)/2)
           } else {
             m <- level/2
           }

           if (level==3) {
             ret[1:2,1] <- 1
             ret[1:level,2] <- 1
             return (ret)
           }
           ret[1:m, 1:(m-1)] <- 1
           ret[(m+1):level, (m+1):(level-1)] <- 1
           ret[,m] <- 1
         },
         binarysplit = {
           if ((level %% 2) != 0) {  # odd
             m <- round((level+1)/2)
             ret[,m] <- 1
             ret[-m,(m-1)] <- 1

             if (level == 3) {
               return(ret)
             }

             for (i in 1:(m-2)) {
               ret[1:(i+1),i] <- 1
             } # left

             for (i in (m+1):(level-1)){
               ret[i:level,i] <- 1
             } # right

           } else {  # even
             m <- level/2
             ret[,m] <- 1
             ret[m+1,m+1] <- 1
             ret[m,m+1] <- 1
             ret[-m,(m-1)] <- 1
             ret[m+1,(m-1)] <- 0

             if (level == 4) {
               return (ret)
             }

             for (i in 1:(m-2)) {
               ret[1:(i+1),i] <- 1
             } # left

             for (i in (m+2):(level-1)){
               ret[i:level,i] <- 1
             } # right

           }
         },
         nonstandard = ret
  )

  return(ret)
}











# helper function for constructing the target sets of three standard models

target <- function(level, ratio="nonstandard",reverse=F) {
  if (level<=2) {
    stop("invalid input for 'level'")
  }

  ret <- matrix(0,level,level-1)  # initialization

  switch(ratio,
         cumulative = {
           if (!reverse) {
             for (i in 1:(level-1)) {
               for (j in (i+1):level) {
                 ret[j,i] <- 1
               }
             }
           } else {
             for (i in 1:(level-1)) {
               for (j in i:(i+1)) {
                 ret[j,i] <- 1
               }
             }
           }
         },
         adjacent = {
           if (!reverse) {
             for (i in 1:(level-1)) {
               ret[i+1,i] <- 1
             }
           } else {
             for (i in 1:(level-1)) {
               ret[i,i] <- 1
             }
           }
         },
         continuation = {
           if (!reverse) {
             for (i in 1:(level-1)) {
               ret[i+1,i] <- 1
             }
           } else {
             for (i in 1:(level-1)) {
               ret[i,i] <- 1
             }
           }
         },
         sequential = {
           if (!reverse) {
             for (i in 1:(level-1)) {
               for (j in (i+1):level) {
                 ret[j,i] <- 1
               }
             }
           } else {
             for (i in 1:(level-1)) {
               for (j in 1:i) {
                 ret[j,i] <- 1
               }
             }
           }
         },
         baseline = {
           if (!reverse){
             for (i in 1:(level-1)){
               ret[level,i] <- 1
             }
           } else {
             for (i in 1:(level-1)){
               ret[1,i] <- 1
             }
           }
         },
         twolevel = {
           m <- 0
           if ((level %% 2) != 0) {
             m <- round((level+1)/2)
           } else {
             m <- level/2
           }

           if (level==3) {
             ret[2, 1] <- 1
             ret[3, 2] <- 1
             return (ret)
           }
           for (i in 1:(m-1)){
             for (j in (i+1):m){
               ret[j,i] <- 1
             }
           } # left

           for (i in m:(level-1)){
             for (j in (i+1):level){
               ret[j,i] <- 1
             }
           } # right
         },
         binarysplit = {
           if ((level %% 2) != 0) { # odd
             m <- round((level+1)/2)
             ret[m,m] <- 1
             ret[(m+1):level,(m-1)] <- 1

             if (level==3) {
               return (ret)
             }

             for (i in 1:(m-2)) {
               ret[i+1,i] <- 1
             } # left

             for (i in (m+1):(level-1)){
               ret[i,i] <- 1
             } # right
           } else {  # even
             m <- level/2
             ret[m,m] <- 1
             ret[m+1,m] <- 1
             ret[m+1,m+1] <- 1
             ret[(m+2):level,(m-1)] <- 1

             if (level==4) {
               return (ret)
             }

             for (i in 1:(m-2)) {
               ret[(i+1),i] <- 1
             } # left

             for (i in (m+2):(level-1)){
               ret[i,i] <- 1
             } # right

           }
         },
         nonstandard = ret
  )

  return(ret)
}








convertW <- function(Z) {  # new
  if (any(Z!=-1&Z!=0&Z!=1)) {
    stop("bad input for argument 'Z'")
  }

  W <- 1*(Z==1|Z==0)
  return(W)
}


convertY <- function(Z) {  # new
  if (any(Z!=-1&Z!=0&Z!=1)) {
    stop("bad input for argument 'Z'")
  }

  Y <- 1*(Z==1)
  return(Y)
}




Amat <- function(W_revised,Y_revised,predicts_slice){
  lev <- nrow(W_revised)
  predicts_slice <- c(predicts_slice,0)
  predicts_slice <- matrix(predicts_slice,nrow=lev,ncol=1)
  A <- matrix(0,nrow=lev,ncol=lev)

  for (i in 1:lev){ # row
    for (j in 1:lev){ # col
      if (W_revised[i,j]==1 & Y_revised[i,j]==1){
        A[i,j] <- 1 - predicts_slice[i]
      } else if (W_revised[i,j]==1 & Y_revised[i,j]==0){
        A[i,j] <- predicts_slice[i] * (-1)
      } else if (W_revised[i,j]==0 & Y_revised[i,j]==0){
        A[i,j] <- 0
      }
    }
  }

  A[nrow(A),] <- 1

  return (A)
}




