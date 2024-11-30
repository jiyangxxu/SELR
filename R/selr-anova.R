#' Internal functions called by other functions.
#' @noRd
#' @name anova.selr
#'
#'
#' @aliases anova.selr
#'
#'
#' @keywords internal

scorefct <- function(o, beta=NULL, testidx=NULL, sas=FALSE) {

  clusters <- unique(o$id)
  nclusters <- length(clusters)
  if (is.null(beta)) {
    beta <- coef(o)
  }

  # Offsets handled correctly?
  y <- o$y

  x <- model.matrix(o)

  linear.predictors <- x%*%beta
  if (!is.null(o$offset))
    linear.predictors <- linear.predictors + o$offset
  mui <- o$family$linkinv(linear.predictors)

  invert <- ginv

  myres <- lapply(clusters, function(cluster) {
    # Individiuals in cluster
    idx <- (o$id == cluster)

    # Cluster size
    csize <- sum(idx)

    # Di is r*k
    #        Di <- t(x[idx,,drop=FALSE]) %*% MESS:::qdiag(o$family$mu.eta(linear.predictors[idx]), nrow=csize)
    Di <- crossprod(x[idx,,drop=FALSE], diag(o$family$mu.eta(linear.predictors[idx]), nrow=csize))
    A  <- diag(sqrt(o$family$variance(mui[idx])), nrow=csize)
    Rmat <- diag(csize)
    Ralpha <- switch(o$corstr,
                     independence = Rmat,
                     exchangeable = matrix(rep(o$geese$alpha, csize^2), csize),
                     ar1 = o$geese$alpha^abs(row(Rmat) - col(Rmat)))
    if (o$corstr=="exchangeable")
      diag(Ralpha) <- 1

    V <- outer(MESS::qdiag(A), MESS::qdiag(A))*Ralpha*o$geese$gamma  # Ok

    # V inverse
    Vinv <- invert(V)
    DiVinv <- tcrossprod(Di, Vinv)
    list(score = DiVinv %*% (y[idx] - mui[idx]),
         DUD   = tcrossprod(DiVinv, Di))
  })

  ## Speed improvement?
  # S <- apply(sapply(myres, function(oo) oo[[1]]), 1, sum)
  S <- rowSums(sapply(myres, function(oo) oo[[1]]))

  Vsand <- Reduce("+", lapply(myres, function(oo) { tcrossprod(oo[[1]])})) # I_1
  VDUD  <- Reduce("+", lapply(myres, function(oo) { oo[[2]] })) #  I_0
  iVDUD <- invert(VDUD)

  if(is.null(testidx)) {
    cat("Should not really be here")
    as.numeric(S %*% invert(invert(VDUD) %*% Vsand %*% invert(VDUD)) %*% S / nclusters)
  } else {
    if (sas) {
      as.numeric(t(S) %*% iVDUD[,testidx] %*% invert( (iVDUD %*%  Vsand  %*% iVDUD)[testidx,testidx] ) %*% iVDUD[testidx,] %*% S)
    }
    else {
      myvar <- Vsand[testidx,testidx] - Vsand[testidx, -testidx] %*% invert(Vsand[-testidx,-testidx]) %*% Vsand[-testidx,testidx]
      as.numeric(t(S[testidx]) %*% invert( myvar ) %*% S[testidx])
    }
  }
}













anovageePrim2 <- function(m1, m2, ..., test=NULL){

  mm1 <- model.matrix(m1)
  mm2 <- model.matrix(m2)

  P1 <- mm1 %*% solve(t(mm1)%*%mm1) %*% t(mm1)
  P2 <- mm2 %*% solve(t(mm2)%*%mm2) %*% t(mm2)
  e2 <- mm2 - P1 %*% mm2
  e1 <- mm1 - P2 %*% mm1

  m2inm1 <- all(apply(e2,2,var) < 1e-10)
  m1inm2 <- all(apply(e1,2,var) < 1e-10)

  if (!any(c(m2inm1,m1inm2)))
    cat("Models not nested\n")
  else
    if (all(c(m2inm1,m1inm2)))
      cat("Models are identical\n")
  else {
    if (m1inm2){
      tmp <- m1
      m1 <- m2
      m2 <- tmp
    }
    ## Now mm2 < mm1
    mm1 <- model.matrix(m1)
    mm2 <- model.matrix(m2)

    mf1 <- paste(paste(formula(m1))[c(2,1,3)],collapse=" ")
    mf2 <- paste(paste(formula(m2))[c(2,1,3)],collapse=" ")

    ## Reparametrize the model
    mm <- cbind(mm2,mm1)
    qmm <- qr(mm)
    qmmq <- qr.Q(qmm)
    nymm1 <- as.data.frame(qmmq[,1:qmm$rank])
    colnames(nymm1) <- paste("parm",1:ncol(nymm1),sep=".")  # bigger model
    nymm2 <- nymm1[,1:ncol(mm2),drop=FALSE]  # smaller model

    formula1 <- formula(paste(formula(m1)[[2]],formula(m1)[[1]],
                              paste(c("-1",colnames(nymm1)),collapse="+"),collapse=""))

    m1call <- m1$call
    m1call$family <- m1$family

    if (m1$family$link != m2$family$link) {
      warning("result might be incorrect due to different link functions")
    }

    ## BUGFIX provided by Stefan Boehringer
    ##nymm1[,paste(formula(m1call)[[2]])] <- m1$y
    nymm1[, paste(formula(m1)[[2]])] <- m1$y

    nymm1[,paste(m1call$id)] <- m1$id
    m1call$offset <- m1$offset
    m1call$weights <- m1$weights
    m1call$formula <- formula1
    m1call$data <- nymm1

    m1ny <- eval(m1call)

    beta  <- coef(m1ny)
    vbeta <- summary(m1ny)$cov.unscaled
    df    <- dim(mm1)[2] - dim(mm2)[2]
    rbeta <- rep(1, length(beta))
    rbeta[1:df] <- 0
    beta0   <- rev(rbeta)
    zeroidx <- beta0 == 0
    ##X2 <- t(beta[zeroidx]) %*% solve(vbeta[zeroidx, zeroidx, drop=FALSE]) %*% beta[zeroidx]
    ##X2 <- t(beta[zeroidx]) %*% solve(vbeta[zeroidx, zeroidx, drop=FALSE], beta[zeroidx])

    if (test=="score"){
      ## Calculate score statistic

      if (any(m1$weights!=1|m2$weights!=1))
      {
        if (any(m1$weights!=1)){
          weights <- m1$weights
          exclude <- which(weights==0)
          call    <- m1$call
          data    <- model.frame(m1)[-exclude,]
          data[,paste(formula(m1)[[2]])] <- m1$y[-exclude]
          data[, paste(call$id)] <- m1$id[-exclude]
          call$offset  <- m1$offset[-exclude]
          call$weights <- c(rep(1,length(m1$y[-exclude])))  # modified
          call$data    <- data
          mod_formula  <- m1$formula
          call$formula <- as.formula(mod_formula)
          call$family  <- m1$family  # added
          mod1 <- eval(call)
        }
        if (any(m2$weights!=1)){
          call <- NULL
          weights <- m2$weights
          exclude <- which(weights==0)
          call    <- m2$call
          data    <- model.frame(m2)[-exclude,]
          data[,paste(formula(m2)[[2]])] <- m2$y[-exclude]
          data[, paste(call$id)] <- m2$id[-exclude]
          call$offset  <- m2$offset[-exclude]
          call$weights <- c(rep(1,length(m2$y[-exclude])))  # modified
          call$data    <- data
          mod_formula  <- m2$formula
          call$formula <- as.formula(mod_formula)
          call$family  <- m2$family
          mod2 <- eval(call)
        }

        m1 <- mod1
        m2 <- mod2
      }

      if ("selr" %in% class(m1) ){
        m1 <- m1$mod
      }

      if ("selr" %in% class(m2) ){
        m2 <- m2$mod
      }

      newbeta <- coef(m1)
      newbeta[which(!zeroidx)] <- coef(m2)
      newbeta[which(zeroidx)] <- 0
      # X2 <- MESS:::scorefct(m1, newbeta, testidx=which(zeroidx))
      X2  <-scorefct(m1, newbeta, testidx=which(zeroidx))

    } else {
      ## Calculate wald statistic
      V0 <- vbeta[zeroidx, zeroidx, drop=FALSE]
      b0 <- beta[zeroidx]
      ##bv <<- list(b0=b0, V0=V0)
      ##X2 <- t(b0) %*% solve(V0, b0)
      X2 <- as.numeric( t(b0) %*% MASS::ginv(V0) %*% b0)
      ev <- eigen(V0, only.values=TRUE)$values
      df.real <- sum(ev > 1e-12)
    }


    ## Make table with results
    topnote <- paste("Model 1", mf1,"\nModel 2", mf2)
    title  <- if(test=="Wald")"Analysis of 'Wald statistic' Table\n"
    else{"Analysis of 'Score statistic' Table\n"}
    table  <- if(test=="Wald"){data.frame(Df=df.real, X2=X2, p=1 - pchisq(X2, df.real))}
    else {data.frame(Df=df, X2=X2, p=1 - pchisq(X2, df))}
    dimnames(table) <- list("1", c("Df", "X2", "P(>|Chi|)"))
    val <- structure(table, heading = c(title, topnote),
                     class = c("anova", "data.frame"))
    return(val)
  }
}








#' @export
anova.selrlist <-
  function (object, ..., dispersion = NULL, test = NULL)
  {
    # object <- object$mod ##
    object <- lapply(object, function(object) object$mod)
    responses <- as.character(lapply(object, function(x) {
      deparse(formula(x)[[2]])
    }))
    sameresp <- responses == responses[1]
    if (!all(sameresp)) {
      object <- object[sameresp]
      warning("Models with response ", deparse(responses[!sameresp]),
              " removed because response differs from ", "model 1")
    }

    ns <- sapply(object, function(x) length(x$residuals))
    if (any(ns != ns[1]))
      stop("models were not all fitted to the same size of dataset")

    objects <- list(object,...)
    m1 <- objects[[1]][[1]]  # ?
    if (length(objects[[1]])>1)
      m2 <- objects[[1]][[2]]
    else
      m2 <- NULL

    value <- anovageePrim2(m1,m2,test=test)
    return(value)
  }






#' @export
anova.selr <- function(object, ..., dispersion = NULL, test=c("score","Wald"))
{
  dotargs <- list(...)
  named <- if (is.null(names(dotargs)))
    rep(FALSE, length(dotargs))
  else (names(dotargs) != "")
  if (any(named))
    warning("The following arguments to anova.glm(..) are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
  dotargs <- dotargs[!named]
  is.glm <- unlist(lapply(dotargs, function(x) inherits(x, "glm")))
  dotargs <- dotargs[is.glm]
  test <- match.arg(test)

  if (length(dotargs) > 0)
    return(anova.selrlist(c(list(object), dotargs), dispersion = dispersion,
                          test = test))

  object <- object$mod  ##
  varlist <- attr(object[["terms"]], "variables")
  ##print(varlist)
  x <- if (n <- match("x", names(object), 0)) {
    object[[n]]
  }
  else {
    model.matrix(object)
  }

  varseq   <- attr(x, "assign")
  nvars    <- max(0, varseq)  # highest hierarchy
  xnames   <- attr(object[["terms"]], "term.labels")  # added
  indices  <- sapply(unique(xnames), function(name) which(xnames==name)[1])  # added
  betaList <- vbetaList <- modelList <- NULL  # revised

  if (nvars > 1) {
    method <- object[["method"]]
    if (!is.function(method))
      method <- get(method, mode = "function", envir = parent.frame())
    for (i in 1:nvars) {  # revised
      ##eprint("calling fit....")
      ##print(length(object$y))
      fit <- method(x = x[, varseq <= i, drop = FALSE],
                    y = object[["y"]], weights = object$prior.weights,
                    corstr = object[["corstr"]],
                    start  = object$start, offset = object$offset, id=object$id,
                    family = object$family, control = object$control)

      # added

      y<-paste(as.character(formula(object)[[2]]),"~-1+")

      mod_formula<-formula(paste(y,paste(names(indices)[indices<=i],collapse="+")))

      weights <- object$weights
      exclude <- which(weights==0)

      call <- object$call

      if(length(exclude)>0){
        df <- model.frame(object)[-exclude,]
        df[,paste(formula(object)[[2]])] <- object$y[-exclude]
        df[, paste(call$id)] <- object$id[-exclude]
        call$offset  <- object$offset[-exclude]
        call$weights <- c(rep(1,length(object$y[-exclude])))

      }  # modified
      else{
        df <- model.frame(object)
        df[,paste(formula(object)[[2]])] <- object$y
        df[, paste(call$id)] <- object$id
        call$offset  <- object$offset
        call$weights <- c(rep(1,length(object$y)))
      }
      call$data    <- df
      call$formula <- as.formula(mod_formula)
      call$family  <- object$family  # added

      mod <- eval(call)

      betaList <- c(betaList, list(fit$beta))
      vbetaList <- c(vbetaList, list(fit$vbeta))
      modelList <- c(modelList, list(mod))  # added
    }
  }

  # betaList <- c(betaList, list( object$geese$beta ))
  # vbetaList <- c(vbetaList, list( object$geese$vbeta ))
  # modelList <- c(modelList, list(object))

  hasIntercept <- (length(grep("(Intercept)", names(betaList[[1]]))) != 0)
  dimVec <- unlist(lapply(betaList, length))  # dimensions for each "hierarchy"

  if (hasIntercept){
    dfVec <- 1
  } else {
    dfVec <- c(dimVec[[1]])
  }

  if (length(dimVec) > 1){
    for (i in 1:length(dimVec))
      dfVec <- c(dfVec, dimVec[i] - dimVec[i-1])
  }  # dimensions for each variable

  X2Vec <- NULL

  for (i in 2:length(dfVec)){

    beta    <- betaList[[i]]
    vbeta   <- vbetaList[[i]]
    beta0   <- rep(1, length(beta))
    beta0[1:dfVec[i]] <- 0
    beta0   <- rev(beta0)
    zeroidx <- beta0 == 0

    if (test == "score"){  # added

      ## Calculate Score statistics
      newbeta <- modelList[[i]]$coefficients
      newbeta[which(!zeroidx)] <- betaList[[i-1]]
      newbeta[which(zeroidx)] <- 0

      ##print(newbeta)

      ## print(modelList[[i]])
      X2     <-scorefct(modelList[[i]], newbeta, testidx=which(zeroidx))

      ##print(X2)
      X2Vec  <-c(X2Vec,X2)

    } else {

      ## Calculate Wald statistics
      X2      <- t(beta[zeroidx]) %*% solve(vbeta[zeroidx, zeroidx, drop=FALSE]) %*% beta[zeroidx]
      X2Vec   <- c(X2Vec, X2)
    }
  }

  resdf <- dfVec[-1]  # modified
  resdev <- X2Vec

  # print(resdf)
  # print(resdev)
  tab <- data.frame(resdf, resdev, 1-pchisq(resdev,resdf))
  colnames(tab) <- c("Df", "X2", "P(>|Chi|)")

  tl <- attr(object$terms, "term.labels")[-1]

  if (length(tl) == 0)
    tab <- tab[1, , drop = FALSE]

  if (length(tl))
    rownames(tab) <- c(tl)

  title <- if (test=="Wald") {
    paste("Analysis of 'Wald statistic' Table",
          "link: ", object$family$link,
          "\nResponse: ", as.character(varlist[-1])[1],
          "\nTerms added sequentially (first to last)\n",
          sep = "")
  } else {
    paste("Analysis of 'Score statistic' Table",
          "link: ", object$family$link,
          "\nResponse: ", as.character(varlist[-1])[1],
          "\nTerms added sequentially (first to last)\n",
          sep = "")
  }

  structure(tab, heading = title, class = c("anova", "data.frame"))
}












