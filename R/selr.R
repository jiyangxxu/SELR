#' @title Simultaneous and Efficient Logistic Regression (SELR)
#'
#' @description The selr function fits simultaneous and efficient logistic
#'      regressions using the 'selr' function of the 'selr' package.
#'      The underlying function is 'geeglm' in 'Geepack' library.
#'      anova "Wald" and "score" are available.
#'
#' @name selr
#'
#'
#' @param formula See corresponding documentation to \code{glm}
#' @param data See corresponding documentation to \code{glm}
#' @param subset See corresponding documentation to \code{glm}
#' @param na.action See corresponding documentation to \code{glm}
#' @param contrasts See corresponding documentation to \code{glm}
#' @param offset See corresponding documentation to \code{glm}
#' @param link See corresponding documentation to \code{geeglm}.
#'     The following are allowed: \code{'identity'}, \code{'logit'},
#'     \code{'probit'}, \code{'cloglog'}. The default link is \code{'logit'}.
#' @param modelType a character string specifying the ratio value.
#'     The following are permitted: \code{'cumulative'},
#'     \code{'adjacent'}, \code{'continuation'}, \code{'sequential'},
#'     \code{'baseline'}, \code{'twolevel'}, \code{'binarysplit'},
#'     \code{'nonstandard'}. Note that if the ratio is \code{'nonstandard'},
#'     specify it through \code{'risk'} and \code{'target'} or 'Z' below.
#'     'risk' and 'target' or 'Z' should not be assigned if ratio is not
#'     \code{'nonstandard'}. Default value is \code{'nonstandard'}.
#' @param risk a user defined \eqn{S*(S-1)} matrix, also
#'      known as 'risk set'.
#' @param target a user defined \eqn{S*(S-1)} matrix, also
#'      known as 'target set'.
#' @param Z a user defined \eqn{S*(S-1)} matrix. 'Z' is using one
#'      matrix to represent both \code{'risk'} and \code{'target'}. Please note
#'      that 'Z' should not be assigned if \code{'risk'} and \code{'target'} are assigned.
#' @param idmod a string of column names of \code{data} clarifying the clustering
#' @param parallel a logical variable. If parallel is true,
#'     selr will compute parallelism case for \code{'levmod'}. Otherwise,
#'     selr will compute nonparallelism case. Default is true.
#' @param reverse a logical variable. Please note that reverse
#'     only apply to ratios of \code{'cumulative'}, \code{'adjacent'},
#'     \code{'continuation'}, \code{'sequential'}, and \code{'baseline'}.
#'     Default is false.
#' @param ordinal a logical variable, indicating whether fitted data is ordinal
#'     or not. Default is true.
#' @param yord a list contained the desired order of response levels used for
#'     comparison from \code{risk} and \code{target}. Please note that
#'     this argument is not necessary to be assigned if ordinal is true.
#' @param elr a logical variable, indicating whether using efficient logistics
#'     regression or not. Default is false, which is using the simultaneous
#'     logistics regression.
#' @param glmnet a \code{glmnet} object. If not null, the model would return an
#'     additional \code{glmnet} model. See more from \code{glmnet}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{selr} returns a list contains three values: a selr object,
#'     the ratio, the value of Z matrix, the original input data set, and
#'     a \code{glmnet} object if the user requests a \code{glmnet} call
#'
#' @section Warning: \code{'elr'} should be used with more care. If the response
#'     is not numeric class but ordinal, please sort the data set based desired
#'     order of response before using \code{'selr'}.
#'
#' @note Default models are \eqn{\text{link}{\Pr(Y_i\geq s)}}, \eqn{\text{link}{\Pr(Y_i=s|Y_i=s-1)}},
#'    \eqn{\text{link}{\Pr(Y_i=s+1|Y_i\leq s+1)}}, \eqn{\text{link}{\Pr(Y_i\leq s|Y_i\leq s+1)}}
#'    for 'cumulative', 'adjacent', 'continuation', 'sequential' ratios.
#'    If reverse is true, models are \eqn{\text{link}{\Pr(Y_i\geq s)}},
#'    \eqn{\text{link}{\Pr(Y_i=s|Y_i=s+1}}, \eqn{\text{link}{\Pr(Y_i=s|Y_i\geq s)}},
#'    \eqn{\text{link}{\Pr(Y_i\geq s|Y_i\geq s-1)}} with the same order.
#'
#' @author Jiyang Xu \email{jiyangx3@@illinois.edu},
#'         Douglas Simpson \email{dgs@@illinois.edu}
#'
#' @seealso \code{\link{geeglm}}, \code{\link{glm}}, \code{\link{anova.glm}}
#'
#' @references
#' Fu, L., & Simpson, D. G. (2002).
#' Conditional risk models for ordinal response data:
#' simultaneous logistic regression analysis and generalized score tests.
#' Journal of Statistical Planning and Inference, 108(1–2), 201–217.
#' https://doi.org/10.1016/s0378-3758(02)00279-3
#'
#' Halekoh, U., Højsgaard, S., & Yan, J. (2006).
#' The R Package geepack for Generalized Estimating Equations.
#' Journal of Statistical Software, 15(2). https://doi.org/10.18637/jss.v015.i02
#'
#' Ekstrøm, C. T. (2022, June 20).
#' Miscellaneous Esoteric Statistical Scripts [R package MESS version 0.5.9].
#' https://cran.r-project.org/web/packages/MESS
#'
#'
#' @keywords models
#' @examples
#'
#' data(mntlhlth)
#' head(mntlhlth)
#'
#' mod1 <- selr(mntlhlth~socioeconomic,data=mntlhlth,modelType="cumulative",parallel=FALSE)
#' summary(mod1)
#'
#' mod2 <- selr(mntlhlth~socioeconomic,data=mntlhlth,modelType="cumulative")
#' summary(mod2)
#'
#' anova(mod1,test="score")
#' anova(mod1,mod2,test="score")
#'
#' QIC(mod1)
#' QIC(mod1,mod2)
#'
#' mod3 <- selr(mntlhlth~socioeconomic,data=mntlhlth,modelType="continuation",parallel=FALSE,elr=TRUE)
#' summary(mod3)
#'
#' mod4 <- selr(mntlhlth~socioeconomic,data=mntlhlth,modelType="continuation",elr=TRUE)
#' summary(mod4)
#'
#' anova(mod4,test="score")
#' anova(mod3,mod4,test="Wald")
#'
#' QIC(mod3)
#' QIC(mod3,mod4)
#'
#' # non-standard 1 (Y,W)
#' # {1} vs. {2}, {1,2} vs. {3,4}, {3} vs. {4}
#' W <- matrix(c(1,1,0,0,1,1,1,1,0,0,1,1), nrow=4, ncol=3)
#' W
#'
#' Y <- matrix(c(0,1,0,0,0,0,1,1,0,0,0,1), nrow=4, ncol=3)
#' Y
#'
#' Z <- matrix(c(0,1,-1,-1,0,0,1,1,-1,-1,0,1), nrow=4, ncol=3)
#' Z
#'
#' mod5 <- selr(mntlhlth~socioeconomic,data=mntlhlth,risk=W,target=Y)
#' summary(mod5)
#'
#' mod6 <- selr(mntlhlth~socioeconomic,data=mntlhlth,Z=Z)
#' summary(mod6)
#'
#'
#' @export selr
selr <- function(formula,data,subset=NULL,na.action, offset, contrasts=NULL,
                 link=c("logit","identity","probit","cloglog"),
                 modelType=c("nonstandard","cumulative","adjacent","continuation",
                             "baseline","sequential","twolevel","binarysplit"),
                 risk=NULL,target=NULL,Z=NULL,idmod=NULL,parallel=T,reverse=F,ordinal=T,
                 yord=NULL,elr=F,glmnet=NULL,...){

  modelType <- match.arg(modelType)
  link <- match.arg(link)

  ratio <- modelType

  if (ratio=="nonstandard"&is.null(Z)&(is.null(risk)|is.null(target))){
    stop("choose appropriate risk and target sets")
  }  # new

  if (is.null(Z)==F&is.null(risk)==F&is.null(target)==F){
    warning("unnecessary to assign the arguments at the same time")
  }  # new

  if (!is.null(ratio)){
    if(ratio!="nonstandard"&is.null(Z)==F&is.null(risk)==F&is.null(target)==F){
      stop("standard and non-standard cannot be assigned together")}
  }

  if (reverse&(ratio!="cumulative")&(ratio!="adjacent")&(ratio!="continuation")&
      (ratio!="sequential")&(ratio!="baseline")) {
    warning("reverse is not applicable to input ratio")
  }

  # if (!is.null(yord)){
  #   if (length(yord)!=length(as.list(data[names(data)==as.character(formula[[2]])]))){
  #     stop("input value for 'yord' is invalid")
  #   }
  # }

  if (!ordinal&is.null(yord)){
    warning("please sort the data based on desired orders of 'risk' and 'target' sets")
  }

  # if (ordinal&(!is.null(yord))){
  #   warning("argument 'yord' is not designed for ordinal data")
  # }

  mf <- match.call(expand.dots = TRUE)
  m <- match(c("formula", "data", "subset", "na.action", "offset", "idmod"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  idmod_name <- deparse(substitute(idmod))

  if ("(idmod)" %in% names(mf)) {
    names(mf)[names(mf) == "(idmod)"] <- idmod_name
  }

  if (idmod_name %in% names(mf)) {
    idmod <- mf[[idmod_name]]
  }

  mt <- attr(mf, "terms")

  if(!is.null(yord)){
    mf <- mf[order(match(mf[[all.vars(formula(mt))[1]]], yord)), ]
  }

  y <- model.response(mf)
  X <- model.matrix(mt, mf)

  # colnames(mf) <- gsub("as.factor\\((.*)\\)", "\\1", colnames(mf))

  if (is.matrix(y)){
    ylist <- colnames(y)
    mf <- data
    mf <- pivot_longer(mf, cols = all_of(ylist), names_to="response",
                       values_to = "value")
    mf <- mf[mf$value==1,]
    mf$value <- NULL
    mf <- mf[order(match(mf$response, ylist)),]
    y <- mf$response
    yname <- "response"
    formula <- formula(mt)
    X_name <- all.vars(formula)
    X_name <- setdiff(X_name, ylist)
  } else {
    ylist <- unique(y)
    formula <- formula(mt)
    yname <- all.vars(formula)[1]
    X_name <- all.vars(formula)[-1]
    t <- gsub(":", "*", attr(terms(formula), "term.labels"))
  }

  terms_info <- attr(terms(formula), "term.labels")
  factor_terms <- grep("as.factor", terms_info, value = TRUE)

  for (factor_term in factor_terms) {
    factor_var <- gsub("as.factor\\((.+)\\)", "\\1", factor_term)
    X_name[X_name == factor_var] <- factor_term
  }

  if (!is.null(idmod)) {
    X_name <- c(X_name, idmod_name)
  }

  level<-length(ylist)
  if(ratio!="nonstandard"){  # standard models
    risk<-risk(level=level,ratio=ratio,reverse=reverse)
    target<-target(level=level,ratio=ratio,reverse=reverse)
  }

  if(!is.null(Z)){
    risk<-convertW(Z)
    target<-convertY(Z)
  }

  # print(target)
  # print(risk)


  df_x <- mf[,X_name,drop=F]
  ##print(df_x)

  ####have to modify so that it could evaluate outside the function??
  df <- preprocess(df_x,y,Y=target,W=risk)

  ##print(colnames(df))

  colnames(df)[colnames(df)=="response"]<-yname
  df$levmod <- as.factor(df$levmod)

  terms_info <- attr(terms(as.formula(formula)), "term.labels")
  factor_terms <- grep("as.factor", terms_info, value = TRUE)

  for (factor_term in factor_terms) {
    factor_var <- gsub("as.factor\\((.+)\\)", "\\1", factor_term)
    if (factor_var %in% colnames(df)) {
      df[[factor_var]] <- as.factor(df[[factor_var]])
    }
  }

  colnames(df) <- gsub("as.factor\\((.+)\\)", "\\1", colnames(df))


  if (parallel==T)
  {formula<-as.formula(paste(yname,"~-1+levmod+", formula)[[3]])}
  else # note: x1+x2*x3->levmod*x1+levmod*x2*x3
  {#list    <- trimws(unlist(strsplit(deparse(formula[[3]]), "\\+")))
    formula <- as.formula(paste(yname,"~ -1+",paste("levmod*", t, collapse = "+")))}  # revised

  weights <- df$weights
  family <- substitute(binomial(link = link))
  lev <- nrow(risk)

  if (!elr)
  {
    if (is.null(idmod))
    {
      call <- as.call(c(quote(geeglm), formula = formula, family = family,
                       data = quote(df), weights = quote(weights), id = quote(id),
                       corstr = quote("independence")))}
    else
    {

      # idmod = df[[idmod]]
      # correlation_mat = diag(length(unique(interaction(idmod,df[["id"]]))))

      # call <- as.call(c(quote(geeglm), formula = formula, family = family,
      #                  data = quote(df), weights = quote(weights), id = quote(interaction(idmod,id)),
      #                  corstr = quote("userdefined"), zcor = quote(correlation_mat)))

      call <- as.call(c(quote(geeglm), formula = formula, family = family,
                        data = quote(df), weights = quote(weights), id = quote(eval(parse(text = idmod_name))),
                        corstr = quote("independence")))
    }
  }
  else if (elr)
  {
    if (level == 3){
      if (is.null(idmod))
      {call <- as.call(c(quote(geeglm), formula = formula, family = family,
                         data = quote(df), weights = quote(weights), id = quote(id),
                         corstr = quote("exchangeable")))}
      else
      {
        # idmod = df[[idmod]]
        # correlation_mat = diag(length(unique(interaction(idmod,df[["id"]]))))

        call <- as.call(c(quote(geeglm), formula = formula, family = family,
                         data = quote(df), weights = quote(weights), id = quote(eval(parse(text = idmod_name))),
                         corstr = quote("exchangeable")))}
    } else {
      if (is.null(idmod))
      {call <- as.call(c(quote(geeglm), formula = formula, family = family,
                         data = quote(df), weights = quote(weights), id = quote(id),
                         corstr = quote("unstructured")))}
      else
      {
        # idmod = df[[idmod]]
        # correlation_mat = diag(length(unique(interaction(idmod,df[["id"]]))))

        call <- as.call(c(quote(geeglm), formula = formula, family = family,
                         data = quote(df), weights = quote(weights), id = quote(eval(parse(text = idmod_name))),
                         corstr = quote("unstructured")))
        }
    }
  }
  mod <- eval(call, envir = df)

  ##print(risk)
  ##print(target)

  ##added output values
  if(is.null(Z)){
    Z<-matrix(0,level,level-1)
    for(i in seq_len(ncol(risk))){
      for(j in seq_len(nrow(risk))){
        if(risk[j,i]==1&target[j,i]==1){
          Z[j,i]<- 1
        } else if(risk[j,i]==1&target[j,i]==0){
          Z[j,i]<- 0
        } else if(risk[j,i]==0&target[j,i]==0){
          Z[j,i]<- -1
        }
      }
    }
  }

  comparison <- vector("list",level-1)
  for(i in 1:ncol(Z)){
    comparison[[i]] <- paste0("levmod", i)
  }
  comparison <- unlist(comparison)

  rownames(Z) <- ylist
  colnames(Z) <- comparison

  # revised glmnet
  matrix_data <- model.matrix(mod)
  levmod_indices <- grep("levmod", colnames(matrix_data))
  selected_columns <- matrix_data[, levmod_indices]  # pivot wider levmods
  df[["levmod"]] <- NULL
  df <- cbind(df, selected_columns)

  if(!is.null(glmnet)){
    call <- match.call()
    call <- call[["glmnet"]]
    addargs <- as.list(call)[-1]

    df[["weights"]] <- NULL
    c <- list("id")
    X <- df[,!(colnames(df) == yname)]
    idx <- match(c, colnames(X))
    levmod_indices <- grep("levmod", colnames(X))
    vec <- rep(1,ncol(X))
    vec[idx] <- 0
    vec[levmod_indices] <- 0

    args <- list(
      x = X,
      y = df[, yname],
      family = "binomial",
      weights = df[["weights"]],
      penalty.factor = vec,
      intercept = F
    )

    addargs[names(addargs) %in% names(args)] <- NULL
    addargs <- addargs[names(addargs) != ""]
    addargs <- addargs[!sapply(addargs, is.null)]
    args <- c(args, addargs)
    glmnet <- do.call(glmnet::glmnet, args)
  }

  if (!is.null(glmnet))
  {x<-list(mod=mod,modelType=modelType,Z=Z,data=data,glmnet=glmnet)}
  else
  {x<-list(mod=mod,modelType=modelType,Z=Z,data=data)}
  attr(x, "class") <- c("selr","geeglm","glm", "lm")
  return (x)
}







#' @export
summary.selr <- function(object,...){

  class(object) <- "selr"

  value <- NULL
  value$call <- object$mod$call
  covmat <- object$mod$geese$vbeta
  value$cov.scaled   <-   value$cov.unscaled <-  covmat

  mean.sum  <- data.frame(estimate = object$mod$geese$beta, std.err=sqrt(diag(covmat)))
  mean.sum$wald <- (mean.sum$estimate / mean.sum$std.err)^2
  mean.sum$p <- 1 - pchisq(mean.sum$wald, df=1)
  names(mean.sum) <- c("Estimate", "Std.err", "Wald", "Pr(>|W|)")
  value$coefficients <- mean.sum
  covmatgam <- object$mod$geese$vgamma
  scale.sum  <- data.frame(Estimate = object$mod$geese$gamma, Std.err=sqrt(diag(covmatgam)))
  if (!is.null(object$mod$geese$zsca.names)) rownames(scale.sum) <- object$mod$geese$zsca.names
  value$dispersion   <-  scale.sum

  covmatalpha <- object$mod$geese$valpha
  corr.sum <- data.frame(Estimate = object$mod$geese$alpha, Std.err=sqrt(diag(covmatalpha)))
  value$corr <- corr.sum

  value$corstr    <- object$mod$geese$model$corstr
  value$scale.fix <- object$mod$geese$model$scale.fix
  value$cor.link  <- object$mod$geese$model$cor.link
  value$clusz <- object$mod$geese$clusz
  value$error <- object$mod$geese$error

  value$modelType <- object$modelType
  value$Z <- object$Z
  class(value) <- "summary.selr"

  return(value)
}





#' @export
print.summary.selr <- function (x,
                                digits = max(3, getOption("digits") - 3),
                                quote = FALSE, prefix = "", ...)
{
  xg <- x$geese
  if (is.null(digits))
    digits <- options()$digits
  else options(digits = digits)
  cat("\nCall:\n");   print(x$call)
  cat("\n Coefficients:\n");
  ##print(as.matrix(x$coef), digits = digits)
  printCoefmat(as.matrix(x$coef), digits = digits)

  cat("\nCorrelation structure =", x$corstr, "\n")
  if (!(x$scale.fix)) {
    cat("Estimated Scale Parameters:\n\n")
    print(x$dispersion[1:2], digits = digits)
  }
  else cat("Scale is fixed.\n\n")

  if (pmatch(x$corstr, "independence", 0) == 0) {
    cat("  Link =", x$cor.link, "\n")
    cat("\nEstimated Correlation Parameters:\n")
    print(x$corr, digits = digits)
  }

  cat("Number of clusters:  ", length(x$clusz),
      " Maximum cluster size:", max(x$clusz), "\n", "\n")

  cat("modelType: ", x$modelType, "\n")
  cat("Z matrix:\n")

  Z_matrix <- as.matrix(x$Z)
  col_names <- colnames(Z_matrix)

  max_col_widths <- apply(Z_matrix, 2, function(x) max(nchar(as.character(x))))
  max_row_name_width <- max(nchar(rownames(Z_matrix)))
  max_name_width <- max(max(max_col_widths), max_row_name_width) + 2

  center_pad <- function(text, width) {
    pad_left <- floor((width - nchar(text)) / 2)
    pad_right <- width - nchar(text) - pad_left
    if (pad_left < 0) {
      pad_left <- 0
    }
    if (pad_right < 0) {
      pad_right <- 0
    }
    paste0(strrep(" ", pad_left), text, strrep(" ", pad_right))
  }

  cat(center_pad("", max_row_name_width))
  for (col_name in col_names) {
    cat(center_pad(col_name, max_name_width))
  }
  cat("\n")

  for (i in 1:nrow(Z_matrix)) {
    cat(center_pad(rownames(Z_matrix)[i], max_row_name_width))

    for (j in 1:ncol(Z_matrix)) {
      cat(center_pad(as.character(Z_matrix[i, j]), max_name_width))
    }
    cat("\n")
  }



  invisible(x)
}





#' @export
print.selr <- function(x, digits = NULL, quote = FALSE, prefix = "", ...){
  xg <- x$geese
  if (is.null(digits))
    digits <- options()$digits
  else options(digits = digits)


  cat("\nCall:\n")
  print(x$mod$call)
  cat("\nCoefficients:\n")
  print(unclass(x$mod$coefficients), digits = digits)

  cat("\nDegrees of Freedom:", length(x$y), "Total (i.e. Null); ",
      x$df.residual, "Residual\n")

  cat("modelType: ", x$modelType)

}












#' @title Predictions for SELR
#'
#' @description Predict function for \code{selr} model.
#'
#' @name predict.selr
#'
#'
#' @param object a \code{selr} object
#' @param newdata a data set containing the explanatory variables
#' @param type a logical variable, available options are \code{'link'}, \code{'response'}
#' @param marginal a logical variable, indicating whether marginal probability
#' should be calculated. If not, conditional risk probbaility would be returned.
#'
#'
#' @return \code{predict.selr} returns a data set containing the predicted values.
#'
#'
#'
#' @author Jiyang Xu \email{jiyangx3@@illinois.edu},
#'         Douglas Simpson \email{dgs@@illinois.edu}
#'
#' @seealso \code{\link{predict.glm}}
#'
#'
#' @keywords predict
#' @examples
#'
#' data(mntlhlth)
#' head(mntlhlth)
#'
#' set.seed(123)
#' mntlhlth$mntlhlth <- factor(mntlhlth$mntlhlth, levels = c("well", "mild", "moderate", "impaired"))
#' mntlhlth$socioeconomic = as.numeric(mntlhlth$socioeconomic)
#' itrain=sample(seq(nrow(mntlhlth)),floor(0.8*nrow(mntlhlth)),replace=FALSE)
#' train=mntlhlth[itrain,]
#' test=mntlhlth[-itrain,]
#' ord <- match(train$mntlhlth, list(1,2,3,4))
#' train <- train[order(ord),]
#' train <- train[order(train$socioeconomic),]
#' ord <- match(test$mntlhlth, list(1,2,3,4))
#' test <- test[order(ord),]
#' test <- test[order(test$socioeconomic),]
#' head(train)
#' head(test)
#'
#' mod1 <- selr(mntlhlth~socioeconomic,data=train,modelType="continuation")
#' summary(mod1)
#'
#' data<-data.frame(test[,-2])
#' colnames(data)<-c("socioeconomic")
#' pred <- predict(object=mod1,newdata=data)
#' print(pred)
#'
#'
#' @export
predict.selr <- function(object,newdata){
  Z_mat <- object[["Z"]]
  Y <- 1*(Z_mat==1)
  W <- 1*(Z_mat==1|Z_mat==0)
  coefs = coef(object$mod)

  if (missingArg(newdata)){
    yname <- all.vars(formula)[1]
    X_mat = unique(object[["data"]][-yname])
    predicts <- selr_pi_optim(X_mat,W,Y,coefs)
  } else {
    X_mat = unique(newdata)
    predicts <- selr_pi_optim(X_mat,W,Y,coefs)
  }

  predicts <- as.data.frame(matrix(predicts, ncol = nrow(Y), byrow = TRUE))
  colnames(predicts) <- c(rownames(Z_mat))
  predicts = cbind(X_mat, predicts)

  return(predicts)
}

# objective function: negative entropy
eval_f0 <- function(param, N, J) {
  param <- matrix(param[1:(N*J)], nrow=N, ncol=J)

  entropy <- 0

  for (i in 1:N) {
    for (j in 1:J){
      entropy <- entropy + param[i, j]*log(param[i, j])
    }
  }
  return (entropy)
}

# gradient of objective function
eval_grad_f0 <- function(param, N, J) {
  param <- matrix(param[1:(N*J)], nrow=N, ncol=J)

  grad <- numeric(length(param))
  for (i in 1:N) {
    for (j in 1:J){
      grad[(i-1)*J + j] <- (1 + log(param[i, j]))
    }
  }

  return (grad)
}


# constraint function
eval_g0 <- function(param, W, Y, X_mat, alpha, beta) {
  N <- nrow(X_mat)
  J <- nrow(Y)
  param <- matrix(param[1:(N*J)], nrow=N, ncol=J)

  # selr constraints
  pi_constraints <- numeric(N*(J-1))
  i = 1
  n = 1
  while (i <= N*(J-1)) {
    diag_eta <- alpha + as.numeric(as.matrix(X_mat[n, ]) %*% beta)
    Z <- W - Y
    diag_eta_mat <- diag(exp(diag_eta))
    pi_constraints[i:(i+J-2)] <- (t(Y) %*% param[n,] - diag_eta_mat %*% t(Z) %*% param[n,])
    i = i+J-1
    n = n+1
  }

  return (pi_constraints)
}


# gradient of constraint function

eval_jac_g_eq <- function(param, W, Y, X_mat, alpha, beta) {
  N <- nrow(X_mat)
  J <- nrow(Y)
  param <- matrix(param[1:(N*J)], nrow=N, ncol=J)

  pi_constraints <- matrix(0, nrow=N*(J-1), ncol=J*N)

  # selr constraints
  row = 1
  for (n in 1:N) {
    diag_eta <- alpha + as.numeric(as.matrix(X_mat[n, ]) %*% beta)
    exp_diag_eta <- exp(as.vector(diag_eta))

    pi_sum <- rep(0, J-1)
    for (s in 1:(J-1)) {
      cols_Y1_W1 <- which(Y[,s] == 1 & W[,s] == 1)
      cols_Y0_W1 <- which(Y[,s] == 0 & W[,s] == 1)
      pi_sum[s] <- sum(param[n, cols_Y0_W1]) # s in R_s\T_s
    }

    for (s in 1:(J-1)) {
      for (j in 1:J) {
        if ( Y[j,s]==1 & W[j,s]==1){
          pi_constraints[row, ((n-1)*J + j)] <- 1
        } else if ( Y[j,s]==0 & W[j,s]==1){
          pi_constraints[row, ((n-1)*J + j)] <- -exp_diag_eta[s]
        }
      }

      row <- row + 1
    }
  }

  return(pi_constraints)
}


selr_pi_optim <- function(X_mat, W, Y, coefs) {
  N <- nrow(X_mat)
  J <- nrow(W)

  init.pi <- matrix(1/J, nrow=N, ncol=J)

  # use geeglm as initial guess

  alpha = coefs[1:(J-1)]
  beta = coefs[J:length(coefs)]

  res <- nloptr(x0 = init.pi,
                eval_f = function(x0) eval_f0(x0, N, J),
                eval_grad_f = function(x0) eval_grad_f0(x0, N, J),
                eval_g_eq = function(x0) eval_g0(x0, W, Y, X_mat, alpha, beta),
                eval_jac_g_eq = function(x0) eval_jac_g_eq(x0, W, Y, X_mat, alpha, beta),
                lb = c(rep(0, N*J)),
                ub = c(rep(1, N*J)),
                opts = list("algorithm" = "NLOPT_LD_SLSQP",
                            "xtol_rel"=1.0e-12,
                            "ftol_rel" = 1.0e-12,
                            "maxeval" = 10000))

  pi_mat <- matrix(res$solution, nrow=N, ncol=J)

  for (i in 1:N) {
    row_sum <- sum(pi_mat[i, ])
    if (row_sum != 0) {
      pi_mat[i, ] <- pi_mat[i, ] / row_sum
    } else {
      pi_mat[i, ] <- 0
    }
  }


  return(pi_mat)
}



#' @export
QIC.selr <- function(object, ..., tol=.Machine$double.eps) {

  ## cat("QIC\n")
  ## print(env)
  # from QIC.geeglm

  if (! (("selr" %in% class(object)) | ("geeglm" %in% class(object$mod)))) {
    stop("QIC requires a selr or geeglm object as input")
  }

  if (length(list(...))) {
    objects <- lapply(list(object, ...), function(object) object$mod)
  } else {
    object <- object$mod
  }

  # Setup functions
  invert <- ginv

  # Missing:
  # Check correct handling of link and family functions

  # Create function to make the computations
  computeqic <- function(object) {
    ## Fitted and observed values for quasi likelihood
    weights <- object$weights  #### added!
    mu <- object$fitted.values[weights != 0]
    y  <- object$y[weights != 0]

    quasi <- sum(y * log(mu / (1 - mu)) + log(1 - mu)) # only "binomial"

    ## Fit model with independence correlation structure
    object$call$corstr <- "independence"
    object$call$zcor <- NULL

    df <- object$data  ########### added!!
    model.indep <- eval(object$call, envir=df) ######### modified!!
    ## model.indep <- eval(object$call, parent.frame()) ## FIXME parent.frame() is wrong...
    ## model.indep <- update(object, corstr="independence",zcorr=NULL)

    # Trace term (penalty for model complexity)
    AIinverse <- invert(model.indep$geese$vbeta.naiv, tol=tol)
    Vr <- object$geese$vbeta
    trace <- sum(diag(AIinverse %*% Vr))
    params <- length(coef(object)) # Mean parameters in the model

    kpm <- params+length(object$geese$alpha)

    # QIC
    QIC  <- -2*(quasi - trace)
    QICu <- -2*(quasi - params)
    QICC <- QIC + (2 * kpm * (kpm + 1)) / (length(unique(object$id)) - kpm - 1)
    output <- c(QIC, QICu, quasi, trace, params, QICC)
    names(output) <- c("QIC", "QICu", "Quasi Lik", "CIC", "params", "QICC")
    output
  }

  if (length(list(...))) {
    # Make the computations
    results <- lapply(objects, computeqic)

    # Check same data size
    check <- sapply(objects, function(x) {
      length(x$y)
    })

    if (any(check != check[1]))
      warning("models are not all fitted to the same number of observations")

    # Merge the results together in a data.matrix
    res <- do.call("rbind", results)

    # Set the row names corresponding to the models
    Call <- match.call()
    Call$k <- NULL
    row.names(res) <- as.character(Call[-1L])
    res
  } else {
    computeqic(object)
  }
}



#' @title  Squared Absolute RPS
#'
#' @description A modified version of ranked probability score for ordinal
#' classification
#'
#' @name sa_RPS
#'
#'
#' @param p a matrix of predicted probabilities, where each column corresponds
#' to one level
#' @param y a matrix of one-hot representation for observed outcomes
#'
#' @return \code{sa_RPS} returns a vector of squared absolute rps for each
#' observation
#'
#'
#' @author Jiyang Xu \email{jiyangx3@@illinois.edu},
#'         Douglas Simpson \email{dgs@@illinois.edu}
#'
#' @seealso \code{\link[verification]{rps}}
#'
#'
#' @references
#' Galdran, A. (2023). Performance Metrics for Probabilistic Ordinal Classifiers.
#' In: Greenspan, H., et al. Medical Image Computing and Computer Assisted
#' Intervention – MICCAI 2023. MICCAI 2023. Lecture Notes in Computer Science,
#' vol 14222. Springer, Cham. https://doi.org/10.1007/978-3-031-43898-1_35
#'
#' @export
sa_RPS <- function(p, y){
  rps_list <- 0
  K <- ncol(p)
  for (i in 1:nrow(p)){
    p_cumsum <- cumsum(unlist(p[i,]))[-ncol(p)]
    y_cumsum <- cumsum(y[i,])[-ncol(y)]
    rps_list[i] <- (sum(abs(p_cumsum - y_cumsum))/(K-1))^2
  }

  return (rps_list)
}

#' @export
logloss.selr <- function(p, y){
  sum <- 0
  for (i in 1:nrow(p)){
    temp <- 0
    for (j in 1:ncol(p)){
      if (y[i,j]!=0){
        temp <- log(p[i,j])
      }
    }
    sum = sum+temp
  }

  return ( (-2)*sum )
}


