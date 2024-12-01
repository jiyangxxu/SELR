% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selr.R
\name{selr}
\alias{selr}
\title{Simultaneous and Efficient Logistic Regression (SELR)}
\usage{
selr(
  formula,
  data,
  subset = NULL,
  na.action,
  offset,
  contrasts = NULL,
  link = c("logit", "identity", "probit", "cloglog"),
  modelType = c("nonstandard", "cumulative", "adjacent", "continuation", "baseline",
    "sequential", "twolevel", "binarysplit"),
  risk = NULL,
  target = NULL,
  Z = NULL,
  idmod = NULL,
  parallel = T,
  reverse = F,
  ordinal = T,
  yord = NULL,
  elr = F,
  glmnet = NULL,
  ...
)
}
\arguments{
\item{formula}{See corresponding documentation to \code{glm}}

\item{data}{See corresponding documentation to \code{glm}}

\item{subset}{See corresponding documentation to \code{glm}}

\item{na.action}{See corresponding documentation to \code{glm}}

\item{offset}{See corresponding documentation to \code{glm}}

\item{contrasts}{See corresponding documentation to \code{glm}}

\item{link}{See corresponding documentation to \code{geeglm}.
The following are allowed: \code{'identity'}, \code{'logit'},
\code{'probit'}, \code{'cloglog'}. The default link is \code{'logit'}.}

\item{modelType}{a character string specifying the ratio value.
The following are permitted: \code{'cumulative'},
\code{'adjacent'}, \code{'continuation'}, \code{'sequential'},
\code{'baseline'}, \code{'twolevel'}, \code{'binarysplit'},
\code{'nonstandard'}. Note that if the ratio is \code{'nonstandard'},
specify it through \code{'risk'} and \code{'target'} or 'Z' below.
'risk' and 'target' or 'Z' should not be assigned if ratio is not
\code{'nonstandard'}. Default value is \code{'nonstandard'}.}

\item{risk}{a user defined \eqn{S*(S-1)} matrix, also
known as 'risk set'.}

\item{target}{a user defined \eqn{S*(S-1)} matrix, also
known as 'target set'.}

\item{Z}{a user defined \eqn{S*(S-1)} matrix. 'Z' is using one
matrix to represent both \code{'risk'} and \code{'target'}. Please note
that 'Z' should not be assigned if \code{'risk'} and \code{'target'} are assigned.}

\item{idmod}{a string of column names of \code{data} clarifying the clustering}

\item{parallel}{a logical variable. If parallel is true,
selr will compute parallelism case for \code{'levmod'}. Otherwise,
selr will compute nonparallelism case. Default is true.}

\item{reverse}{a logical variable. Please note that reverse
only apply to ratios of \code{'cumulative'}, \code{'adjacent'},
\code{'continuation'}, \code{'sequential'}, and \code{'baseline'}.
Default is false.}

\item{ordinal}{a logical variable, indicating whether fitted data is ordinal
or not. Default is true.}

\item{yord}{a list contained the desired order of response levels used for
comparison from \code{risk} and \code{target}. Please note that
this argument is not necessary to be assigned if ordinal is true.}

\item{elr}{a logical variable, indicating whether using efficient logistics
regression or not. Default is false, which is using the simultaneous
logistics regression.}

\item{glmnet}{a \code{glmnet} object. If not null, the model would return an
additional \code{glmnet} model. See more from \code{glmnet}.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
\code{selr} returns a list contains three values: a selr object,
    the ratio, the value of Z matrix, the original input data set, and
    a \code{glmnet} object if the user requests a \code{glmnet} call
}
\description{
The selr function fits simultaneous and efficient logistic
     regressions using the 'selr' function of the 'selr' package.
     The underlying function is 'geeglm' in 'Geepack' library.
     anova "Wald" and "score" are available.
}
\note{
Default models are \eqn{\text{link}{\Pr(Y_i\geq s)}}, \eqn{\text{link}{\Pr(Y_i=s|Y_i=s-1)}},
   \eqn{\text{link}{\Pr(Y_i=s+1|Y_i\leq s+1)}}, \eqn{\text{link}{\Pr(Y_i\leq s|Y_i\leq s+1)}}
   for 'cumulative', 'adjacent', 'continuation', 'sequential' ratios.
   If reverse is true, models are \eqn{\text{link}{\Pr(Y_i\geq s)}},
   \eqn{\text{link}{\Pr(Y_i=s|Y_i=s+1}}, \eqn{\text{link}{\Pr(Y_i=s|Y_i\geq s)}},
   \eqn{\text{link}{\Pr(Y_i\geq s|Y_i\geq s-1)}} with the same order.
}
\section{Warning}{
 \code{'elr'} should be used with more care. If the response
    is not numeric class but ordinal, please sort the data set based desired
    order of response before using \code{'selr'}.
}

\examples{

data(mntlhlth)
head(mntlhlth)

mod1 <- selr(mntlhlth~socioeconomic,data=mntlhlth,modelType="cumulative",parallel=FALSE)
summary(mod1)

mod2 <- selr(mntlhlth~socioeconomic,data=mntlhlth,modelType="cumulative")
summary(mod2)

anova(mod1,test="score")
anova(mod1,mod2,test="score")

QIC(mod1)
QIC(mod1,mod2)

mod3 <- selr(mntlhlth~socioeconomic,data=mntlhlth,modelType="continuation",parallel=FALSE,elr=TRUE)
summary(mod3)

mod4 <- selr(mntlhlth~socioeconomic,data=mntlhlth,modelType="continuation",elr=TRUE)
summary(mod4)

anova(mod4,test="score")
anova(mod3,mod4,test="Wald")

QIC(mod3)
QIC(mod3,mod4)

# non-standard 1 (Y,W)
# {1} vs. {2}, {1,2} vs. {3,4}, {3} vs. {4}
W <- matrix(c(1,1,0,0,1,1,1,1,0,0,1,1), nrow=4, ncol=3)
W

Y <- matrix(c(0,1,0,0,0,0,1,1,0,0,0,1), nrow=4, ncol=3)
Y

Z <- matrix(c(0,1,-1,-1,0,0,1,1,-1,-1,0,1), nrow=4, ncol=3)
Z

mod5 <- selr(mntlhlth~socioeconomic,data=mntlhlth,risk=W,target=Y)
summary(mod5)

mod6 <- selr(mntlhlth~socioeconomic,data=mntlhlth,Z=Z)
summary(mod6)


}
\references{
Fu, L., & Simpson, D. G. (2002).
Conditional risk models for ordinal response data:
simultaneous logistic regression analysis and generalized score tests.
Journal of Statistical Planning and Inference, 108(1–2), 201–217.
https://doi.org/10.1016/s0378-3758(02)00279-3

Halekoh, U., Højsgaard, S., & Yan, J. (2006).
The R Package geepack for Generalized Estimating Equations.
Journal of Statistical Software, 15(2). https://doi.org/10.18637/jss.v015.i02

Ekstrøm, C. T. (2022, June 20).
Miscellaneous Esoteric Statistical Scripts [R package MESS version 0.5.9].
https://cran.r-project.org/web/packages/MESS
}
\seealso{
\code{\link{geeglm}}, \code{\link{glm}}, \code{\link{anova.glm}}
}
\author{
Jiyang Xu \email{jiyangx3@illinois.edu},
        Douglas Simpson \email{dgs@illinois.edu}
}
\keyword{models}
