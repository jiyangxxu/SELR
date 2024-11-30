



#'Mental Health Status
#'
#' @format The \code{mntlhlth} data frame has 1660 rows and 2 columns.
#'
#' @details Data of a study on the relationship between an individual’s mental health status
#' and the socioeconomic status of his or her parents. 'socioeconomic' has 0 to 5 levels,
#' from 'lowest' to 'highest' and 'mntlhlth' has four levels: 'well', 'mild', 'moderate', 'impaired'
#'
#'
#'
#' \describe{
#' \item{socioeconomic}{The socioeconomic status of the individual's parents (0~5)}
#' \item{mntlhlth}{an individual’s mental health status ('well', 'mild', 'moderate', 'impaired')}
#' }
#'
#' @references
#' Fu, L., & Simpson, D. G. (2002).
#' Conditional risk models for ordinal response data:
#' simultaneous logistic regression analysis and generalized score tests.
#' Journal of Statistical Planning and Inference, 108(1–2), 201–217.
#' https://doi.org/10.1016/s0378-3758(02)00279-3
#'
#'
#' @keywords datasets
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
#' anova(mod1,mod2,test="score")
#' anova(mod1,test="score")
#'
#' QIC(mod1)
#' QIC(mod1,mod2)
#'
#' mod3 <- selr(mntlhlth~socioeconomic,data=mntlhlth,modelType="cumulative",parallel=FALSE,elr=TRUE)
#' summary(mod3)
#'
"mntlhlth"







#' White Wine Bitterness
#'
#' @format The \code{bitterness} data frame has 72 rows and 3 columns.
#'
#' @details Data of a study on the factors influencing white wine bitterness. Bitterness
#' is measured by 1~5 level from "least bitter" to "most bitter". This data frame
#' has already been reformatted for SELR package usage.
#'
#'
#'
#' \describe{
#' \item{temperature}{The temperature controlled during crushing the grapes ('warm' or 'cold')}
#' \item{contact}{The contact controlled during crushing the grapes ('yes' or 'no')}
#' \item{bitterness}{The bitterness of the white wine (1~5)}
#' }
#'
#' @references
#' LRandall, J. H. (1989). The Analysis of Sensory Data by Generalized Linear Model.
#' Biometrical Journal, 31(7), 781–793. https://doi.org/10.1002/bimj.4710310703
#'
#' @keywords datasets
#' @examples
#'
#' data(bitterness)
#' head(bitterness)
#'
#' mod1 <- selr(bitterness~temperature*contact,data=bitterness,modelType="cumulative")
#' summary(mod1)
#'
#' mod2 <- selr(bitterness~temperature+contact,data=bitterness,modelType="cumulative",parallel = FALSE)
#' summary(mod2)
#'
#' mod3 <- selr(bitterness~temperature+contact,data=bitterness,modelType="continuation",elr=TRUE)
#' summary(mod3)
#'
#' anova(mod1,test="score")
#'
#' QIC(mod1)
#' QIC(mod1,mod2)
#'
"bitterness"




#' Preterm Birth
#'
#' @format The \code{ptb} data frame has 424 rows and 17 columns.
#'
#' @details Data of a study on the factors influencing preterm birth. The response
#' are nominal factors: "FTB", "sPTB", "mPTB", where "FTB" for Full-Term Birth,
#' "sPTB" for Spontaneous Preterm Birth, "mPTB" for Moderate Preterm Birth.
#'
#'
#' @references
#'
#'
#' @keywords datasets
#' @examples
#'
#' data(ptb)
#' head(ptb)
#'
"ptb"





