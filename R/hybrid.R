#' @export
#' @rdname hybrid
hybrid <- function(x,...)
  UseMethod("hybrid")
#' Feature selection with the hybrid algorithm
#'
#' hybrid is an all relevant random forest feature selection wrapper algorithm that uses the corrected impurity importance (Nembrini et al. 2019) as variable importance measure (VIM);
#' Analogously to to Janitza et al. 2018, variables negative impurity importance are interpreted unimportant and used to generate the null distribution and to calculate the p-values.
#' That is, this implementation is fully based on the original implementation of the Vita variable selection algorithm available in \code{\link[ranger]{ranger}} and that of \code{\link[Boruta]{Boruta}}.
#' The method performs a top-down search for relevant features by comparing original attributes' importance with importance achievable at random, estimated using their permuted copies, and progressively eliminating irrelevant features to stabilise that test.
#' @rdname hybrid
#' @method hybrid default
#'
#' @param x data frame of predictors.
#' @param y response vector; factor for classification, numeric vector for regression, \code{Surv} object for survival (supports depends on importance adapter capabilities).
#' @param getImp function used to obtain attribute importance.
#' The default is get_imp_ranger, which runs random forest from the \code{ranger} package and gathers Z-scores of corrected impurities.
#' It should return a numeric vector of a size identical to the number of columns of its first argument, containing importance measure of respective attributes.
#' Any order-preserving transformation of this measure will yield the same result.
#' It is assumed that more important attributes get higher importance. +-Inf are accepted, NaNs and NAs are treated as 0s, with a warning.
#' @param pValue confidence level. Default value should be used.
#' @param mcAdj if set to \code{TRUE}, a multiple comparisons adjustment using the Bonferroni method will be applied. Default value should be used; older (1.x and 2.x) versions of hybrid were effectively using \code{FALSE}.
#' @param maxRuns maximal number of importance source runs.
#' You may increase it to resolve attributes left Tentative.
#' @param holdHistory if set to \code{TRUE}, the full history of importance is stored and returned as the \code{ImpHistory} element of the result.
#' Can be used to decrease a memory footprint of hybrid in case this side data is not used, especially when the number of attributes is huge; yet it disables plotting of such made \code{hybrid} objects and the use of the \code{\link{TentativeRoughFix}} function.
#' @param doTrace verbosity level. 0 means no tracing, 1 means reporting decision about each attribute as soon as it is justified, 2 means the same as 1, plus reporting each importance source run, 3 means the same as 2, plus reporting of hits assigned to yet undecided attributes.
#' @param alpha significance threshold used by Vita
#' @param seed Seed
#' @param ... additional parameters passed to \code{getImp}.
#'
#' @return An object of class \code{hybrid}, which is a list with the following components:
#' \item{finalDecision}{a factor of three value: \code{Confirmed}, \code{Rejected} or \code{Tentative}, containing final result of feature selection.}
#' \item{ImpHistory}{a data frame of importances of attributes gathered in each importance source run.
#' Beside predictors' importances, it contains maximal, mean and minimal importance of shadow attributes in each run.
#' Rejected attributes get \code{-Inf} importance.
#' Set to \code{NULL} if \code{holdHistory} was given \code{FALSE}.}
#' \item{timeTaken}{time taken by the computation.}
#' \item{impSource}{string describing the source of importance, equal to a comment attribute of the \code{getImp} argument.}
#' \item{call}{the original call of the \code{hybrid} function.}
#' @details hybrid iteratively compares importances of attributes with importances of shadow attributes, created by shuffling original ones.
#' Attributes that have significantly worst importance than shadow ones are being consecutively dropped.
#' On the other hand, attributes that are significantly better than shadows are admitted to be Confirmed.
#' Shadows are re-created in each iteration.
#' Algorithm stops when only Confirmed attributes are left, or when it reaches \code{maxRuns} importance source runs.
#' If the second scenario occurs, some attributes may be left without a decision.
#' They are claimed Tentative.
#' You may try to extend \code{maxRuns} or lower \code{pValue} to clarify them, but in some cases their importances do fluctuate too much for hybrid to converge.
#' Instead, you can use \code{\link{TentativeRoughFix}} function, which will perform other, weaker test to make a final decision, or simply treat them as undecided in further analysis.
#' @references
#' Nembrini, S., Koenig, I. R. & Wright, M. N. (2018). The revival of the Gini Importance? Bioinformatics. https://doi.org/10.1093/bioinformatics/bty373.
#' Janitza, S, Celik, E, Boulesteix, AL. (2018). A computationally fast variable importance test for random forests for high-dimensional data. Adv Data Anal Classif.; doi.org: 10.1007/s11634-016-0276-4
#' Kursa, M. B. and Rudnicki, W. R. (2010). Feature Selection with the Boruta Package. Journal of Statistical Software. \emph{Journal of Statistical Software, 36(11)}, p. 1-13. URL: \url{http://www.jstatsoft.org/v36/i11/}.
#'
#' @examples
#' set.seed(777)
#'
#' \dontrun{
#' #hybrid on the "small redundant XOR" problem; read ?srx for details
#' data(srx)
#' hybrid(Y~.,data=srx)->hybrid.srx
#'
#' #Results summary
#' print(hybrid.srx)
#'
#' #Result plot
#' plot(hybrid.srx)
#'
#' #Attribute statistics
#' attStats(hybrid.srx)
#'
#' #Using alternative importance source, rFerns
#' hybrid(Y~.,data=srx,getImp=getImpFerns)->hybrid.srx.ferns
#' print(hybrid.srx.ferns)
#'
#' #Verbose
#' hybrid(Y~.,data=srx,doTrace=2)->hybrid.srx
#'}
#' \dontrun{
#' #hybrid on the iris problem extended with artificial irrelevant features
#' #Generate said features
#' iris.extended<-data.frame(iris,apply(iris[,-5],2,sample))
#' names(iris.extended)[6:9]<-paste("Nonsense",1:4,sep="")
#' #Run hybrid on this data
#' hybrid(Species~.,data=iris.extended,doTrace=2)->hybrid.iris.extended
#' #Nonsense attributes should be rejected
#' print(hybrid.iris.extended)
#' }
#'
#' \dontrun{
#' #hybrid on the HouseVotes84 data from mlbench
#' library(mlbench); data(HouseVotes84)
#' na.omit(HouseVotes84)->hvo
#' #Takes some time, so be patient
#' hybrid(Class~.,data=hvo,doTrace=2)->Bor.hvo
#' print(Bor.hvo)
#' plot(Bor.hvo)
#' plotImpHistory(Bor.hvo)
#' }
#' \dontrun{
#' #hybrid on the Ozone data from mlbench
#' library(mlbench); data(Ozone)
#' library(randomForest)
#' na.omit(Ozone)->ozo
#' hybrid(V4~.,data=ozo,doTrace=2)->Bor.ozo
#' cat('Random forest run on all attributes:\n')
#' print(randomForest(V4~.,data=ozo))
#' cat('Random forest run only on confirmed attributes:\n')
#' print(randomForest(ozo[,getSelectedAttributes(Bor.ozo)],ozo$V4))
#' }
#' \dontrun{
#' #hybrid on the Sonar data from mlbench
#' library(mlbench); data(Sonar)
#' #Takes some time, so be patient
#' hybrid(Class~.,data=Sonar,doTrace=2)->Bor.son
#' print(Bor.son)
#' #Shows important bands
#' plot(Bor.son,sort=FALSE)
#' }
#' @export
hybrid.default <- function(x, y,
                           pValue = 0.01,
                           mcAdj = TRUE,
                           maxRuns = 100,
                           doTrace = 0,
                           holdHistory = TRUE,
                           getImp,
                           alpha = 0.05,
                           seed, ...){
  #Timer starts... now!
  timeStart<-Sys.time()

  #Extract the call to store in output
  cl<-match.call()
  cl[[1]]<-as.name('hybrid')

  #Convert x into a data.frame
  if(!is.data.frame(x))
    x<-data.frame(x)

  ##Some checks on x & y
  # if(length(grep('^shadow',names(x)))>0)
  #   stop('Attributes with names starting from "shadow" are reserved for internal use. Please rename them.')
  if(any(c(is.na(x),is.na(y))))
    stop('Cannot process NAs in input. Please remove them.')
  if(maxRuns<11)
    stop('maxRuns must be greater than 10.')

  ##Expands the information system with newly built random attributes and calculates importance
  addShadowsAndGetImp<-function(decReg,runs){
    #Notifying user of our progress
    if(doTrace>1)
      message(sprintf(' %s. run of importance source...',runs))

    #Calling importance source; "..." can be used by the user to pass rf attributes (for instance ntree)
    # impRaw<-getImp(cbind(x[,decReg!="Rejected"],xSha),y,...)
    impRaw <- getImp(cbind(x[,decReg!="Rejected"]),y, seed = seed + runs, ...)
    if(!is.numeric(impRaw))
      stop("getImp result is not a numeric vector. Please check the given getImp function.")
    # if(length(impRaw)!=sum(decReg!="Rejected")+ncol(xSha))
    if(length(impRaw)!=sum(decReg!="Rejected"))
      stop("getImp result has a wrong length. Please check the given getImp function.")
    if(any(is.na(impRaw)|is.nan(impRaw))){
      impRaw[is.na(impRaw)|is.nan(impRaw)]<-0
      warning("getImp result contains NA(s) or NaN(s); replacing with 0(s), yet this is suspicious.")
    }

    #Importance must have Rejected attributes put on place and filled with -Infs
    # imp<-rep(-Inf,nAtt+nSha);names(imp)<-c(attNames,names(xSha))
    imp<-rep(-Inf,nAtt);names(imp)<-c(attNames)
    # impRaw->imp[c(decReg!="Rejected",rep(TRUE,nSha))]
    impRaw->imp[c(decReg!="Rejected")]
    # shaImp<-imp[(nAtt+1):length(imp)];imp[1:nAtt]->imp
    if((!any(imp < 0)) & (runs == 1)){
      warning("No negative importance. Shadow variables may be required.")
    }
    shaImp<- c(0, abs(imp[imp <= 0])); imp[1:nAtt]->imp
    return(list(imp=imp,shaImp=shaImp))
  }

  ##Assigns hits
  assignHits<-function(hitReg,curImp){
    curImp$imp>=quantile(curImp$shaImp, 1 - alpha)->hits
    if(doTrace>2){
      uncMask<-decReg=="Tentative"
      intHits<-sum(hits[uncMask])
      if(intHits>0)
        message(sprintf("Assigned hit to %s attribute%s out of %s undecided.",sum(hits[uncMask]),if(intHits==1) "" else "s",sum(uncMask)))
      else
        message("None of undecided attributes scored a hit.")
    }
    hitReg[hits]<-hitReg[hits]+1
    return(hitReg)
  }

  ##Checks whether number of hits is significant
  doTests<-function(decReg,hitReg,runs){
    pAdjMethod<-ifelse(mcAdj[1],'bonferroni','none')
    #If attribute is significantly more frequent better than shadowMax, its claimed Confirmed
    toAccept<-stats::p.adjust(stats::pbinom(hitReg-1,runs,0.5,lower.tail=FALSE),method=pAdjMethod)<pValue
    (decReg=="Tentative" & toAccept)->toAccept

    #If attribute is significantly more frequent worse than shadowMax, its claimed Rejected (=irrelevant)
    toReject<-stats::p.adjust(stats::pbinom(hitReg,runs,0.5,lower.tail=TRUE),method=pAdjMethod)<pValue
    (decReg=="Tentative" & toReject)->toReject

    #Update decReg
    decReg[toAccept]<-"Confirmed";"Rejected"->decReg[toReject]

    #Report progress
    if(doTrace>0){
      nAcc<-sum(toAccept)
      nRej<-sum(toReject)
      nLeft<-sum(decReg=="Tentative")
      if(nAcc+nRej>0)
        message(sprintf("After %s iterations, +%s: ",runs,format(difftime(Sys.time(),timeStart),digits=2)))
      if(nAcc>0)
        message(sprintf(" confirmed %s attribute%s: %s",
                        nAcc,ifelse(nAcc==1,'','s'),.attListPrettyPrint(attNames[toAccept])))
      if(nRej>0)
        message(sprintf(" rejected %s attribute%s: %s",
                        nRej,ifelse(nRej==1,'','s'),.attListPrettyPrint(attNames[toReject])))
      if(nAcc+nRej>0)
        if(nLeft>0){
          message(sprintf(" still have %s attribute%s left.\n",
                          nLeft,ifelse(nLeft==1,'','s')))
        }else{
          if(nAcc+nRej>0) message(" no more attributes left.\n")
        }
    }
    return(decReg)
  }

  ##Creating some useful constants
  nAtt<-ncol(x); nrow(x)->nObjects
  attNames<-names(x); c("Tentative","Confirmed","Rejected") -> confLevels

  ##Initiate state
  decReg<-factor(rep("Tentative",nAtt),levels=confLevels)
  hitReg<-rep(0,nAtt);names(hitReg)<-attNames
  impHistory<-list()
  runs<-0

  ##Main loop

  while(any(decReg=="Tentative") && (runs+1->runs)<maxRuns){
    curImp<-addShadowsAndGetImp(decReg,runs)
    hitReg<-assignHits(hitReg,curImp)
    decReg<-doTests(decReg,hitReg,runs)

    #If needed, update impHistory with scores obtained in this iteration
    if(holdHistory){
      imp<-c(curImp$imp,
             shadowMax=max(curImp$shaImp),
             shadowMean=mean(curImp$shaImp),
             shadowMin=min(curImp$shaImp))
      impHistory<-c(impHistory,list(imp))
    }
  }

  ##Building result
  impHistory<-do.call(rbind,impHistory)
  names(decReg)<-attNames
  ans<-list(finalDecision=decReg,ImpHistory=impHistory,
            pValue=pValue,maxRuns=maxRuns,light=TRUE,mcAdj=mcAdj,
            timeTaken=Sys.time()-timeStart,roughfixed=FALSE,call=cl,
            impSource=comment(getImp))

  "hybrid" -> class(ans)
  return(ans)
}

.attListPrettyPrint<-function(x,limit=5){
  x<-sort(x)
  if(length(x)<limit+1)
    return(sprintf("%s;",paste(x,collapse=", ")))
  sprintf("%s and %s more;",paste(utils::head(x,limit),collapse=", "),length(x)-limit)
}

#' @rdname hybrid
#' @method hybrid formula
#' @param formula alternatively, formula describing model to be analysed.
#' @param data in which to interpret formula.
#' @export
hybrid.formula<-function(formula,data=.GlobalEnv,...){
  ##Grab and interpret the formula
  stats::terms.formula(formula,data=data)->t
  x<-eval(attr(t,"variables"),data)
  apply(attr(t,"factors"),1,sum)>0->sel
  nam<-rownames(attr(t,"factors"))[sel]
  data.frame(x[sel])->df;names(df)<-nam
  x[[attr(t,"response")]]->dec

  ##Run hybrid
  ans<-hybrid.default(df,dec,...)
  ans$call<-match.call()
  ans$call[[1]]<-as.name('hybrid')
  formula->ans$call[["formula"]]
  return(ans)
}

#' Print hybrid object
#'
#' Print method for the hybrid objects.
#' @method print hybrid
#' @param x an object of a class hybrid.
#' @param ... additional arguments passed to \code{\link{print}}.
#' @return Invisible copy of \code{x}.
#' @export
print.hybrid <- function(x,...){
  if(class(x)!='hybrid') stop("This is NOT a hybrid object!")
  cat(paste('hybrid performed ',dim(x$ImpHistory)[1],' iterations in ',format(x$timeTaken),'.\n',sep=''))
  if(x$roughfixed) cat(paste('Tentatives roughfixed over the last ',x$averageOver,' iterations.\n',sep=''))
  if(sum(x$finalDecision=='Confirmed')==0){
    cat(' No attributes deemed important.\n')} else {
      writeLines(strwrap(paste(sum(x$finalDecision=='Confirmed'),' attributes confirmed important: ',
                               .attListPrettyPrint(names(x$finalDecision[x$finalDecision=='Confirmed']))),indent=1))
    }
  if(sum(x$finalDecision=='Rejected')==0){
    cat(' No attributes deemed unimportant.\n')} else {
      writeLines(strwrap(paste(sum(x$finalDecision=='Rejected'),' attributes confirmed unimportant: ',
                               .attListPrettyPrint(names(x$finalDecision[x$finalDecision=='Rejected']))),indent=1))
    }
  if(sum(x$finalDecision=='Tentative')!=0){
    writeLines(strwrap(paste(sum(x$finalDecision=='Tentative'),' tentative attributes left: ',
                             .attListPrettyPrint(names(x$finalDecision[x$finalDecision=='Tentative']))),indent=1))
  }
  invisible(x)
}
