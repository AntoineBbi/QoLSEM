#' Estimation of a QoLSEM model at each visit
#'
#' This function is the generalisation of \code{QoLSEM} function where sereval visits/measurement times are considered
#' in the follow-up of patients (statistic units).
#'
#' @param data (dataframe/matrix) dataset
#' @param col_y (interger) index of the colum associated with the response variable
#' @param col_X1 (vector) index of the colum of variables for X1, by default =0 for without covaraites
#' @param col_X2 (vector) index of the colum of variables for X2, by default =0 for without covaraites
#' @param col_Ty (interger/vector) index of the colum of covariates for y, by default =0 for without covaraites
#' @param col_T1 (interger/vector) index of the colum of covariates for X1, by default =0 for without covaraites
#' @param col_T2 (interger/vector) index of the colum of covariates for X2, by default =0 for without covaraites
#' @param col_visit (interger) index of the colum associated with the visits, by default =0 if only one visit
#'
#' @return A list with the following elements:
#'    \describe{
#'   \item{\code{output.data}}{matrix of data including predicted factors.}
#'   \item{\code{C}}{matrix of the 2 estimated parameters c1 and c2, rows corresponds to the visits/datasets}
#'   \item{\code{A1}}{matrix of parameters associated with the factor in the first variable block}
#'   \item{\code{A2}}{matrix of parameters associated with the factor in the second variable block}
#'   \item{\code{D}}{matrix of parameters associated with the covariates in the structural equation}
#'   \item{\code{D1}}{matrix of parameters associated with the covariates in block 1}
#'   \item{\code{D2}}{matrix of parameters associated with the covariates in block 2}
#'   \item{\code{output.sigma2}}{matrix of variance parameters}
#'    }
#'
#' @author Antoine Barbieri, Myriam Tami
#'
#' @seealso \code{\link{generation.QoLSEM}}
#'
#' @references Barbieri A, Tami M, Bry X, Azria D, Gourgou S, Mollevi C, Lavergne C. (2018) \emph{EM algorithm estimation of a structural equation model for the longitudinal study of the quality of life}. Statistics in Medicine. 37(6) :1031-1046.
#'
#' @examples
#' ## test avec N>1
#' test20 <- generation.QoLSEM(N=20,I=150,
#'                             c=as.matrix(c(2,-2)),
#'                             a1=matrix(c(1:7),nrow=7),
#'                             a2=matrix(c(1:12),nrow=12),
#'                             d=c(80),
#'                             D1=matrix(seq(2,14,2)+70,ncol=7),
#'                             D2=matrix(seq(2,24,2)+20,ncol=12),
#'                             sigma.y=10,sigma.X1=7,sigma.X2=5)
#'
#' ## Estimation step
#' simu20 <- QoLSEM_longitudinal(data=test20,
#'                               col_y=4,
#'                               col_X1=seq(5,5+6),
#'                               col_X2=seq(5+7,5+7+11),
#'                               col_Ty=0,
#'                               col_T1=0,
#'                               col_T2=0,
#'                               col_visit=1)
#'
#' ## Estimation of variance parameters for the 20 visits
#' simu20$output.sigma2
#'
#' ## Estimation of intercept parameters for the 20 visits associated with block 1
#' ## Only the intercept is considered because no covariate is taking into account
#' simu20$output.D1
#'
#' @export
QoLSEM_longitudinal <- function(data,col_y,col_X1,col_X2,col_Ty=0,col_T1=0,col_T2=0,col_visit=0){

  #### intialization
  output.data   <- NULL
  output.C      <- NULL
  output.A1     <- NULL
  output.A2     <- NULL
  output.D      <- NULL
  output.D1     <- NULL
  output.D2     <- NULL
  output.sigma2 <- NULL
  name <- NULL

  if(col_Ty==0){
    Ty <- NULL
  } else {
    Ty <- as.matrix(data[,col_Ty])
  }

  if(col_T1==0){
    T1 <- NULL
  } else {
    T1 <- col_T1
  }

  if(col_T2==0){
    T2 <- NULL
  } else {
    T2 <- col_T2
  }

  if(col_visit!=0){
    visit <- unique(data[,col_visit])
    for(i in 1:length(visit)){
      data_tmp <- na.omit(data[which(data[,col_visit]==visit[i]),])

      SEM_tmp <- QoLSEM(y=as.matrix(data_tmp[,col_y]),
                        X1=as.matrix(data_tmp[,col_X1]),
                        X2=as.matrix(data_tmp[,col_X2]),
                        Ty=as.matrix(data_tmp[,col_Ty]),
                        T1=as.matrix(data_tmp[,col_T1]),
                        T2=as.matrix(data_tmp[,col_T2]),
                        convergence=0.001,nb_it=500)


      output.data <- rbind(output.data,
                           cbind(data_tmp,SEM_tmp$Factors))
      output.C <- rbind(output.C,
                        c(visit[i],SEM_tmp$C))
      output.A1 <- rbind(output.A1,SEM_tmp$A1)
      output.A2 <- rbind(output.A2,SEM_tmp$A2)
      output.D  <- rbind(output.D,SEM_tmp$D)
      output.D1 <- rbind(output.D1,SEM_tmp$D1)
      output.D2 <- rbind(output.D2,SEM_tmp$D2)
      output.sigma2 <- rbind(output.sigma2,SEM_tmp$estimate_sigma2)

      name <- c(name,paste("Visit_", visit[i], sep = ""))

      cat('SEM ', visit[i], '\n',
          'Nb of iterations: ', SEM_tmp$Iteration, '\n',
          'Convergence     : ', SEM_tmp$Convergence, '\n')

    }

  } else {

    SEM_tmp <- SEM_questionnaire_QoL(y=as.matrix(data[,]),
                                     X1=as.matrix(data[,col_X1]),
                                     X2=as.matrix(data[,col_X2]),
                                     Ty=as.matrix(data[,col_Ty]),
                                     T1=as.matrix(data[,col_T1]),
                                     T2=as.matrix(data[,col_T2]),
                                     convergence=0.001,nb_it=500)

    output.data <- rbind(output.data,
                         cbind(data_tmp,SEM_tmp$Factors))
    output.C <- rbind(output.C,
                      c(visit[i],SEM_tmp$C))
    output.A1 <- rbind(output.A1,SEM_tmp$A1)
    output.A2 <- rbind(output.A2,SEM_tmp$A2)
    output.D  <- rbind(output.D,SEM_tmp$D)
    output.D1 <- rbind(output.D1,SEM_tmp$D1)
    output.D2 <- rbind(output.D2,SEM_tmp$D2)
    output.sigma2 <- rbind(output.sigma2,SEM_tmp$estimate_sigma2)

    cat('\n',
        'Nb of iteration:', SEM_tmp$Iteration, '\n',
        'Convergence:', SEM_tmp$Convergence, '\n')

  }

  rownames(output.C) <- name
  rownames(output.A1) <- name
  rownames(output.A2) <- name
  rownames(output.sigma2) <- name

  return(list(output.data=output.data,
              output.C=output.C,
              output.A1=output.A1,
              output.A2=output.A2,
              output.D=output.D,
              output.D1=output.D1,
              output.D2=output.D2,
              output.sigma2=output.sigma2
              )
         )

}
