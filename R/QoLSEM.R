#' Quality of Life analysis using a structural equation model (QoLSEM)
#'
#' This function fits a structural equation model for the analysis of the Quality of life
#' data from EORTC questionnaires. The estimation is achieved through the EM algorithm.
#' More details can be found in the article cited reference section.
#'
#' @param y a numeriacl vector of observations of the response variable
#' @param X1 matrix of a first block of data associated with the first factor
#' @param X2 matrix of a second block of data associated with the second factor
#' @param Ty matrix/vecteur of covariable(s) associated with \code{Y}
#' @param T1 matrix/vecteur covariable(s) associated with \code{X1}
#' @param T2 matrix/vecteur covariable(s) associated with \code{X2}
#' @param convergence a numerial scalar denoting the criteria of convergence to stop the iteration; =0.001 by default
#' @param nb_it a numerical scalar denoting the maximal number of iterations of the EM algorithm; =500 by default
#' @param trace if TRUE, the function shows for each iteration the convergence of algorithm (default=FALSE)
#'
#' @return A list with the following elements:
#'    \describe{
#'   \item{\code{Factors}}{matrix of size (nx2) of the predicted factors. The first (reps. second) column is associated with the first (resp. second) factor}
#'   \item{\code{C}}{Vector of 2 scalars corresponding to the estimated parameters c1 and c2}
#'   \item{\code{A1}}{vector of parameters associated with the factor in the first variable block}
#'   \item{\code{A2}}{vector of parameters associated with the factor in the second variable block}
#'   \item{\code{D}}{vector of parameters associated with the covariates in the structural equation}
#'   \item{\code{D1}}{matrix of parameters associated with the covariates in block 1}
#'   \item{\code{D2}}{matrix of parameters associated with the covariates in block 2}
#'   \item{\code{estimate_sigma2}}{vector of variance parameters}
#'   \item{\code{Convergence}}{argument \code{diff}}
#'   \item{\code{Iteration}}{Number of iterations until the convergence}
#'    }
#'
#' @import expm graphics FactoMineR MASS
#'
#' @seealso \code{\link{generation.QoLSEM}}
#'
#' @author Antoine Barbieri, Myriam Tami
#'
#' @references Barbieri A, Tami M, Bry X, Azria D, Gourgou S, Mollevi C, Lavergne C. (2018) \emph{EM algorithm estimation of a structural equation model for the longitudinal study of the quality of life}. Statistics in Medicine. 37(6):1031-1046.
#'
#' @examples
#' ## generation of a dataset
#' test <- generation.QoLSEM(N=1,
#'                           I=150,
#'                           c=as.matrix(c(2,-2)),
#'                           a1=matrix(c(1:7),nrow=7),
#'                           a2=matrix(c(1:12),nrow=12),
#'                           d=c(80),
#'                           D1=matrix(seq(2,14,2)+70,ncol=7),
#'                           D2=matrix(seq(2,24,2)+20,ncol=12),
#'                           sigma.y=10,
#'                           sigma.X1=10,
#'                           sigma.X2=10)
#'
#' ## Estimation of the model
#' simu <- QoLSEM(test$y,
#'                test$X1,
#'                test$X2,
#'                Ty=NULL,
#'                T1=NULL,
#'                T2=NULL,
#'                convergence=0.001,
#'                nb_it=500,
#'                trace=FALSE)
#'
#' ## Predicted subject-speficic factors
#' simu$Factors
#'
#' ## Estimation parameters associated with the intercept of the block 1
#' simu$D1
#'
#' @export
QoLSEM <- function(y, X1, X2, Ty=NULL, T1=NULL, T2=NULL, convergence=0.001, nb_it=500, trace=FALSE)
  {

  ##### used object
  qy  <- ncol(as.matrix(y))
  qX1 <- ncol(as.matrix(X1))
  qX2 <- ncol(as.matrix(X2))
  n   <- nrow(as.matrix(y))
  Ty  <- cbind(rep(1,n),Ty)
  T1  <- cbind(rep(1,n),T1)
  T2  <- cbind(rep(1,n),T2)
  qT  <- ncol(Ty)
  qT1 <- ncol(T1)
  qT2 <- ncol(T2)
  if(n!=nrow(X1) | n!=nrow(X2)) {
    cat("Warning, row number is different between Y and X1 or/and X2 !")
  }


  #######################################################################################
  #### Initialisation of parameters through ACP and regression  #########################

  #- ACP to initialize the factors F1 and F2
  #- Regression to initialize the parameters

  #### Concerning X1
  # give the fist composant from FactomineR package (factor F1)
  tmp1 <- PCA(X1, scale.unit=FALSE,graph = FALSE)$ind$coord[,1]
  tmp1 <- cr(tmp1)

  # Initiation of D1 (regression coefficients associated with T1)
  # Remarque: D1[1] is the intercept (or the "mean")
  # estimation of parameter muX1 and vector of parameter A1_sim
  regX1_Fact <-lm(as.matrix(X1)~T1+tmp1-1)
  D1 <- coef(regX1_Fact)[1:qT1,]
  A1 <- coef(regX1_Fact)[(qT1+1),]

  # initialisation of sigma2_X1 parameter
  sig_lm_qX1_Fact <- rep(NA,qX1)
  for(i in 1:qX1){
    sig_lm_qX1_Fact[i] <- summary(lm(as.matrix(X1)~T1+tmp1-1))[[i]]$sigma
  }
  sigma2_X1 <- mean(sig_lm_qX1_Fact^2)


  #### Concerning X2
  # give the fist composant from FactomineR package (factor F2)
  tmp2 <- PCA(X2, scale.unit=FALSE,graph = FALSE)$ind$coord[,1]
  tmp2 <- cr(tmp2)

  # Initiation of D2 (regression coefficients associated with T2)
  # Remarque: D2[1] is the intercept (or the "mean")
  # estimation of parameter muX1 and vector of parameter A2
  regX2_Fact <-lm(as.matrix(X2)~T2+tmp2-1)
  D2 <- coef(regX2_Fact)[1:qT2,]
  A2 <- coef(regX2_Fact)[(qT2+1),]

  # initialisation of sigma2X parameter
  sig_lm_qX2_Fact <- rep(NA,qX2)
  for(i in 1:qX2){
    sig_lm_qX2_Fact[i] <- summary(lm(as.matrix(X2)~T2+tmp2-1))[[i]]$sigma
  }
  sigma2_X2 <- mean(sig_lm_qX2_Fact^2)

  #### Concerning y
  model <- lm(y~Ty + tmp1 + tmp2-1)
  d <- model$coefficients[1:qT]
  C1 <- model$coefficients[qT+1]
  C2 <- model$coefficients[qT+2]
  sigma2_y <- (summary(model)$sigma)^2

  ## needed to calculate: M-i and Sigma_i for distribution
  #which gives F_tilde, G_tilde,Phi_tilde et Gamma_tilde
  Psi_X1 <- diag(sigma2_X1, qX1)
  Psi_X2 <- diag(sigma2_X2, qX2)
  #Psi_Y <- diag(sigma2_y, qy)

  #### initialization concerning iterations
  it <- 0
  diff <- 1
  diffgraph <- NULL
  detE2 <- NULL

  while( (diff > convergence) && (it < nb_it) )  {
    # using &&: only the first argument is assessed. And if TRUE, the second is check

    d_it <- d
    D1_it<- D1
    D2_it<- D2
    A1_it <- A1
    A2_it <- A2
    C1_it <- C1
    C2_it <- C2
    sigma2_y_it <- sigma2_y
    sigma2_X1_it <- sigma2_X1
    sigma2_X2_it <- sigma2_X2

    #Calcule pour les n individus i de : phi_tilde, gama_tilde, f_tilde et g_tilde :
    #Besoin dabord de difinir M_i et GAMMA_i resp. Le vecteur moy et la matrice de var-cov de la loi normale de H_i|Z_i

    # size: (2x[1+qX1+Qx2])
    E1 <- as.matrix(cbind(c(C1,C2),
                          rbind(A1,rep(0,qX1)),
                          rbind(rep(0,qX2),A2)
    )
    )

    # size: ([1+qX1+Qx2]x[1+qX1+Qx2])
    E2 <- as.matrix(rbind(c(C1^2+C2^2+sigma2_y,C1*A1,C2*A2),
                          cbind(C1*A1,A1%*%t(A1)+Psi_X1,matrix(0,nrow=qX1,ncol=qX2)),
                          cbind(C2*A2,matrix(0,nrow=qX2,ncol=qX1),A2%*%t(A2)+Psi_X2)
    )
    )
    detE2 <- c(detE2,det(E2))

    # VOIR avec les listes (plus facile ?!)
    # size: ( [1+qX1+Qx2] x 1 x n )
    #     E3 <- array( 1:((qy+qX1+qX2)*1*n), dim = c(qy+qX1+qX2, 1, n))
    #     for(i in 1:n){
    #       E3[,,i] = t(c( y[i] - t(as.matrix(d)) %*% Ty[i,],
    #                      X1[i,] - t(as.matrix(D1) %*% as.matrix(T1[i,])),
    #                      X2[i,] - t(as.matrix(D2) %*% as.matrix(T2[i,]))
    #                      )
    #                   )
    #     }
    #     E3_bis <- matrix( rep(1,((qy+qX1+qX2)*n)), nrow=n, ncol=(qy+qX1+qX2))
    #     for(i in 1:n){
    #       E3_bis[i,] = t(c( y[i] - t(d) %*% Ty[i,],
    #                         X1[i,] - T1[i,] %*% D1,
    #                         X2[i,] - T2[i,] %*% D2
    #                         )
    #                      )
    #     }

    E3_y  <- y- Ty %*% d
    E3_X1 <- X1 - T1 %*% D1
    E3_X2 <- X2 - T2 %*% D2
    E3 <- cbind(E3_y,E3_X1,E3_X2)

    # size: (2x2)
    E4 <- diag(1,2)

    # E12: working matrix
    # GAMMA size: (2x2)
    E12 <- E1%*%ginv(E2)
    GAMMA <- E4-E12%*%t(E1)

    # M_i is E1%*%ginv(E2)%*%E3[,,i] which is E12%*%E3[,,i]
    #factor <- prod_tab_mat(E12 , E3)
    # factor is an array of size [1+Qx1+Qx2,1,n] de i=1 ` n iliments factor[,,i]
    # factor[,1,i] is the vector of M_i concerning the individual i
    factor <- E12%*%t(E3)


    ############# Step E: expectation values

    F1_tilde <- factor[1,]
    F2_tilde <- factor[2,]
    phi11 <- t(F1_tilde)%*%F1_tilde+n*GAMMA[1,1]
    phi22 <- t(F2_tilde)%*%F2_tilde+n*GAMMA[2,2]
    phi12 <- t(F1_tilde)%*%F2_tilde+n*GAMMA[1,2]
    phi21 <- t(F2_tilde)%*%F1_tilde+n*GAMMA[2,1]
    phi1_tilde <- F1_tilde^2+GAMMA[1,1]
    phi2_tilde <- F2_tilde^2+GAMMA[2,2]

    ############# Step M: Value of parameters

    #     F1_tilde_bar <- mean(F1_tilde)
    #     F2_tilde_bar <- mean(F2_tilde)
    #     phi1_tilde_bar <- mean(phi1_tilde)
    #     phi2_tilde_bar <- mean(phi2_tilde)
    #
    #     F1_tilde_T1_bar <- apply(F1_tilde*T1,2,mean)
    #     inv_T1_prim_T1_bar <- solve(t(T1)%*%T1/n)
    #     denom_A1 <- mean(phi1_tilde)-t(F1_tilde_T1_bar)%*%inv_T1_prim_T1_bar%*%F1_tilde_T1_bar
    #
    #     F2_tilde_T2_bar <- apply(F2_tilde*T2,2,mean)
    #     inv_T2_prim_T2_bar <- solve(t(T2)%*%T2/n)
    #     denom_A2 <- mean(phi2_tilde)-t(F2_tilde_T2_bar)%*%inv_T2_prim_T2_bar%*%F2_tilde_T2_bar
    #
    #     X1_bar <- apply(X1, MARGIN=2, mean)
    #     X2_bar <- apply(X2, MARGIN=2, mean)
    #     y_bar  <- mean(y)

#     M_T <- solve(t(T1)%*%T1)
#     X1_T1_prim <- t(X1)%*%T1
#     F1_tilde_T1 <- apply(F1_tilde*T1,2,sum)
#     A1 <- ((t(X1)%*%F1_tilde-X1_T1_prim%*%M_T%*%F1_tilde_T1)
#            %*%solve(phi11-t(F1_tilde_T1)%*%M_T%*%F1_tilde_T1))
#
#
#         X2_T2_prim_bar <- t(X2)%*%T2/n
#         A2 <- t(apply(F2_tilde*X2,2,mean)
#                 -X2_T2_prim_bar%*%inv_T2_prim_T2_bar%*%F2_tilde_T2_bar
#         )/(as.numeric(denom_A2))



    # pour gagner du temps de calcul (- au + rapide) :
    # (matrice1x1) t(T1)%*%T1 == sum(T1^2) (numeric)

    T1T1 <- solve(t(T1)%*%T1)
    tmp_A1 <- t(F1_tilde)%*%T1%*%T1T1%*%t(T1)
    A1 <- (solve(phi11-tmp_A1%*%F1_tilde)%*%(t(F1_tilde)%*%X1-tmp_A1%*%X1))

    T2T2 <- solve(t(T2)%*%T2)
    tmp_A2 <- t(F2_tilde)%*%T2%*%T2T2%*%t(T2)
    A2 <- (solve(phi22-tmp_A2%*%F2_tilde)%*%(t(F2_tilde)%*%X2-tmp_A2%*%X2))

    D1 <- (T1T1)%*%(t(T1)%*%X1-t(T1)%*%F1_tilde%*%A1)
    D2 <- (T2T2)%*%(t(T2)%*%X2-t(T2)%*%F2_tilde%*%A2)

    tmp  <- Ty%*%solve(t(Ty)%*%Ty)%*%t(Ty)
    tmp1 <- t(F1_tilde)%*%tmp
    tmp2 <- solve(tmp1%*%F1_tilde+phi11)
    tmp3 <- t(F2_tilde)%*%tmp
    tmp4 <- tmp2%*%t(F1_tilde)%*%y-tmp2%*%tmp1%*%y
    C2 <- (solve(tmp3%*%F2_tilde
                 +phi22
                 -tmp3%*%F1_tilde%*%tmp2%*%tmp1%*%F2_tilde
                 +tmp3%*%F1_tilde%*%tmp2%*%phi12
                 -phi21%*%tmp2%*%tmp1%*%F2_tilde
                 -phi21%*%tmp2%*%phi12)
           %*%(t(F2_tilde)%*%y-tmp3%*%y-tmp3%*%F1_tilde%*%tmp4-phi21%*%tmp4)
    )

    C1 <- (tmp2%*%t(F1_tilde)%*%y
           -tmp2%*%tmp1%*%y
           -tmp2%*%tmp1%*%F2_tilde%*%C2
           -tmp2%*%phi12%*%C2)

    d <- solve(t(Ty)%*%Ty)%*%(t(Ty)%*%y-t(Ty)%*%F1_tilde%*%C1-t(Ty)%*%F2_tilde%*%C2)


    tmp1 <- (y-Ty%*%d)
    sigma2_y <- ((t(tmp1)%*%tmp1-t(tmp1)%*%F1_tilde*C1-t(tmp1)%*%F2_tilde*C2
                  -C1*t(F1_tilde)%*%tmp1+C1^2*phi11+C1*C2*phi12
                  -C2*t(F2_tilde)%*%tmp1+C1*C2*phi21+C2^2*phi22)/(n-qT))



    #
    #
    #     y_T_prim_bar <- t(y)%*%Ty/n
    #     inv_Ty_prim_Ty_bar <- solve(t(Ty)%*%Ty/n)
    #     F1_tilde_Ty_bar <- apply(F1_tilde*Ty,2,mean)
    #     F2_tilde_Ty_bar <- apply(F2_tilde*Ty,2,mean)
    #
    #     tmp <- t(y_T_prim_bar)%*%inv_Ty_prim_Ty_bar
    #     tmp1 <- (tmp%*%F1_tilde_Ty_bar)
    #     tmp2 <- (tmp%*%F2_tilde_Ty_bar)
    #
    #     C_a <- (t(F1_tilde_Ty_bar)%*%inv_Ty_prim_Ty_bar%*%F1_tilde_Ty_bar-phi1_tilde_bar)
    #     C_b <- (t(F2_tilde_Ty_bar)%*%inv_Ty_prim_Ty_bar%*%F2_tilde_Ty_bar-phi2_tilde_bar)
    #     C_c <- (t(F2_tilde_Ty_bar)%*%inv_Ty_prim_Ty_bar%*%F1_tilde_Ty_bar
    #             -(GAMMA[1,2]-mean(F1_tilde*F2_tilde)))
    #     C_d <- (t(F1_tilde_Ty_bar)%*%inv_Ty_prim_Ty_bar%*%F2_tilde_Ty_bar
    #             -(GAMMA[2,1]-mean(F2_tilde*F1_tilde)))
    #     C_e <- (tmp1-apply(F1_tilde*y,2,mean))
    #     C_f <- (tmp2-apply(F2_tilde*y,2,mean))
    #     denom_C <- solve(C_a*C_b-C_c*C_d)
    #
    #     C1 <- (denom_C%*%(C_e*C_b-C_f*C_d))
    #     C2 <- (denom_C%*%(C_a*C_f-C_c*C_e))
    #
    #     # doute ici
    #     D1 <- t( (X1_T1_prim_bar - t(A1)%*%t(F1_tilde_T1_bar))
    #              %*%inv_T1_prim_T1_bar)
    #     D2 <- t( (X2_T2_prim_bar - t(A2)%*%t(F2_tilde_T2_bar))
    #              %*%inv_T2_prim_T2_bar)
    #     d <- (y_T_prim_bar-C1*F1_tilde_Ty_bar-C2*F2_tilde_Ty_bar)%*%inv_Ty_prim_Ty_bar
    #
    #     mat <- (y-Ty%*%d)
    #     mat_t <- t(mat)
    #
    #     esp_f1f1 <- GAMMA[1,1]+t(F1_tilde)%*%F1_tilde
    #     esp_f2f2 <- GAMMA[2,2]+t(F2_tilde)%*%F2_tilde
    #     esp_f1f2 <- GAMMA[1,2]+t(F1_tilde)%*%F2_tilde
    #     esp_f2f1 <- GAMMA[2,1]+t(F2_tilde)%*%F1_tilde
    #     sigma2_y <- (mat_t%*%mat-C1*mat_t%*%F1_tilde-C2*mat_t%*%F2_tilde
    #                  -C1*t(F1_tilde)%*%mat+C1^2*esp_f1f1+C1*C2*esp_f1f2
    #                  -C2*t(F2_tilde)%*%mat+C1*C2*esp_f2f1+C2^2*esp_f2f2)


    #     sigma2_y <- (t(tmp1)%*%tmp1-t(tmp1)%*%F1_tilde*C1-t(tmp1)%*%F2_tilde*C2
    #                  -c1*t(F1_filde)%*%tmp1+c1^2*phi11+c1*c2*phi12
    #                  -c2*t(F2_tilde)%*%tmp1+c1*c2*phi21+c2^2*phi22)

    #sigma2y <- sum((y- Ty%*%d - C1*F1_tilde - C2*F2_tilde)^2)/n


    res <- (X1-T1%*%D1)
    sigma2_X1 <- sum(apply((res)^2,1,sum)
                     +sum(A1^2)*phi1_tilde
                     -2*diag(res%*%t(A1)%*%t(as.matrix(c(F1_tilde))))
    )/(n*qX1)

    res <- (X2-T2%*%D2)
    sigma2_X2 <- sum(apply((res)^2,1,sum)
                     +sum(A2^2)*phi2_tilde
                     -2*diag(res%*%t(A2)%*%t(as.matrix(c(F2_tilde))))
    )/(n*qX2)

    #Psi_Y <- sigma2_Y*diag(1,qY)
    Psi_X1 <- sigma2_X1*diag(1,qX1)
    Psi_X2 <- sigma2_X2*diag(1,qX2)

    D1_f <- D1
    D2_f <- D2
    d_f  <- d
    A1_f <- A1
    A2_f <- A2
    C1_f <- C1
    C2_f <- C2

    D1 <- as.matrix(D1)
    D2 <- as.matrix(D2)
    d  <- as.matrix(d)
    A1 <- as.numeric(A1)
    A2 <- as.numeric(A2)
    C1 <- as.numeric(C1)
    C2 <- as.numeric(C2)
    sigma2_y <- as.numeric(sigma2_y)

    # measure the change between iterations
    diff = sum(
      sum(abs( d -  d_it)/abs(d)) ,
      sum(abs(D1 -  D1_it)/abs(D1)),
      sum(abs(D2 -  D2_it)/abs(D2)),
      sum(abs( A1 -  A1_it)/abs(A1)),
      sum(abs( A2 -  A2_it)/abs(A2)) ,
      (abs( C1 -  C1_it)/abs(C1)) ,
      (abs( C2 -  C2_it)/abs(C2)) ,
      (abs(sigma2_y - sigma2_y_it)/abs(sigma2_y)) ,
      (abs(sigma2_X1 - sigma2_X1_it)/abs(sigma2_X1)) ,
      (abs(sigma2_X2 - sigma2_X2_it)/abs(sigma2_X2))
    )


    diffgraph <- c(diffgraph,diff)

    if(trace){
      cat("Nb of iteration:", it, "\n","Convergence:", diff, "\n")
    }

    it <- it+1

  } #end while


  namesTy <- c("intercept")
  namesT1 <- c("intercept")
  namesT2 <- c("intercept")
  if(qT1>1){
    for(i in 2:qT1){
      name <- paste("Var", i-1, sep = "")
      namesT1 <- c(namesT1,name)
    }
  }
  if(qT2>1){
    for(i in 2:qT2){
      name <- paste("Var", i-1, sep = "")
      namesT2 <- c(namesT2,name)
    }
  }
  if(qT>1){
    for(i in 2:qT){
      name <- paste("Var", i-1, sep = "")
      namesTy <- c(namesTy,name)
    }
  }
  rownames(d_f) <- namesTy
  rownames(D1_f) <- namesT1
  rownames(D2_f) <- namesT2

  Factors <- cbind(F1_tilde,F2_tilde)
  colnames(Factors) <- c("F1_predict","F2_predict")

  return(list(Factors=Factors,
              C=cbind(C1_f,C2_f),
              A1=A1_f,
              A2=A2_f,
              D=d_f,
              D1=D1_f,
              D2=D2_f,
              estimate_sigma2=cbind(sigma2_y,sigma2_X1,sigma2_X2),
              Convergence=diff,
              Iteration=it-1
  )
  )

}
