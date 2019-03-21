#' Generation dataset from a QoLSEM model
#'
#' This function is used to generate dataset(s) from the model proposed by Barbieri, Tami et al. (2018).
#' It is used to propose an example to the users.
#'
#' @param N number of dataset to generate
#' @param I number of observations (statistic units) by datasets; a pair number are expected
#' @param c a numerical vector denoting the two scalars c=[c1,c2]: parameters associated with the two factors to explain the response varaible Y
#' @param a1 a numerical vector (of size q_1) of parameters associated with the factor f_1
#' @param a2 a numerical vector (of size q_2) of parameters associated with the factor f_2
#' @param d vector of explanatory variables associated with y
#' @param D1 matrix of parameters associated with the block 1
#' @param D2 matrix of parameters associated with the block 2
#' @param sigma.y standard deviation of the response variable y
#' @param sigma.X1 average standard deviation of the block 1
#' @param sigma.X2 average standard deviation of the block 2
#'
#' @return A dataframe for N>1 and a list for N=1
#'
#' @author Antoine Barbieri, Myriam Tami
#'
#' @references Barbieri A, Tami M, Bry X, Azria D, Gourgou S, Mollevi C, Lavergne C. (2018) \emph{EM algorithm estimation of a structural equation model for the longitudinal study of the quality of life}. Statistics in Medicine. 37(6) :1031-1046.
#'
#' @examples
#' ## test avec N=1
#' test1 <- generation.QoLSEM(N=1,I=150,
#'                            c=as.matrix(c(2,-2)),
#'                            a1=matrix(c(1:7),nrow=7),
#'                            a2=matrix(c(1:12),nrow=12),
#'                            d=c(80),
#'                            D1=matrix(seq(2,14,2)+70,ncol=7),
#'                            D2=matrix(seq(2,24,2)+20,ncol=12),
#'                            sigma.y=10,sigma.X1=10,sigma.X2=10)
#'
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
#' @export
generation.QoLSEM <- function(N=500,
                              I=150,
                              c=as.matrix(c(2,-2)),
                              a1=c(1:7),
                              a2=c(1:12),
                              d=as.matrix(c(80)),
                              D1=matrix(seq(0:14,2)+70,ncol=7),
                              D2=matrix(seq(0:24,2)+20,ncol=12),
                              sigma.y=10,
                              sigma.X1=10,
                              sigma.X2=10){
	#### WARNING
	# the number of column for D_b (b=1,2) has to be equal to the length of a_b : the number of variable X_b

#
# 	#### fonction(s) necessaires(s)
# 	cr <- function(v_test) {
# 		if (class(v_test) != "numeric" ) cat( "Attention : l'argument n'est pas un vecteur il faut une class() numeric")
# 		esp <- mean(v_test)
# 		ectyp <- sqrt(var(v_test)*(length(v_test)-1)/length(v_test))
# 		v_test_centre_redui <- (v_test-rep(esp,length(v_test)))/rep(ectyp, length(v_test))
# 		return(v_test_centre_redui)
# 	}

  simulated_data <- NULL

  for(n in 1:N){


  	#### Creation des matrices des erreurs
  	#
  	varepsilon.y <- as.matrix(rnorm(I,0,sigma.y))
  	#
  	if(length(sigma.X1)!=1){
  		l <- length(sigma.X1)
  		varepsilon.X1 <- rnorm(I,0,sigma.X1[1])
  		for(i in 2:l){varepsilon.X1 <- cbind(varepsilon.X1,rnorm(I,0,sigma.X1[l]))}
  	}else{varepsilon.X1 <- matrix(rnorm(I*length(a1),0,sigma.X1),ncol=length(a1))}
  	#
  	if(length(sigma.X2)!=1){
  		l <- length(sigma.X2)
  		varepsilon.X2 <- rnorm(I,0,sigma.X2[1])
  		for(i in 2:l){varepsilon.X2 <- cbind(varepsilon.X2,rnorm(I,0,sigma.X2[l]))}
  	}else{varepsilon.X2 <- matrix(rnorm(I*length(a2),0,sigma.X2),ncol=length(a2))}
  	#------------------

  	#### Creation des facteurs
  	f1 <- rnorm(I)
  	f2 <- rnorm(I)
  	#Normalisation des facteurs pour avoir une variance de 1 exactement
  	f1 <- as.numeric(f1)
  	f1 <-cr(f1)
  	f1 <- as.matrix(f1)
  	f2 <- as.numeric(f2)
  	f2 <-cr(f2)
  	f2 <- as.matrix(f2)

  	#CREATION DE T, T1 et T2
    t  <- as.matrix(rep(1,I))
    T1 <- matrix(1,ncol=1,nrow=I)
  	T2 <- matrix(1,ncol=1,nrow=I)


    #Modele QoLeG
    X1 <- T1 %*% D1 + f1%*%t(a1) + varepsilon.X1
    X2 <- T2 %*% D2 + f2%*%t(a2) + varepsilon.X2
    y  <- t %*%  d  + f1*c[1,] + f2*c[2,] + varepsilon.y


    simulated_data <- rbind(simulated_data,
                            cbind(rep(n,I),f1,f2,y,X1,X2))


  }

	if(N==1){
	  return(list("y"=y,
	              "X1"=X1,
	              "X2"=X2,
	              "f1"=f1,
	              "f2"=f2)
           )
	} else {
    return(simulated_data)
	}
}
