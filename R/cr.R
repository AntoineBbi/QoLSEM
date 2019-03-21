cr <- function(v_test) {
  if (class(v_test) != "numeric" ) cat( "Attention : l'argument n'est pas un vecteur il faut une class() numeric")
  esp <- mean(v_test)
  ectyp <- sqrt(var(v_test)*(length(v_test)-1)/length(v_test))
  v_test_centre_redui <- (v_test-rep(esp,length(v_test)))/rep(ectyp, length(v_test))
  return(v_test_centre_redui)
}
