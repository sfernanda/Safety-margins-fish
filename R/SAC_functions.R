## ============ Functions for Species Accumulation Curve  ================= ##

curves <- function(k, z, x, type) {
  if(type=='gleason'){y <- k+log(x)*z}
  if(type=='gitay'){y <- k+z*(log(x))^2}
  if(type=='arrhenius'){y <- k*x^z}
  if(type=='michaelis-menten'){y <- k*x/(z+x)}
  return(y)
}

first.der <- function(k, z, x, type) {
  if(type=='gleason'){y <- z/x}
  if(type=='gitay'){y <- 2*z/x*log(x)}
  if(type=='arrhenius'){y <- k*z*x^(z-1)}
  if(type=='michaelis-menten'){y <- k*z/((z+x)^2)}
  return(y)
}
u=4