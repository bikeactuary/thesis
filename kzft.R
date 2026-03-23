coeff <- function(m, k){as.vector(polynom::polynomial(rep(1,m))^k/m^k)}

kzft <- function(x, f=0, m=1, k=5)
{
  h <- floor((m-1)/2)
  t  <- ceiling((m-1)/2)
  x <- append(x, rep(NA,t*k), after = length(x) )
  x <- append(x, rep(NA,h*k), after = 0)
  
  z <- as.complex(x)
  
  for(i in 1:k) {
    s     <- rep(0, length(z))
    count <- rep(0, length(z))
    
    count <- count + !is.na(z)
    z[is.na(z)] <- 0
    s <- as.complex(s) + z
    
    for(j in 1:h) {
      # add offset
      z   <- c(rep(NA, j), x[1:(length(x)-j)])
      z <- z * exp(complex(real = 0, imaginary = 2*pi*f*j))
      count <- count + !is.na(z)
      z[is.na(z)] <- 0
      s <- s + z
    }
    for(j in 1:t) {
      z   <- c(x[(j+1):length(x)], rep(NA, j))
      z <- z * exp(complex(real = 0, imaginary = -2*pi*f*j))
      count <- count + !is.na(z)
      z[is.na(z)] <- 0
      s <- s + z
    }
    z <- s/count
    x<-z
  }
  z <- z[(t*k+1):(length(z)-(h*k))]
  return(z);
}

kzp <- function(y, m=length(y), k=1, f_range = NULL, ncores = 1)
{
  M<-(m-1)*k+1
  n=length(y)
  
  ## default to full spectrum as determined by m
  ## but allow user-specified range (provided in frequency-bin)
  if (is.null(f_range)) {
    f_est <- c(0,m-1)
  } else {
    f_est <- f_range
  }
  
  ## parallelized over the frequencies to be evaluate
  z_tmp <- parallel::mclapply(
    X = f_est[1]:f_est[2],
    FUN = function(i, x, m, k) kzft(x,f=i/m,m,k=k),
    x = y, m = m, k=k,
    mc.cores = ncores
  )
  z <- matrix(unlist(z_tmp), nrow = n, ncol = f_est[2]-f_est[1]+1)
  
  d <- apply(z, 2, function(z) {(abs(z)^2)*M})
  a <- colMeans(d, na.rm=TRUE)
  
  if(is.null(f_range)) {
    a <- a[1:(m/2)]
    f_range <- c(0, length(a)-1)
  } else {
    a <- a[1:(f_range[2] - f_range[1] + 1)]
  }
  b <- 2*sqrt(a/M)
  c <- b^2
  
  structure(list(
    periodogram = a, ## original
    energy = b, ## estimated energy
    power = c, ## estimated power
    window=m,
    f_range=f_range, ## for user-specified frequency ranges
    k=k,
    var=var(y),
    smooth_periodogram=NULL,
    smooth_method=NULL,
    call=match.call()
  ),
  class = "kzp")
}


plot.kzp <- function(x, remove0 = TRUE, scale = NA, ylab = '', ...) {
  if (!is.na(scale)) {
    stopifnot("scale must be one of 'original', 'power', or 'energy'" = scale %in% c("original", "power", "energy")) 
  }
  
  if(is.na(scale) | scale == "original") {
    p <- x$periodogram
  } else if (scale == "energy") {
    p <- x$energy
  } else if (scale == "power") {
    p <- x$power
  }
  
  if(x$f_range[1] == 0 & remove0) {
    dz <- p[-1]
    omega <- ((x$f_range[1]+1):x$f_range[2])/x$window
    plot(omega, dz, type="l", xlab="Frequency", ylab=ylab)
  } else  {
    dz <- p
    omega <- (x$f_range[1]:x$f_range[2])/x$window
    plot(omega, dz, type="l", xlab="Frequency", ylab=ylab)
  }
}

