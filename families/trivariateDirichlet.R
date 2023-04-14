DirichletAlpha1 = function(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL) {
  
  y_check = function(y) {
    if ((is.matrix(y) && NCOL(y)!=3))
      stop("response should be a three-column matrix (y1, y2 and y3) for this Dirichlet family")
  }
  
  loss = function(alpha2 , alpha3 , y , f = f) {
    y3 = y[,3]
    y2 = y[,2]
    y1 = y[,1]
    
    - (lgamma(exp(f) + alpha2 + alpha3) - (lgamma(exp(f)) + lgamma(alpha2) + lgamma(alpha3))
    + ((exp(f) - 1) * log(y1) + (alpha2 - 1) * log(y2) + (alpha3 - 1) * log(y3)))
    
  }
  
  
  
  
  ngradient = function(y, f, w = 1) {
    y3 = y[,3]
    y2 = y[,2]
    y1 = y[,1]
    
    ngr = exp(f) * (digamma(exp(f) + alpha2 + alpha3) -  digamma(exp(f)) +  log(y1))
    return(ngr)
  }
  
  
  
  
  risk = function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, alpha2 = alpha2, alpha3 = alpha3))
  }
  
  
  
  
  offset = function(y, w) {
    if (!is.null(alpha1)) {
      RET = log(alpha1)
    }
    else {
      RET = min(y[,1])
    }
    return(RET)
  }
  Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f),
         offset = offset, name = "Trivariate Dirichlet distribution: alpha1 (log link)")
}



DirichletAlpha2 = function(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL) {
  
  loss = function(alpha1 , alpha3 , y , f = f) {
    y3 = y[,3]
    y1 = y[,1]
    y2 = y[,2]
    
    - (lgamma(alpha1 + exp(f) + alpha3) - (lgamma(alpha1) + lgamma(exp(f)) + lgamma(alpha3))
      + ((alpha1 - 1) * log(y1) + (exp(f) - 1) * log(y2) + (alpha3 - 1) * log(y3)))
    
  }
  
  
  
  
  ngradient = function(y, f, w = 1) {
    y3 = y[,3]
    y1 = y[,1]
    y2 = y[,2]
    
    ngr = exp(f) * (digamma(alpha1 + exp(f) + alpha3) - digamma(exp(f)) + log(y2))
    return(ngr)
  }
  
  
  
  
  risk = function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, alpha1 = alpha1, alpha3 = alpha3))
  }
  
  
  
  
  offset = function(y, w) {
    if (!is.null(alpha2)) {
      RET = log(alpha2)
    }
    else {
      RET = min(y[,2])
    }
    return(RET)
  }
  Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f),
         offset = offset, name = "Trivariate Dirichlet distribution: alpha2 (log link)")
}



DirichletAlpha3 = function(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL) {
  
  loss = function(alpha1 , alpha2 , y , f = f) {
    y1 = y[,1]
    y2 = y[,2]
    y3 = y[,3]
    
    - (lgamma(alpha1 + alpha2 + exp(f)) - (lgamma(alpha1) + lgamma(alpha2) + lgamma(exp(f)))
      + ((alpha1 - 1) * log(y1) + (alpha2 - 1) * log(y2) + (exp(f) - 1) * log(y3)))
    
  }
  
  
  
  
  ngradient = function(y, f, w = 1) {
    y1 = y[,1]
    y2 = y[,2]
    y3 = y[,3]
    
    ngr = exp(f) * (digamma(alpha1 + alpha2 + exp(f)) - digamma(exp(f)) + log(y3))
    return(ngr)
  }
  
  
  
  
  risk = function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, alpha1 = alpha1, alpha2 = alpha2))
  }
  
  
  
  
  offset = function(y, w) {
    if (!is.null(alpha3)) {
      RET = log(alpha3)
    }
    else {
      RET = min(y[,3])
    }
    return(RET)
  }
  Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f),
         offset = offset, name = "Trivariate Dirichlet distribution: alpha3 (log link)")
}


DirichletTV = function(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL) {
  
  if ((!is.null(alpha1) && alpha1 <= 0)) 
    stop(sQuote("alpha1"), " must be greater than zero")
  if ((!is.null(alpha2) && alpha2 <= 0)) 
    stop(sQuote("alpha2"), " must be greater than zero")
  if ((!is.null(alpha3) && alpha3 <= 0)) 
    stop(sQuote("alpha3"), " must be greater than zero")
  
  Families(alpha1 = DirichletAlpha1(alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3),
           alpha2 = DirichletAlpha2(alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3),
           alpha3 = DirichletAlpha3(alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3),
           name = "DirichletTV")
  
  
}



