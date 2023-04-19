DirichletAlpha1 = function(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL, alpha4 = NULL, alpha5 = NULL, alpha6 = NULL, alpha7 = NULL) {
  
  y_check = function(y) {
    if ((is.matrix(y) && NCOL(y)!=7))
      stop("response should be a three-column matrix (y1, y2, y3, y4, y5, y6, y7) for this Dirichlet family")
  }
  
  loss = function(alpha2, alpha3, alpha4, alpha5, alpha6, alpha7, y, f = f) {
    y7 = y[,7]
    y6 = y[,6]
    y5 = y[,5]
    y4 = y[,4]
    y3 = y[,3]
    y2 = y[,2]
    y1 = y[,1]
    
    - (lgamma(exp(f) + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7) - (lgamma(exp(f)) + lgamma(alpha2) + lgamma(alpha3) + lgamma(alpha4) + lgamma(alpha5) + lgamma(alpha6) + lgamma(alpha7))
       + ((exp(f) - 1) * log(y1) + (alpha2 - 1) * log(y2) + (alpha3 - 1) * log(y3) + (alpha4 - 1) * log(y4) + (alpha5 - 1) * log(y5) + (alpha6 - 1) * log(y6) + (alpha7 - 1) * log(y7)))
    
  }
  
  
  
  
  ngradient = function(y, f, w = 1) {
    y7 = y[,7]
    y6 = y[,6]
    y5 = y[,5]
    y4 = y[,4]
    y3 = y[,3]
    y2 = y[,2]
    y1 = y[,1]
    
    ngr = exp(f) * (digamma(exp(f) + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + alpha7) -  digamma(exp(f)) +  log(y1))
    return(ngr)
  }
  
  
  
  
  risk = function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4, alpha5 = alpha5, alpha6 = alpha6, alpha7 = alpha7))
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
         offset = offset, name = "Septivariate Dirichlet distribution: alpha1 (log link)")
}


DirichletAlpha2 = function(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL, alpha4 = NULL, alpha5 = NULL, alpha6 = NULL, alpha7 = NULL) {

  
  loss = function(alpha1, alpha3, alpha4, alpha5, alpha6, alpha7, y, f = f) {
    y7 = y[,7]
    y6 = y[,6]
    y5 = y[,5]
    y4 = y[,4]
    y3 = y[,3]
    y1 = y[,1]
    y2 = y[,2]
  
    
    - (lgamma(alpha1 + exp(f) + alpha3 + alpha4 + alpha5 + alpha6 + alpha7) - (lgamma(alpha1) + lgamma(exp(f)) + lgamma(alpha3) + lgamma(alpha4) + lgamma(alpha5) + lgamma(alpha6) + lgamma(alpha7))
       + ((alpha1 - 1) * log(y1) + (exp(f) - 1) * log(y2) + (alpha3 - 1) * log(y3) + (alpha4 - 1) * log(y4) + (alpha5 - 1) * log(y5) + (alpha6 - 1) * log(y6) + (alpha7 - 1) * log(y7)))
    
  }
  
  
  
  
  ngradient = function(y, f, w = 1) {
    y7 = y[,7]
    y6 = y[,6]
    y5 = y[,5]
    y4 = y[,4]
    y3 = y[,3]
    y1 = y[,1]
    y2 = y[,2]
    
    ngr = exp(f) * (digamma(alpha1 + exp(f) + alpha3 + alpha4 + alpha5 + alpha6 + alpha7) -  digamma(exp(f)) +  log(y2))
    return(ngr)
  }
  
  
  
  
  risk = function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, alpha1 = alpha1, alpha3 = alpha3, alpha4 = alpha4, alpha5 = alpha5, alpha6 = alpha6, alpha7 = alpha7))
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
         offset = offset, name = "Septivariate Dirichlet distribution: alpha2 (log link)")
}



DirichletAlpha3 = function(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL, alpha4 = NULL, alpha5 = NULL, alpha6 = NULL, alpha7 = NULL) {
  
  
  loss = function(alpha1, alpha2, alpha4, alpha5, alpha6, alpha7, y, f = f) {
    y7 = y[,7]
    y6 = y[,6]
    y5 = y[,5]
    y4 = y[,4]
    y1 = y[,1]
    y2 = y[,2]
    y3 = y[,3]
    
    
    - (lgamma(alpha1 + alpha2 + exp(f) + alpha4 + alpha5 + alpha6 + alpha7) - (lgamma(alpha1) + lgamma(alpha2) + lgamma(exp(f)) + lgamma(alpha4) + lgamma(alpha5) + lgamma(alpha6) + lgamma(alpha7))
       + ((alpha1 - 1) * log(y1) + (alpha2 - 1) * log(y2) + (exp(f) - 1) * log(y3) + (alpha4 - 1) * log(y4) + (alpha5 - 1) * log(y5) + (alpha6 - 1) * log(y6) + (alpha7 - 1) * log(y7)))
    
  }
  
  
  
  
  ngradient = function(y, f, w = 1) {
    y7 = y[,7]
    y6 = y[,6]
    y5 = y[,5]
    y4 = y[,4]
    y1 = y[,1]
    y2 = y[,2]
    y3 = y[,3]
    
    ngr = exp(f) * (digamma(alpha1 + alpha2 + exp(f) + alpha4 + alpha5 + alpha6 + alpha7) -  digamma(exp(f)) +  log(y3))
    return(ngr)
  }
  
  
  
  
  risk = function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, alpha1 = alpha1, alpha2 = alpha2, alpha4 = alpha4, alpha5 = alpha5, alpha6 = alpha6, alpha7 = alpha7))
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
         offset = offset, name = "Septivariate Dirichlet distribution: alpha3 (log link)")
}



DirichletAlpha4 = function(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL, alpha4 = NULL, alpha5 = NULL, alpha6 = NULL, alpha7 = NULL) {
  
  
  loss = function(alpha1, alpha2, alpha3, alpha5, alpha6, alpha7, y, f = f) {
    y7 = y[,7]
    y6 = y[,6]
    y5 = y[,5]
    y1 = y[,1]
    y2 = y[,2]
    y3 = y[,3]
    y4 = y[,4]
    
    
    - (lgamma(alpha1 + alpha2 + alpha3 + exp(f) + alpha5 + alpha6 + alpha7) - (lgamma(alpha1) + lgamma(alpha2) + lgamma(alpha3) + lgamma(exp(f)) + lgamma(alpha5) + lgamma(alpha6) + lgamma(alpha7))
       + ((alpha1 - 1) * log(y1) + (alpha2 - 1) * log(y2) + (alpha3 - 1) * log(y3) + (exp(f) - 1) * log(y4) + (alpha5 - 1) * log(y5) + (alpha6 - 1) * log(y6) + (alpha7 - 1) * log(y7)))
    
  }
  
  
  
  
  ngradient = function(y, f, w = 1) {
    y7 = y[,7]
    y6 = y[,6]
    y5 = y[,5]
    y1 = y[,1]
    y2 = y[,2]
    y3 = y[,3]
    y4 = y[,4]
    
    ngr = exp(f) * (digamma(alpha1 + alpha2 + alpha3 + exp(f) + alpha5 + alpha6 + alpha7) -  digamma(exp(f)) +  log(y4))
    return(ngr)
  }
  
  
  
  
  risk = function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha5 = alpha5, alpha6 = alpha6, alpha7 = alpha7))
  }
  
  
  
  
  offset = function(y, w) {
    if (!is.null(alpha4)) {
      RET = log(alpha4)
    }
    else {
      RET = min(y[,4])
    }
    return(RET)
  }
  Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f),
         offset = offset, name = "Septivariate Dirichlet distribution: alpha4 (log link)")
}



DirichletAlpha5 = function(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL, alpha4 = NULL, alpha5 = NULL, alpha6 = NULL, alpha7 = NULL) {
  
  
  loss = function(alpha1, alpha2, alpha3, alpha4, alpha6, alpha7, y, f = f) {
    y7 = y[,7]
    y6 = y[,6]
    y1 = y[,1]
    y2 = y[,2]
    y3 = y[,3]
    y4 = y[,4]
    y5 = y[,5]
    
    
    - (lgamma(alpha1 + alpha2 + alpha3 + alpha4 + exp(f) + alpha6 + alpha7) - (lgamma(alpha1) + lgamma(alpha2) + lgamma(alpha3) + lgamma(alpha4) + lgamma(exp(f)) + lgamma(alpha6) + lgamma(alpha7))
       + ((alpha1 - 1) * log(y1) + (alpha2 - 1) * log(y2) + (alpha3 - 1) * log(y3) + (alpha4 - 1) * log(y4) + (exp(f) - 1) * log(y5) + (alpha6 - 1) * log(y6) + (alpha7 - 1) * log(y7)))
    
  }
  
  
  
  
  ngradient = function(y, f, w = 1) {
    y7 = y[,7]
    y6 = y[,6]
    y1 = y[,1]
    y2 = y[,2]
    y3 = y[,3]
    y4 = y[,4]
    y5 = y[,5]
    
    ngr = exp(f) * (digamma(alpha1 + alpha2 + alpha3 + alpha4 + exp(f) + alpha6 + alpha7) -  digamma(exp(f)) +  log(y5))
    return(ngr)
  }
  
  
  
  
  risk = function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4, alpha6 = alpha6, alpha7 = alpha7))
  }
  
  
  
  
  offset = function(y, w) {
    if (!is.null(alpha5)) {
      RET = log(alpha5)
    }
    else {
      RET = min(y[,5])
    }
    return(RET)
  }
  Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f),
         offset = offset, name = "Septivariate Dirichlet distribution: alpha5 (log link)")
}





DirichletAlpha6 = function(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL, alpha4 = NULL, alpha5 = NULL, alpha6 = NULL, alpha7 = NULL) {
  
  
  loss = function(alpha1, alpha2, alpha3, alpha4, alpha5, alpha7, y, f = f) {
    y7 = y[,7]
    y1 = y[,1]
    y2 = y[,2]
    y3 = y[,3]
    y4 = y[,4]
    y5 = y[,5]
    y6 = y[,6]
    
    
    - (lgamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + exp(f) + alpha7) - (lgamma(alpha1) + lgamma(alpha2) + lgamma(alpha3) + lgamma(alpha4) + lgamma(alpha5) + lgamma(exp(f)) + lgamma(alpha7))
       + ((alpha1 - 1) * log(y1) + (alpha2 - 1) * log(y2) + (alpha3 - 1) * log(y3) + (alpha4 - 1) * log(y4) + (alpha5 - 1) * log(y5) + (exp(f) - 1) * log(y6) + (alpha7 - 1) * log(y7)))
    
  }
  
  
  
  ngradient = function(y, f, w = 1) {
    y7 = y[,7]
    y1 = y[,1]
    y2 = y[,2]
    y3 = y[,3]
    y4 = y[,4]
    y5 = y[,5]
    y6 = y[,6]
    
    ngr = exp(f) * (digamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + exp(f) + alpha7) -  digamma(exp(f)) +  log(y6))
    return(ngr)
  }
  
  
  
  
  risk = function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4, alpha5 = alpha5, alpha7 = alpha7))
  }
  
  
  
  
  offset = function(y, w) {
    if (!is.null(alpha6)) {
      RET = log(alpha6)
    }
    else {
      RET = min(y[,6])
    }
    return(RET)
  }
  Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f),
         offset = offset, name = "Septivariate Dirichlet distribution: alpha6 (log link)")
}




DirichletAlpha7 = function(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL, alpha4 = NULL, alpha5 = NULL, alpha6 = NULL, alpha7 = NULL) {
  
  
  loss = function(alpha1, alpha2, alpha3, alpha4, alpha5, alpha6, y, f = f) {
    y1 = y[,1]
    y2 = y[,2]
    y3 = y[,3]
    y4 = y[,4]
    y5 = y[,5]
    y6 = y[,6]
    y7 = y[,7]
    
    
    - (lgamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + exp(f)) - (lgamma(alpha1) + lgamma(alpha2) + lgamma(alpha3) + lgamma(alpha4) + lgamma(alpha5) + lgamma(alpha6) + lgamma(exp(f)))
       + ((alpha1 - 1) * log(y1) + (alpha2 - 1) * log(y2) + (alpha3 - 1) * log(y3) + (alpha4 - 1) * log(y4) + (alpha5 - 1) * log(y5) + (alpha6 - 1) * log(y6) + (exp(f) - 1) * log(y7)))
    
  }
  
  
  
  
  ngradient = function(y, f, w = 1) {
    y1 = y[,1]
    y2 = y[,2]
    y3 = y[,3]
    y4 = y[,4]
    y5 = y[,5]
    y6 = y[,6]
    y7 = y[,7]
    
    ngr = exp(f) * (digamma(alpha1 + alpha2 + alpha3 + alpha4 + alpha5 + alpha6 + exp(f)) -  digamma(exp(f)) +  log(y7))
    return(ngr)
  }
  
  
  
  
  risk = function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4, alpha5 = alpha5, alpha6 = alpha6))
  }
  
  
  
  
  offset = function(y, w) {
    if (!is.null(alpha7)) {
      RET = log(alpha7)
    }
    else {
      RET = min(y[,7])
    }
    return(RET)
  }
  Family(ngradient = ngradient, risk = risk, loss = loss, response = function(f) exp(f),
         offset = offset, name = "Septivariate Dirichlet distribution: alpha7 (log link)")
}




DirichletSV = function(alpha1 = NULL, alpha2 = NULL, alpha3 = NULL, alpha4 = NULL, alpha5 = NULL, alpha6 = NULL, alpha7 = NULL) {
  
  if ((!is.null(alpha1) && alpha1 <= 0)) 
    stop(sQuote("alpha1"), " must be greater than zero")
  if ((!is.null(alpha2) && alpha2 <= 0)) 
    stop(sQuote("alpha2"), " must be greater than zero")
  if ((!is.null(alpha3) && alpha3 <= 0)) 
    stop(sQuote("alpha3"), " must be greater than zero")
  if ((!is.null(alpha4) && alpha4 <= 0)) 
    stop(sQuote("alpha4"), " must be greater than zero")
  if ((!is.null(alpha5) && alpha5 <= 0)) 
    stop(sQuote("alpha5"), " must be greater than zero")
  if ((!is.null(alpha6) && alpha6 <= 0)) 
    stop(sQuote("alpha6"), " must be greater than zero")
  if ((!is.null(alpha7) && alpha7 <= 0)) 
    stop(sQuote("alpha7"), " must be greater than zero")
  
  Families(alpha1 = DirichletAlpha1(alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4, alpha5 = alpha5, alpha6 = alpha6, alpha7 = alpha7),
           alpha2 = DirichletAlpha2(alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4, alpha5 = alpha5, alpha6 = alpha6, alpha7 = alpha7),
           alpha3 = DirichletAlpha3(alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4, alpha5 = alpha5, alpha6 = alpha6, alpha7 = alpha7),
           alpha4 = DirichletAlpha4(alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4, alpha5 = alpha5, alpha6 = alpha6, alpha7 = alpha7),
           alpha5 = DirichletAlpha5(alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4, alpha5 = alpha5, alpha6 = alpha6, alpha7 = alpha7),
           alpha6 = DirichletAlpha6(alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4, alpha5 = alpha5, alpha6 = alpha6, alpha7 = alpha7),
           alpha7 = DirichletAlpha7(alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4, alpha5 = alpha5, alpha6 = alpha6, alpha7 = alpha7),
           name = "DirichletSV")
}
