DirichletAlpha = function(k=k,
                          a1=NULL,a2=NULL,a3=NULL,a4=NULL,a5=NULL,a6=NULL,a7=NULL,a8=NULL,a9=NULL,
                          K=K) {
  
  arg = c("y",paste0("a",1:K))
  arg[arg %in% paste0("a",k)] = "f"
  
  # loss
  
  eval(parse(text=paste0("loss = function(",paste(arg,collapse  = ","),"){
       A =cbind(",gsub("f","exp(f)",paste0(grep("a|f",arg,value=TRUE),collapse=",")),")
       return ( - ( lgamma(rowSums(A)) - rowSums(lgamma(A)) + rowSums((A-1)*log(y)) ))}")))
  
  # ngradient
  
  eval(parse(text=paste0("ngradient = function(y,f,w=1){
       A =cbind(",gsub("f","exp(f)",paste0(grep("a|f",arg,value=TRUE),collapse=",")),")

       ngr = A[,k] * (digamma(rowSums(A)) - digamma(A[,k]) + log(y[,k]))

       return(ngr)}")))
  
  
  # risk
  
  eval(parse(text=paste0("risk = function(y,f,w=1){
       sum(w * loss(",paste(arg,arg,sep = "=",collapse=","),"))}")))
  
  
  offset = function(y, w) {
    
    if (!is.null(get(paste0("a",k)))) {
      RET = log(get(paste0("a",k)))
    }
    else {
      RET = min(y[,k])
    }
    return(RET)
  }
  Family(ngradient = ngradient, risk = risk, loss = loss, response =
           function(f) exp(f),
         offset = offset, name = paste0("Dirichlet Distribution: a",k," 
(log link)"))
}

Dirichlet =
  function(K,a1=NULL,a2=NULL,a3=NULL,a4=NULL,a5=NULL,a6=NULL,a7=NULL,a8=NULL,a9=NULL,
           
           m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL,m8=NULL,m9=NULL,
           
           s1=NULL,s2=NULL,s3=NULL,s4=NULL,s5=NULL,s6=NULL,s7=NULL,s8=NULL,s9=NULL,
           ...){
    
    pars = paste(paste0("a",1:K),collapse = ",")
    
    a = paste(sapply(1:K,function(k){
      paste0("a",k,"=DirichletAlpha(k=",k,",",pars,",K=K)")
    }),collapse=",")
    
    
    eval(parse(text=paste0("Families(",paste(a,sep=","),")")))
    
    
  }
