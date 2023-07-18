library(matrixdist)
library(tidyverse)
library(highcharter)
library(Rcpp)


palette <- colorRampPalette(c("#261104","#A62454","#FA8F55","#FEE978","#80FF90","#80BEFF"))(7)
mon_theme <- hc_theme_merge(hc_theme_google(),hc_theme(colors=palette))

###############################################################################################################################################################
#Code C++
sourceCpp("C:/Users/aper/Documents/Visual Studio 2012/Projects/PH.cpp")

###############################################################################################################################################################
#Fonctions
random_structure <- function(dimension,structure){
  if (structure=="coxienne"){
    alpha <- c(1,rep(0,dimension-1))
    d1 <- -runif(dimension,0,2)
    d2 <- pmin(runif(dimension-1,0,1),-d1[-length(d1)])
    S <- matrix(0,nrow=dimension,ncol=dimension)
    diag(S) <- d1
    S[cbind(1:(nrow(S)-1), 2:ncol(S))] <- d2
  }
  if (structure=="coxienne généralisée"){
    alpha <- runif(dimension)
    alpha <- alpha/sum(alpha)
    d1 <- -runif(dimension,0,2)
    d2 <- pmin(runif(dimension-1,0,1),-d1[-length(d1)])
    S <- matrix(0,nrow=dimension,ncol=dimension)
    diag(S) <- d1
    S[cbind(1:(nrow(S)-1), 2:ncol(S))] <- d2
  }
  if (structure=="hyper-exponentielle"){
    alpha <- runif(dimension)
    alpha <- alpha/sum(alpha)
    d1 <- -runif(dimension,0,2)
    S <- matrix(0,nrow=dimension,ncol=dimension)
    diag(S) <- d1
  }
  if (structure=="erlang généralisée"){
    alpha <- c(1,rep(0,dimension-1))
    d1 <- -runif(dimension,0,2)
    d2 <- -d1[-length(d1)]
    S <- matrix(0,nrow=dimension,ncol=dimension)
    diag(S) <- d1
    S[cbind(1:(nrow(S)-1), 2:ncol(S))] <- d2
  }
  if (structure=="générale"){
    alpha <- runif(dimension)
    alpha <- alpha/sum(alpha)
    S <- matrix(runif(dimension^2),nrow=dimension,ncol=dimension)
    diag(S) <- -((apply(S,1,function(x) sum(x)-diag(S)) %>% diag) + rexp(dimension,2))
  }
  
  
  return(list(alpha=alpha,S=S))
  
}



#####################################
exit <- function(S){
  return(Exit(S))
}





#####################################
# EM_step <- function(x,y,weight){
#   
#   S <- x@pars$S
#   s <- exit(S)
#   alpha <- x@pars$alpha
#   
#   
#   ExpS <- function(y) expm(S*y)
# 
#   
# 
#   expS <- lapply(y,ExpS) #13/10000 sec
# 
#   
# 
#   
#   dimension <- length(alpha)
#   n <- length(y)
# 
#   
#   
#   #Pour le calcul de l'intégrale de l'exponentielle de matrice (Van Loan)
#   # G <- function(x){
#   #   expm((cbind(S,s %*% t(alpha)) %>% rbind(cbind(matrix(0,nrow=dimension,ncol=dimension),S)))*x)[(1:dimension),(dimension+1):(2*dimension)]
#   # }
#   # 
#   # Gy <- lapply(y, G) #64/10000 sec
# 
# 
#   Gy <- lapply(y, Van_loan,S=S,alpha=alpha,s=s) #34/10000 sec Version C++
# 
#   
#   # Etape E
# 
#   denom <- c(denominateur(expS,s,alpha))
# 
#   
#   
#   A <- lapply(1:dimension, function(k) {
# 
#     alpha[k] * sum(weight*(sapply(expS,function(x) x[k,] %*%s))/(denom))
#   }) #16/10000 sec
# 
# 
#   # Abis <- c(A_EM(alpha,weight,s,expS,denom))
#   
# 
#   B <- lapply(1:dimension, function(k) {
#     sum(weight*(sapply(Gy,function(x) x[k,k]))/denom)
#   }) #9/10000 sec
# 
#   # Bbis <- c(B_EM(dimension,Gy,weight,denom))
# 
#   
# 
#   C <- matrix(NA,nrow=dimension,ncol=dimension) 
# 
#   C <- sapply(1:dimension, function(k) {
#           sapply(1:dimension, function(l) {
#             C[k,l] <- sum(weight*S[k,l]*(sapply(Gy,function(x) x[l,k]))/denom)
#                                         })
#                                         
#           }) #33/10000 sec
#    
#   
# 
#   
# 
#   
# 
#   
#   diag(C) <- sapply(1:dimension, function(k) {
#     sum(weight* s[k] * (sapply(expS,function(x) alpha %*% x[,k])/denom))
#   }) #22/10000 sec
# 
#   
# 
#   # Cbis <- C_EM(dimension,S,Gy,weight,denom,s,expS,alpha)
# 
#   
# 
#   
#   A <- unlist(A)
#   B <- unlist(B)
# 
# 
#   
#   #Etape M
#   alpha <- A/n
#   alpha <- alpha/sum(alpha) 
#   s <- diag(C)/B 
#   S <- apply(C,1,function(x) x/B)
#   diag(S) <- 0
#   diag(S) <- -s+apply(S,1,function(x) -sum(x))
# 
#   return(PH(alpha,S))
#   
# }

EM_step <- function(x,y,weight,scale=1){
  y <- y/scale
  S <- x@pars$S
  s <- exit(S)
  alpha <- x@pars$alpha
  
  
  ExpS <- function(y) expm(S*y)
  
  
  
  expS <- lapply(y,ExpS) #13/10000 sec
  
  
  
  
  dimension <- length(alpha)
  n <- length(y)
  
  
  
  #Pour le calcul de l'intégrale de l'exponentielle de matrice (Van Loan)
  # G <- function(x){
  #   expm((cbind(S,s %*% t(alpha)) %>% rbind(cbind(matrix(0,nrow=dimension,ncol=dimension),S)))*x)[(1:dimension),(dimension+1):(2*dimension)]
  # }
  # 
  # Gy <- lapply(y, G) #64/10000 sec
  
  
  Gy <- lapply(y, Van_loan,S=S,alpha=alpha,s=s) #34/10000 sec Version C++
  
  
  # Etape E
  
  denom <- c(denominateur(expS,s,alpha))
  
  
  
  A <- c(A_EM(alpha,weight,s,expS,denom))
  
  
  
  B <- c(B_EM(dimension,Gy,weight,denom))
  
  
  
  
  
  C <- C_EM(dimension,S,Gy,weight,denom,s,expS,alpha)
  
  
  #Etape M
  alpha <- A/n
  alpha <- alpha/sum(alpha) 
  s <- diag(C)/B 
  S <- apply(C,1,function(x) x/B)
  diag(S) <- 0
  diag(S) <- -s+exit(S)
  if (class(x)=="PH"){
    return(PH(alpha,S/scale))}
  
  else if(class(x)=="IPH"){
    return(IPH(PH(alpha,S/scale),gfun=x@gfun$name,gfunpars = x@gfun$pars/scale))
  }
  
}
EM_step(obj3,y,weights)
#####################################
PMLE <- function(x){
  lambda=250
  q=4
  alpha <- abs(x[1:q])
  S <- matrix(x[(q+1):length(x)],ncol=q)
  
  
  
  G <- alpha %*% solve(-S) 
  G <- sort(G,decreasing = T)
  gamma <- 1/(1+sum((cumsum(G))/lead(G)-1,na.rm=T))
  PH <- ph(alpha,S)
  if ((sum(alpha)-1)<0.00001 & all(exit(S)>=0) & all((c(S)>=0)==c(diag(-1,q)>=0))){
    sum(log(dens(PH,data)))-lambda*gamma
  }
  else{
    -Inf
  }
}

###############################################################################################################################################################
#Définition du constructeur de classe PH
setClass("PH",slots = list(name="character",pars="list",structure="character",dimension="numeric",fit="list"))

PH <- function (alpha = NULL, S = NULL, structure = NULL, dimension = 3) 
{
  
  if (any(is.null(alpha)) & any(is.null(S)) & is.null(structure)) {
    stop("Entrez une matrice d'intensité et un vecteur pour la loi initiale ou entrez une structure et une dimension.")
  }
  if(is.null(structure) ) {
    if (dim(S)[1] != dim(S)[2]) {
      stop("S doit être une matrice carré.")
    }
    if (length(alpha) != dim(S)[1]) {
      stop("Dimensions incompatibles.")
    }
    structure <- "personnalisée"
  }
  else {
    if (!structure %in% c("coxienne","coxienne généralisée", "hyper-exponentielle", "erlang généralisée", "générale")){
      stop("Structure non-gérée.")
    }
    rs <- random_structure(dimension=dimension, structure = structure)
    alpha <- rs[[1]]
    S <- rs[[2]]
  }
  
  methods::new("PH", name = paste(structure, " PH(", length(alpha), 
                                  ")", sep = ""), pars = list(alpha = alpha, S = S),dimension=dimension,structure=structure)
}


###############################################################################################################################################################
#Définition du constructeur de classe IPH
setClass("IPH",slots = list(gfun="list"),contains = "PH")

IPH <- function (PH=NULL, gfun="gompertz",gfunpars=c(beta=1)) 
{
  
  if (is.null(PH)) {
    stop("Il faut définir une loi phase-type avant de créer un objet IPH.")
  }
  if(!gfun %in% c("gompertz") ) {
      stop("Fonction d'intensité non-gérée.")
  }
  if(gfun=="gompertz"){
    if (gfunpars<=0 | length(gfunpars)!=1){
      stop("Le paramètre de l'IPH devrait être strictement positif et de longueur 1.")
    }
    else{
      methods::new("IPH",gfun=list(name=gfun,pars=gfunpars,inverse=function(t,beta) (exp(t*beta)-1)/beta, inverse_prime= function(t,beta) exp(t*beta)*(t*beta-1)/beta^2,intensite=function(t,beta) exp(t*beta), intensite_prime=function(t,beta) t*exp(t*beta)),name = paste(PH@structure, " IPH(", length(PH@pars$alpha), 
                                                                              ")", sep = ""),PH)
    }
  }
}



###############################################################################################################################################################
#Méthodes PH
setMethod("show",
          "PH",
          function(object) {
            cat("Objet de la classe PH.", "\n", "\n")
            cat("Nom : Loi phase-type",object@name, "\n", "\n")
            cat("Paramètres :", "\n","@pars", "\n", "\n")
            cat("Loi initiale", "\n", "$alpha", "\n", "\n",object@pars$alpha, "\n", "\n")
            cat("Matrice d'intensité", "\n", "$S", "\n")
            print(object@pars$S)
          
          }
)

setMethod("fit",
          "PH",
          function(x,y,weight,stepsEM=500,every=10,scale=1) {
              obj <- x
              y <- y/scale
              WY <- list(weights=weight,obs=y)
              cat(format(Sys.time(),'%H:%M:%S')," : Début de l'algorithme EM","\n")
              for (i in 1:stepsEM){
                obj <- EM_step(obj,y,weight)
                obj@fit <- WY
                
                if (i %% every==0){
                  cat("\r",format(Sys.time(),'%H:%M:%S')," Etape",i," ")
                  cat("Log-Vraisemblance :",logLik(obj))}
                  flush.console()
              }
            obj@fit[["loglik"]] <- logLik(obj)
            cat("\n",format(Sys.time(),'%H:%M:%S')," : Fin de l'algorithme EM")
            obj <- scale(obj,scale=scale)
            return(obj)
          }
)

setMethod("scale",
          "PH",
          function(x,scale){
            obj <- x
            obj@pars$S <- x@pars$S/scale
            if (length(x@fit)>0){
              obj@fit$obs <- x@fit$obs*scale
            }
            return(obj)
          })

setMethod("logLik",
          "PH",
          function(object){
            weights <- object@fit$weights
            obs <- object@fit$obs
            S <- object@pars$S
            alpha <- object@pars$alpha
            s <- exit(object@pars$S)
            return(logvraisemblance(weights,obs,S,alpha,s))
            # sum(object@fit$weights*log((lapply(object@fit$obs,function(x) expm(object@pars$S*x)) %>% sapply(function(x) object@pars$alpha %*% x %*% exit(object@pars$S)))))
          })

setMethod("dens",
          "PH",
          function(x,y){
            S <- x@pars$S
            alpha <- x@pars$alpha
            s <- exit(S)
            Densite(alpha,S,s,y)
          })

setMethod("plot",
          "PH",
          function(x){
              df <- tibble(X=x@fit$obs,y=dens(x,X),weights=x@fit$weights)
              highchart() %>% 
                hc_add_series(df,"scatter",hcaes(x=X,y=weights),name="Observé") %>% 
                hc_add_series(df,"line",hcaes(x=X,y=y),name="Ajusté",marker=list(enabled=F)) %>% 
                hc_title(text=paste("Ajustement d'une distribution phase-type",x@name)) %>% 
                hc_add_theme(mon_theme)
            
          })

setMethod("coef",
          "PH",
          function(object){
            object@pars
          })

setMethod("pmle",
          "PH",
          function(x){
              opt <- optim(par=c(coef(x)$alpha, (coef(x)$S)), fn=PMLE, method="Nelder-Mead", control=list(maxit = 5000, reltol = 1e-6,fnscale=-1))
              q <- length(coef(x)$alpha)
              alpha <- abs(opt$par[1:q])
              S <- matrix(opt$par[(q+1):length(opt$par)],ncol=q)
              piU <- alpha %*% solve(-S)
              rew <- ifelse(piU/sum(piU)>=0.5/q,1,0)
              res <- matrixdist::TVR(ph(alpha,S),rew)
              res <- PH(res@pars$alpha,res@pars$S)
              res@pars$alpha <- res@pars$alpha/sum(res@pars$alpha)
              res@fit <- x@fit
              return(res)
              })

setMethod("exit",
          "PH",
          function(x){
            return(c(Exit(x@pars$S)))
          })

###############################################################################################################################################################
#Méthodes IPH

setMethod("show",
          "IPH",
          function(object) {
            cat("Objet de la classe IPH.", "\n", "\n")
            cat("Nom : Loi phase-type inhomogène",object@name, "\n", "\n")
            cat("Paramètres :", "\n","@pars", "\n", "\n")
            cat("Loi initiale", "\n", "$alpha", "\n", "\n",object@pars$alpha, "\n", "\n")
            cat("Matrice d'intensité", "\n", "$S", "\n")
            print(object@pars$S)
            cat("Paramètres de la fonction g :", "\n","@gfun", "\n", "\n")
            cat("Fonction d'intensité de type", object@gfun$name, "\n", "\n")
            cat("Paramètres de g", "\n", "$pars", "\n")
            print(object@gfun$pars)
          }
)

setMethod("scale",
          "IPH",
          function(x,scale){
            obj <- x
            obj@pars$S <- x@pars$S/scale
            if (length(x@fit)>0){
              obj@fit$obs <- x@fit$obs*scale
            }
            obj@gfun$pars <- x@gfun$pars/scale
            return(obj)
          })

setMethod("dens",
          "IPH",
          function(x,y){
            S <- x@pars$S
            alpha <- x@pars$alpha
            s <- exit(S)
            beta <- x@gfun$pars
            if (x@gfun$name=="gompertz"){
              DensiteGompertz(alpha,S,s,y,beta)
            }
            # x@gfun$intensite(y,pars) * x@pars$alpha %*% expm(S* x@gfun$inverse(y,pars)) %*% exit(S)
          })

setMethod("fit",
          "IPH",
          function(x,y,weight,stepsEM=500,every=10,scale=1) {
            obj <- x
            y <- y/scale
            WY <- list(weights=weight,obs=y)
            cat(format(Sys.time(),'%H:%M:%S')," : Début de l'algorithme EM","\n")
            for (i in 1:stepsEM){
              y_trans <- obj@gfun$inverse(y,obj@gfun$pars)
              obj2 <- IPH(EM_step(obj,y_trans,weight))
              if (obj@gfun$name=="gompertz"){
                famaxi <- function(beta){sum(weight*log(obj2@gfun$intensite(y,beta) *  sapply(y,function(x) obj2@pars$alpha %*% expm(obj2@pars$S* obj2@gfun$inverse(x,beta)) %*% exit(obj2@pars$S)))) }
                beta <- optim(obj@gfun$pars,famaxi,control = list(fnscale=-1))
              }
              else{
                stop("Fonction d'intensité non-gérée.")
              }
              
              obj <- obj2
              obj@gfun$pars <- beta$par
              obj@fit <- WY

              if (i %% every==0){
                cat("\r",format(Sys.time(),'%H:%M:%S')," Etape",i," ")
                cat("Log-Vraisemblance :",logLik(obj))}
              flush.console()
            }
            obj@fit[["loglik"]] <- logLik(obj)
            cat("\n",format(Sys.time(),'%H:%M:%S')," : Fin de l'algorithme EM")
            obj <- scale(obj,scale=scale)
            return(obj)
          }
)



###############################################################################################################################################################
#Test

obj <- PH(alpha=c(0.2,0.3,0.5),S=matrix(c(-1,0,0,0.2,-1,0,0,0.3,-2),nrow=3))
obj

y <- 1:100
weights <- dexp(1:100,1/10)


obj <- fit(obj,y,weights,100,every=10,scale=1)



obj2 <- ph(alpha=c(0.2,0.3,0.5),S=matrix(c(-1,0,0,0.2,-1,0,0,0.3,-2),nrow=3))
obj2 <- fit(obj2,y,weight=weights,stepsEM=100)


plot(y,weights/sum(weights))
curve(dens(obj,x),add=T,col='red')
curve(dens(obj2,x),add=T,col="green")

plot(obj)



obj3 <- PH(structure = "erlang généralisée",dimension=4)
y <- 1:500
weights <- dgamma(1:500,3,1/25)

obj3 <- fit(obj3,y,weights,1000,every=100)

plot(obj3)

obj4 <- pmle(obj3)
plot(obj4)

obj5 <- IPH(PH(structure = "coxienne",dimension=2))
y <- 1:30
# y <- sample(1:30,size=30,replace=F)
weights <- actuar::dgumbel(y,2,3)


obj6 <- fit(obj5,y,weights,stepsEM=1000,every=1,scale=100)
plot(obj6)

obj9 <- iph(ph(structure="general",dimension=10),gfun="gompertz")
obj10 <- fit(obj9,y/100,weights,stepsEM=1000)

obj10@pars$S <- obj10@pars$S/100
obj10@gfun$pars <- obj10@gfun$pars/100 

plot(y,weights)
curve(dens(obj10,x),col='red',add=T)

m <- 17
set.seed(1)
x1 <- IPH(PH(structure = "coxienne généralisée" , dimension = m),gfun="gompertz")
weights <- (donnees_ajust %>% pull(dx))/100000
y <- (donnees_ajust %>% pull(Age)) + 0.5
set.seed(1)
x_fit1 <- fit(x1, y, weight=weights, stepsEM = 200,every=10,scale=100)
plot(x_fit1)

#benchmark              
library(microbenchmark)
bench <- microbenchmark(fit(obj,y,weights,stepsEM=10,scale=10),fit(obj2,y,weight=weights,stepsEM=10,every=10))
bench
plot(bench)

bench2 <- microbenchmark(EM_step(obj,y,weights),fit(obj2,y,weight=weights,stepsEM=1),times=500)
bench2
plot(bench2)
