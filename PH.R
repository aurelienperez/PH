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
  if (dimension==1){
    alpha <- 1
    S <- matrix(-runif(1))
  }
  else{
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
    diag(S) <- -((apply(S,1,function(x) sum(x)-diag(S)) %>% diag) + runif(dimension))
  }}
  
  
  return(list(alpha=alpha,S=S))
  
}



#####################################
exit <- function(x){
  return(Exit(x))
}


#####################################
surv <- function(x,y,...){
  Survie(x$alpha,x$S,y)
}




#####################################
EM_step <- function(x,y,weight,scale=1){
  y <- y/scale
  S <- x@pars$S
  alpha <- x@pars$alpha
  
  
  dimension <- length(alpha)
  n <- length(y)
  

  
  # Etape E
  
  ABC <- ABC(alpha,S,y,weight)
  
  
  A <- ABC[[1]]
  
  B <- ABC[[2]]
  
  C <- t(ABC[[3]][,1:dimension])
  
  diag(C) <- ABC[[3]][,(dimension+1)]

  
  #Etape M
  alpha <- A/n
  alpha <- alpha/sum(alpha) 
  s <- diag(C)/B 
  S <- apply(C,1,function(x) x/B)
  if (dimension==1){
    S <- matrix(-s)
  }
  else{
    diag(S) <- 0
    diag(S) <- -s+exit(S)
  }
  
  if (class(x)=="PH"){
    return(PH(alpha,S/scale))}
  
  else if(class(x)=="IPH"){
    return(IPH(PH(alpha,S/scale),gfun=x@gfun$name,gfunpars = x@gfun$pars/scale))
  }
  else if(class(x)=="SPH"){
    return(SPH(IPH=IPH(PH=PH(alpha,S/scale),gfun=x@gfun$name,gfunpars = x@gfun$pars/scale),regpars=x@reg$pars))
  }
  
}


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
#Définition du constructeur de classe SPH
setClass("SPH",slots = list(reg="list"),contains = "IPH")

SPH <- function(IPH=NULL,regpars=numeric(0)){
  methods::new("SPH",reg=list(pars=regpars),name = paste(IPH@structure, " SPH(", length(IPH@pars$alpha), 
                                                         ")", sep = ""), IPH)
                                                                                                                                                                                                                                                                         
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
                  cat("Log-Vraisemblance :",logLik(scale(obj,scale=scale)))}
                  flush.console()
              }
            obj@fit[["loglik"]] <- logLik(scale(obj,scale=scale))
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
            return(sum(weights*log(dens(object,obs))))
            
          })

setMethod("dens",
          "PH",
          function(x,y){
            S <- x@pars$S
            alpha <- x@pars$alpha
            s <- exit(S)
            c(Densite(alpha,S,s,y))
          })

setMethod("surv",
          "PH",
          function(x,y){
            S <- x@pars$S
            alpha <- x@pars$alpha
            Survie(alpha,S,y)
          })

setMethod("plot",
          "PH",
          function(x,type="dens",f="identity",hch=highchart()){
            if (type %in% c("dens","surv","haz")){
              df <- tibble(X=x@fit$obs,y=get(f)(get(type)(x,X)),weights=x@fit$weights)
              hc <- hch %>% 
                hc_title(text=paste("Ajustement d'une distribution phase-type",x@name)) %>% 
                hc_add_theme(mon_theme)
              if (type=="dens"){
                hc <- hc %>%  hc_add_series(df,"scatter",hcaes(x=X,y=weights),name="Observé") 
              }
              hc <- hc %>% 
                hc_add_series(df,"line",hcaes(x=X,y=y),name="Ajusté",marker=list(enabled=F)) 
                }
            
            hc
          })

setMethod("coef",
          "PH",
          function(object){
            object@pars
          })

# setMethod("pmle",
#           "PH",
#           function(x){
#               opt <- optim(par=c(coef(x)$alpha, (coef(x)$S)), fn=PMLE, method="Nelder-Mead", control=list(maxit = 5000, reltol = 1e-6,fnscale=-1))
#               q <- length(coef(x)$alpha)
#               alpha <- abs(opt$par[1:q])
#               S <- matrix(opt$par[(q+1):length(opt$par)],ncol=q)
#               piU <- alpha %*% solve(-S)
#               rew <- ifelse(piU/sum(piU)>=0.5/q,1,0)
#               res <- matrixdist::TVR(ph(alpha,S),rew)
#               res <- PH(res@pars$alpha,res@pars$S)
#               res@pars$alpha <- res@pars$alpha/sum(res@pars$alpha)
#               res@fit <- x@fit
#               return(res)
#               })

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
              c(DensiteGompertz(alpha,S,s,y,beta))
            }
          })

setMethod("surv",
          "IPH",
          function(x,y){
            S <- x@pars$S
            alpha <- x@pars$alpha
            beta <- x@gfun$pars
            if (x@gfun$name=="gompertz"){
              c(SurvieGompertz(alpha,S,y,beta))
            }
          })

setMethod("haz",
          "IPH",
          function(x,y){
            dens(x,y)/surv(x,y)
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
              obj2 <- IPH(EM_step(obj,y_trans,weight),gfunpars = obj@gfun$pars)
              if (obj@gfun$name=="gompertz"){
                LL <- function(beta,alpha,S,s,y,weight){ 
                  sum(weight*c(logDensiteGompertz(alpha,S,s,y,beta)))}
                beta <- optim(par=obj@gfun$pars,fn=LL,alpha=obj2@pars$alpha,S=obj2@pars$S,s=exit(obj2),weight=weights,y=y,control = list(fnscale=-1,warn.1d.NelderMead=F))
                }
              else{
                stop("Fonction d'intensité non-gérée.")
              }
              
              obj <- obj2
              obj@gfun$pars <- beta$par
              obj@fit <- WY

              if (i %% every==0){
                cat("\r",format(Sys.time(),'%H:%M:%S')," Etape",i," ")
                cat("Log-Vraisemblance :",logLik(scale(obj,scale=scale)))}
              flush.console()
            }
            obj@fit[["loglik"]] <- logLik(scale(obj,scale=scale))
            cat("\n",format(Sys.time(),'%H:%M:%S')," : Fin de l'algorithme EM","\n")
            obj <- scale(obj,scale=scale)
            return(obj)
          }
)



setGeneric("reg", function(object,formule,data,weight,stepsEM=500,every=10,scale=1) {
  standardGeneric("reg")
})

setMethod("reg",
          "IPH",
          function(object,formule,data,weight,stepsEM=500,every=10,scale=1){
            y <- data[[as.character(attr(terms(formule), which = "variables")[[2]])]]
            y <- y/scale
            X <- data %>% select(attr(terms(formule), which = "term.labels")) %>% as.matrix()
            obj <- SPH(object,rep(1,dim(X)[2]))
            WY <- list(weights=weight,obs=y)
            cat(format(Sys.time(),'%H:%M:%S')," : Début de l'algorithme EM","\n")
            for (i in 1:stepsEM){
              y_trans <- exp(X %*% obj@reg$pars) * obj@gfun$inverse(y,obj@gfun$pars)
              obj2 <- EM_step(obj,y_trans,weight)
              if (obj@gfun$name=="gompertz"){
                LL <- function(params,alpha,S,y,weight) {
                  beta <- params[1]
                  theta <- params[2:length(params)]
                  mx <- exp(X %*% theta)
                  sum(weight*log(c(DensiteGompertzs(alpha,S,y,beta,mx))))
                  
                }
                LL2 <- function(params,alpha,S,y,weight) {
                  beta <- params[1]
                  theta <- params[2:length(params)]
                  l_mx <- X %*% theta
                  sum(weight*c(logDensiteGompertzs(alpha,S,y,beta,l_mx)))
                  
                }
                beta <- optim(par=c(obj@gfun$pars,obj@reg$pars),fn=LL2,alpha=obj2@pars$alpha,S=obj2@pars$S,weight=weights,y=y,control = list(fnscale=-1,warn.1d.NelderMead=F))
                }
              
              else{
                stop("Fonction d'intensité non-gérée.")
              }
              
              obj <- obj2
              obj@gfun$pars <- beta$par[1]
              obj@reg$pars <- beta$par[2:length(beta$par)]
              # obj@gfun$pars <- c(beta=beta$maximum)
              obj@fit <- WY
              
              if (i %% every==0){
                cat("\r",format(Sys.time(),'%H:%M:%S')," Etape",i," ")
                cat("Log-Vraisemblance :",beta$value)}
              flush.console()
            }
            # obj@fit[["loglik"]] <- logLik(scale(obj,scale=scale))
            cat("\n",format(Sys.time(),'%H:%M:%S')," : Fin de l'algorithme EM")
            obj <- scale(obj,scale=scale)
            return(obj)
          })

###############################################################################################################################################################
#Méthodes SPH

setMethod("show",
          "SPH",
          function(object) {
            cat("Objet de la classe SPH.", "\n", "\n")
            cat("Nom : Loi phase-type inhomogène",object@name, "\n", "\n")
            cat("Paramètres :", "\n","@pars", "\n", "\n")
            cat("Loi initiale", "\n", "$alpha", "\n", "\n",object@pars$alpha, "\n", "\n")
            cat("Matrice d'intensité", "\n", "$S", "\n")
            print(object@pars$S)
            cat("Paramètres de la fonction g :", "\n","@gfun", "\n", "\n")
            cat("Fonction d'intensité de type", object@gfun$name, "\n", "\n")
            cat("Paramètres de g", "\n", "$pars", "\n")
            print(object@gfun$pars)
            cat("Paramètres de la régression", "\n","@reg","\n", "$pars", "\n")
            print(object@reg$pars)
          }
)

setMethod("scale",
          "SPH",
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
          "SPH",
          function(x,y,X){
            S <- x@pars$S
            alpha <- x@pars$alpha
            s <- exit(S)
            beta <- x@gfun$pars
            theta <- x@reg$pars
            mx <- exp(X %*% theta)
            if (x@gfun$name=="gompertz"){
              c(DensiteGompertzs(alpha,S,y,beta,mx))
            }
          })

setMethod("surv",
          "SPH",
          function(x,y,X){
            S <- x@pars$S
            alpha <- x@pars$alpha
            s <- exit(S)
            beta <- x@gfun$pars
            theta <- x@reg$pars
            mx <- exp(X %*% theta)
            if (x@gfun$name=="gompertz"){
              c(SurvieGompertzs(alpha,S,y,beta,mx))
            }
          })

setMethod("haz",
          "SPH",
          function(x,y,X){
            dens(x,y,X)/surv(x,y,X)
          })