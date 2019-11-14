#########################################################################
##  SJIVE.R
##  
##  This version: 1/8/2015
##
##
##  Federico Crudu
##  Instituto de Estadística
##  Facultad de Ciencias
##  Pontificia Universidad Católica de Valparaiso
##  Avenida Errazuriz 2734
##  federico.crudu@gmail.com
##  https://sites.google.com/site/federicocrudu/
##  
##  This is an R program.
##  It will compute the SJIVE/SJEF estimator and the associated variance 
##  covariance matrix.
##  
#########################################################################


SJIVEfitF<-function(y=NULL,X1=NULL,X2=NULL,Z1=NULL,Z2 = NULL,full=2)
{
  
  X2<-Z2 # included instruments
  X<-cbind(X1,X2) # X is your nxg matrix of regressors 
  Z<- cbind(Z1,Z2) # this is your nxk matrix of instruments
  
  ZZ<-t(Z)%*%Z
  ZZ.inv<-solve(ZZ)
  
  
  Pz<-Z%*%(ZZ.inv)%*%t(Z) # projection matrix
  
  
  Dp<-diag(diag(Pz)) # diagonal matrix; the diagonal elements are the 
  # same as Pz
  
  I.n<-diag(1,length(y)) # nxn identity matrix
  
  Mz<- I.n-Pz # projection matrix
  G<-eiMu((Dp),diag(1/diag(I.n-Dp)))
  Delta <- eiMu(eiMu(Pz,G),Pz)-(eiMu(Pz,G)+eiMu(G,Pz))/2
  A<- Pz + Delta
  B<- eiMu(eiMu(Mz,G),Mz)
  C<- A-B
  
  ### minimum eigenvalue problem ##########################################

  X2X2<-eiMu(t(X2),X2)
  X2X2.inv<-solve(X2X2)
  Px2<-eiMu(eiMu(X2,(X2X2.inv)),t(X2)) # projection matrix
  
  C.star<- C-eiMu(eiMu(A,Px2),A)
  yX1<-cbind(y,X1)
  Q1<-eiMu(eiMu(t(yX1),B),yX1)
  Q2<-eiMu(eiMu(t(yX1),C.star),yX1)
  Q3<-eiMu(solve(Q1),Q2)
  
  min.lambda<- min(eigen(Q3)$values) # minimum eigenvalue
  #########################################################################
  
  
  ### SJIVE/SJEF estimator ################################################
  
  XCX<-eiMu(eiMu(t(X),C),X)
  XCy<-eiMu(eiMu(t(X),C),y)
  XBX<-eiMu(eiMu(t(X),B),X)
  XBy<-eiMu(eiMu(t(X),B),y)
  alpha<-full # Fuller parameter; alpha=0 computes the SJIVE estimator
  traceB<-sum(diag(B))
  lambda.hat<-min.lambda-alpha/traceB
  
  # beta.hat is the SJIVE/SJEF estimator
  beta.hat<-eiMu(solve(XCX-lambda.hat*XBX),(XCy-lambda.hat*XBy))
  
  #########################################################################
  
  C.hat<-C-lambda.hat*B
  XC.hatX<-eiMu(eiMu(t(X),C.hat),X)
  I.1plusg<-diag(1,(1+ncol(X))) # (1+g)x(1+g) identity matrix
  I0<-I.1plusg[,-1]
  b0<-c(1,-beta.hat)
  b1<-cbind(b0,I0)
  yX<-cbind(y,X)
  Omega.hat<-eiMu(t(yX),B)%*%yX/ncol(Z)
  Sigma.hat<-eiMu(eiMu(t(b1),Omega.hat),b1)
  sigma11.hat<-as.numeric(Sigma.hat[1,1])
  sigma12.hat<-Sigma.hat[1,2:(1+ncol(X))]
  epsilon.hat<- y-X%*%beta.hat
  
  X.tilde<- X-sigma11.hat*epsilon.hat%*%sigma12.hat
  C2<-C*C
  D.eps.hat<-diag(as.vector(epsilon.hat))
  D.eps.hat2<-eiMu(D.eps.hat,D.eps.hat)
  
  ### variance covariance estimator #######################################
  
  Meat1<-eiMu(eiMu(C,D.eps.hat2),C)
  Meat2<-eiMu(eiMu(D.eps.hat,C2),D.eps.hat)
  Meat<-Meat1+Meat2
  Bread1<-eiMu(eiMu(t(X),C.hat),X)
  Bread<-eiMu(X.tilde,solve(Bread1))
  
  # estimator of the variance covariance matrix
  Var.beta.hat<-eiMu(eiMu(t(Bread),Meat),Bread )
  return(list (beta = beta.hat,variance = Var.beta.hat))
}





