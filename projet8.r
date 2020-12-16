# a= coef de reproduction des proies
# b= coef de mortalité des proies
# c= coef de repro des predateurs
# d= coef de mortalité des prédateurs


###fonction équation Lotka-Volterra
  
modelLV<-function(z,a,b,c,d){
  
  xdt<-a*z$prey - b*z$prey*z$pred  #Proies
  ydt<-c*z$prey*z$pred - d*z$pred  #Prédateurs
  return(data.frame(xdt=xdt,ydt=ydt))
} 


#méthode Euler

methodEuler<-function(z,a,b,c,d,h,n){
  Z<-data.frame(t=0)
  Z<-cbind(Z,z)  #nb de proies et prédateurs à t=0
  for(i in 1:n){
    prey<-Z$prey[i]
    pred<-Z$pred[i]
    if(prey<1||pred<1){return (Z)} #condition d'arrêt si une population s'éteint
    res<-modelLV(data.frame(prey=prey,pred=pred),a,b,c,d)
    new_prey<-prey+res$xdt*h
    new_pred<-pred+res$ydt*h
    Z<-rbind(Z,data.frame(t=i,prey=new_prey,pred=new_pred))
  }
  return(Z)
}


#méthode RK4

methodRK4<-function(z,a,b,c,d,h,n){
  Z<-data.frame(t=0)
  Z<-cbind(Z,z)
  for(i in 1:n){
    if(Z$prey[i]<1||Z$pred[i]<1){return (Z)} #condition d'arrêt si une population s'éteind
    res<-data.frame(prey=Z$prey[i],pred=Z$pred[i])
    k1<-modelLV(res,a,b,c,d)
    
    res<-data.frame(prey=Z$prey[i]+h/2*k1[[1]],pred=Z$pred[i]+h/2*k1[[2]])
    k2<-modelLV(res,a,b,c,d)
    
    res<-data.frame(prey=Z$prey[i]+h/2*k2[[1]],pred=Z$pred[i]+h/2*k2[[2]])
    k3<-modelLV(res,a,b,c,d)
    
    res<-data.frame(prey=Z$prey[i]+h*k3[[1]],pred=Z$pred[i]+h*k3[[2]])
    k4<-modelLV(res,a,b,c,d)
    
    prey<-Z$prey[i]+ h/6*(k1[[1]]+2*k2[[1]]+2*k3[[1]]+k4[[1]])
    pred<-Z$pred[i]+ h/6*(k1[[2]]+2*k2[[2]]+2*k3[[2]]+k4[[2]])
    Z<-rbind(Z,data.frame(t=i,prey=prey,pred=pred))
  }
  return (Z)
}


#méthode RK4 avec chocs

methodRK4_choc<-function(z,a,b,c,d,h,n,nExp,r,nNorm,SD){
  Z<-data.frame(t=0)
  Z<-cbind(Z,z)
  for(i in 1:n){
    choc<- calculChoc(nExp,r,nNorm,SD)
    if(Z$prey[i]<1||Z$pred[i]<1){return (Z)} #condition d'arrêt si une population s'éteint
    res<-data.frame(prey=Z$prey[i],pred=Z$pred[i])
    k1<-modelLV(res,a,b,c,d)
    
    res<-data.frame(prey=Z$prey[i]+h/2*k1[[1]],pred=Z$pred[i]+h/2*k1[[2]])
    k2<-modelLV(res,a,b,c,d)
    
    res<-data.frame(prey=Z$prey[i]+h/2*k2[[1]],pred=Z$pred[i]+h/2*k2[[2]])
    k3<-modelLV(res,a,b,c,d)
    
    res<-data.frame(prey=Z$prey[i]+h*k3[[1]],pred=Z$pred[i]+h*k3[[2]])
    k4<-modelLV(res,a,b,c,d)
    
    prey<-Z$prey[i]+ h/6*(k1[[1]]+2*k2[[1]]+2*k3[[1]]+k4[[1]])-choc*Z$prey[i]
    pred<-Z$pred[i]+ h/6*(k1[[2]]+2*k2[[2]]+2*k3[[2]]+k4[[2]])
    Z<-rbind(Z,data.frame(t=i,prey=prey,pred=pred))
  }
  return (Z)
}


#calcul du choc

calculChoc<-function(nExp,r,nNorm,SD){
  choc<-0
  myT<- rexp(nExp, rate=r) #temps restant avant le prochain choc
  #plus r sera grand et moins il aura de chance de se produire (et inversement)
  if(myT>(1/r)*2+1){ 
    choc<-abs(rnorm(n=nNorm,sd=SD))
  }
  return (choc)
}


#Données utilisées pour les différentes simulations
#Partie Euler:
zo<-data.frame(prey=900,pred=800)
donnees<-methodEuler(zo,0.2,0.001,0.001,0.5,0.1,2000)
donnees[,(2:3)]<-apply(donnees[,(2:3)],2,floor) #obtenir des val entiers pour le nb proies et prédateurs
colnames(donnees)<-c("temps","proies","predateur")

#graphique de l'évolution des populations de proies et prédateurs en fonction du temps
matplot(x=donnees[,1],y=(donnees[,2:3]),type="l",xlab = "temps",ylab="Taille pop",main="Titre",col=c(1,2))
legend(x='topleft', legend=c('Prey','Predator'),col=c(1,2),lwd=2)
#graphique des orbites des solutions
plot(donnees[,2],donnees[,3],type="l",xlab = "Prey",ylab = "Predator")

#Partie RK4:
donnees<-methodRK4(zo,0.2,0.001,0.001,0.5,0.1,2000) 
donnees[,(2:3)]<-apply(donnees[,(2:3)],2,round) #obtenir des val entiers pour le nb proies et prédateurs
colnames(donnees)<-c("temps","proies","predateur")

matplot(x=donnees[,1],y=(donnees[,2:3]),type="l",xlab = "temps",ylab="Taille pop",main="Titre",col=c(1,2))
legend(x='topleft', legend=c('Prey','Predator'),col=c(1,2),lwd=2)
plot(donnees[,2],donnees[,3],type="l",xlab = "Prey",ylab = "Predator")

zo<-data.frame(prey=580,pred=120)
donnees<-methodRK4(zo,0.2,0.001,0.001,0.5,0.1,2000) 
donnees[,(2:3)]<-apply(donnees[,(2:3)],2,round) #obtenir des val entiers pour le nb proies et prédateurs
colnames(donnees)<-c("temps","proies","predateur")
points(donnees[,2],donnees[,3],type="l", col="yellow")

zo<-data.frame(prey=1000,pred=100)
donnees<-methodRK4(zo,0.2,0.001,0.001,0.5,0.1,2000)
donnees[,(2:3)]<-apply(donnees[,(2:3)],2,round) #obtenir des val entiers pour le nb proies et prédateurs
colnames(donnees)<-c("temps","proies","predateur")
points(donnees[,2],donnees[,3],type="l", col="blue")

zo<-data.frame(prey=500,pred=200) #point d'équilibre
donnees<-methodRK4(zo,0.2,0.001,0.001,0.5,0.1,2000)
donnees[,(2:3)]<-apply(donnees[,(2:3)],2,round) #obtenir des val entiers pour le nb proies et prédateurs
colnames(donnees)<-c("temps","proies","predateur")
points(donnees[,2],donnees[,3],type="l", col="red")

#partie RK4 avec chocs
zo<-data.frame(prey=900,pred=250)
donnees<-methodRK4_choc(zo,0.2,0.001,0.001,0.5,0.2,1000,1,1/7,1,0.1)
donnees[,(2:3)]<-apply(donnees[,(2:3)],2,round)
matplot(x=donnees[,1],y=(donnees[,2:3]),type="l",xlab = "temps",ylab="Taille pop",main="Evolution des populations en présence de choc sur les proies",col=c(1,2))
legend(x='topleft', legend=c('Prey','Predator'),col=c(1,2),lwd=2)