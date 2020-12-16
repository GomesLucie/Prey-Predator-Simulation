#Methode Runge-Kutta d'ordre 4
methodRK4<-function(z,a,b,c,d,h,n){
  Z<-data.frame(t=0)
  Z<-cbind(Z,z)
  for(i in 1:n){
    if(Z$prey[i]<1||Z$pred[i]<1){return (Z)} #condition d'arrêt si une population s'éteint
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