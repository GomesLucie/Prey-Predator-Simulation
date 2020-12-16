####méthode explicite d'Euler 

# a= coef de reproduction des proies
# b= coef de mortalité des proies
# c= coef de repro des predateurs
# d= coef de mortalité des prédateurs
# n= nb itérations 
# h= pas de temps 

methodEuler<-function(z,a,b,c,d,h,n){
  Z<-data.frame(t=0)
  Z<-cbind(Z,z)  #nb de proies et prédateurs à t=0
  for(i in 1:n){
    prey<-Z$prey[i]
    pred<-Z$pred[i]
    if(prey<1||pred<1){#condition d'arrêt si une population s'éteint
      return (Z)
    } 
    res<-modelLV(data.frame(prey=prey,pred=pred),a,b,c,d)
    new_prey<-prey+res$xdt*h
    new_pred<-pred+res$ydt*h
    Z<-rbind(Z,data.frame(t=i,prey=new_prey,pred=new_pred))
  }
  return(Z)
}