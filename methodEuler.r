####m�thode explicite d'Euler 

# a= coef de reproduction des proies
# b= coef de mortalit� des proies
# c= coef de repro des predateurs
# d= coef de mortalit� des pr�dateurs
# n= nb it�rations 
# h= pas de temps 

methodEuler<-function(z,a,b,c,d,h,n){
  Z<-data.frame(t=0)
  Z<-cbind(Z,z)  #nb de proies et pr�dateurs � t=0
  for(i in 1:n){
    prey<-Z$prey[i]
    pred<-Z$pred[i]
    if(prey<1||pred<1){#condition d'arr�t si une population s'�teint
      return (Z)
    } 
    res<-modelLV(data.frame(prey=prey,pred=pred),a,b,c,d)
    new_prey<-prey+res$xdt*h
    new_pred<-pred+res$ydt*h
    Z<-rbind(Z,data.frame(t=i,prey=new_prey,pred=new_pred))
  }
  return(Z)
}