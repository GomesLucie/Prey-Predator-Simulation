#Calcul du choc à utiliser avec methodRK4_choc

calculChoc<-function(nExp,r,nNorm,SD){
  choc<-0
  myT<- rexp(nExp, rate=r) #temps restant avant le prochain choc
  #plus r sera grand et moins il aura de chance de se produire (et inversement)
  if(myT>(1/r)*2+1){ 
    choc<-abs(rnorm(n=nNorm,sd=SD))
  }
  return (choc)
}
#r=1/2 -> 5
#r=1/3 -> 7 ou 8
#r=1/5 -> 12 ou 13