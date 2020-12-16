###fonction �quation Lotka-Volterra

# a= coef de reproduction des proies
# b= coef de mortalit� des proies
# c= coef de repro des predateurs
# d= coef de mortalit� des pr�dateurs

modelLV<-function(z,a,b,c,d){
  
  xdt<-a*z$prey - b*z$prey*z$pred  #Proies
  ydt<-c*z$prey*z$pred - d*z$pred  #Pr�dateurs
  return(data.frame(xdt=xdt,ydt=ydt))
} 