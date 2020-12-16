###fonction équation Lotka-Volterra

# a= coef de reproduction des proies
# b= coef de mortalité des proies
# c= coef de repro des predateurs
# d= coef de mortalité des prédateurs

modelLV<-function(z,a,b,c,d){
  
  xdt<-a*z$prey - b*z$prey*z$pred  #Proies
  ydt<-c*z$prey*z$pred - d*z$pred  #Prédateurs
  return(data.frame(xdt=xdt,ydt=ydt))
} 