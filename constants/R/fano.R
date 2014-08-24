`fano` <-
function(p,x){
y0<- p[1]
A1<- abs(p[2])
w1<- abs(p[3])
x01<-p[4]
q<-p[5]
B1<-(2*A1/pi)*((q*w1+x-x01)^2/(4*(x-x01)^2+w1^2))
y0+B1
}

