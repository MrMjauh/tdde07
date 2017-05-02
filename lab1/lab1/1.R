inc = c(14,25,45,25,30,33,19,50,34,67)
mmu = 3.5



p <- function(y,mmu,sigmasq){
  return(1/(y*sqrt(2*pi*sigmasq))*exp(-1/(2*sigmasq)*(log(y)-mmu)^2))
}
tausq = sum((log(inc)-mmu)^2)/length(inc)

#2a
#install.packages("geoR")
sim = geoR::rinvchisq(10000,length(inc),tausq)

x = c()
y = c()
for(i in seq(0.01,1, by = 0.0100)){
  x = c(x,i)
  y = c(y,c(geoR::dinvchisq(i,length(inc),tausq)))
}

h = hist(sim,breaks=100)

h$density = h$counts/max(h$counts)

plot(h, freq=FALSE)
lines(x,y/max(y),col='blue')

#2b

G <- function(x){
  return(2*pnorm(x/sqrt(2))-1)
}

h = hist(G(sqrt(sim)),breaks=100)


v = density(G(sqrt(sim)))
plot(v)
#2c

a1 = quantile(G(sqrt(sim)),2.5/100)
b1 = quantile(G(sqrt(sim)),97.5/100)

dn = cumsum(v$y)/sum(v$y)
li = which(dn>=0.025)[1]
ui = which(dn>=0.975)[1]
a2 = v$x[li]
b2 = v$x[ui]


z = HDInterval::hdi(v)


