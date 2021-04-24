###calculating the risk ratio by hand 
salm <- read.table(file = "C:\\Users\\Moon\\Desktop\\Salmonella data.csv", 
                   sep = ",", header = TRUE)


#extract counts from Cecum corresponding to each group 
A1 <- sum(salm$Cecum == 1 & salm$Group == "A")
B1 <- sum(salm$Cecum == 1 & salm$Group == "B")
C1 <- sum(salm$Cecum == 1 & salm$Group == "C")
A0 <- sum(salm$Cecum == 0 & salm$Group == "A")
B0 <- sum(salm$Cecum == 0 & salm$Group == "B")
C0 <- sum(salm$Cecum == 0 & salm$Group == "C")


#Create 2*2 contingency table 
AvsC<-array(data = c(A1, C1, A0, C0), 
            dim = c(2,2), 
            dimnames = list(First = c("Group A", "Group C"),
                            Second = c("Cecum = 1", "Cecum = 0"))) 
AvsC

BvsC<-array(data = c(B1, C1, B0, C0), 
            dim = c(2,2), 
            dimnames = list(First = c("Group B", "Group C"),
                            Second = c("Cecum = 1", "Cecum = 0"))) 
BvsC


#Function for Relative Risk and CI for Relative Risk
RR<-function(df){
  
  pi.hat.table<-df/rowSums(df)
  pi.hat.table
  
  pi.hat1<-pi.hat.table[1,1]
  pi.hat2<-pi.hat.table[2,1]
  
  # Relative risk where success is "Cecum = 1"
  RR<-round(pi.hat1/pi.hat2, 4)  #RR
  
  alpha<-0.05
  n1<-sum(df[1,])
  n2<-sum(df[2,])
  
  # Wald confidence interval
  var.log.rr<-(1-pi.hat1)/(n1*pi.hat1) + (1-pi.hat2)/(n2*pi.hat2)
  ci<-exp(log(pi.hat1/pi.hat2) + qnorm(p = c(alpha/2, 1-alpha/2)) * sqrt(var.log.rr))
  RR.CI<-round(ci, 4)  #RR CI
  
  data.frame(RR, lowerCI=RR.CI[1], upperCI=RR.CI[2])

}
 
RR(AvsC)
#The estimated probability of positive (1) for colonization if the challenge strain of SE was re-isolated 
 #is only 0.4314 times (or 43.14%) as large for the vaccinated group A than for control group C 
1-RR(AvsC)#prevented fraction and its confidence interval
#56.86% is the proportion of incidents in the control group C that could be prevented by vaccinated group A

RR(BvsC)
#The estimated probability of positive (1) for colonization if the challenge strain of SE was re-isolated 
 #is only 0.4231 times (or 42.31%) as large for the vaccinated group B than for control group C 
1-RR(BvsC)#prevented fraction and its confidence interval
#57.69% is the proportion of incidents in the control group C that could be prevented by vaccinated group B


###calculating the risk ratio with riskratio()
#install.packages('fmsb')
library(fmsb)
fmsb::riskratio(11,25,51,50, conf.level=0.95, p.calc.by.independence=TRUE)#same as AvsC

