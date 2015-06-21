#Day 3 Rafael
femaleMiceWeights <- read.csv("~/Desktop/CSHL Course Statistics Course/Day3/femaleMiceWeights.csv")
View(femaleMiceWeights)
dat <- read.csv("~/Desktop/CSHL Course Statistics Course/Day3/femaleMiceWeights.csv")
View(femaleMiceWeights)
control &lt;- dat[1:12,2]
control <- dat[1:12,2]
treatment <- dat[13:24,2]
print(mean(treatment))
print(mean(control))
diff <- mean(treatment)-mean(control)
print(diff)
femaleControlsPopulation <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/femaleControlsPopulation.csv", quote="\"")
View(femaleControlsPopulation)
dat2 <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/femaleControlsPopulation.csv", quote="\"")
View(femaleControlsPopulation)
dat2 <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/femaleControlsPopulation.csv", quote="\"")
View(femaleControlsPopulation)
dat2 <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/femaleControlsPopulation.csv", quote="\"")
View(femaleControlsPopulation)
mice <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/mice.csv", header=TRUE, quote="\"")
View(mice)
n <- 100000
null <- vector("numeric",n)
for(i in 1:n){
control <- sample(population[,1],12)
treatment <- sample(population[,1],12)
null[i] <- mean(treatment) - mean(control)
}
population <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/mice.csv", header=TRUE, quote="\"")
View(mice)
n <- 100000
null <- vector("numeric",n)
for(i in 1:n){
control <- sample(population[,1],12)
for(i in 1:n){
for(i in 1:n){
for(i in 1:n){
control <- sample(population[,1],12)
treatment <- sample(population[,1],12)
null[i] <- mean(treatment) - mean(control)
}
warnings()
population
population <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/mice.csv", header=TRUE, quote="\"")
View(mice)
population <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/mice.csv", header=TRUE, quote="\"")
View(mice)
mice2 <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/mice2.csv", header=TRUE, quote="\"")
View(mice2)
population <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/mice2.csv", header=TRUE, quote="\"")
View(mice2)
n <- 100000
null <- vector("numeric",n)
for(i in 1:n){
control <- sample(population[,1],12)
treatment <- sample(population[,1],12)
null[i] <- mean(treatment) - mean(control)
}
warnings()
pnrom(2)
pnorm(2)
mean(null>diff)
1-pnorm(diff,mean(null,sd(null)))
mice3 <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/mice3.csv", header=TRUE, quote="\"")
View(mice3)
mice3 <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/mice3.csv", header=TRUE, quote="\"")
View(mice3)
mice_pheno <- read.csv("~/Desktop/CSHL Course Statistics Course/Day3/mice_pheno.csv")
View(mice_pheno)
dat <- read.csv("~/Desktop/CSHL Course Statistics Course/Day3/mice_pheno.csv")
View(mice_pheno)
controlPopulation <- dat[dat$Sex == "F" & dat$Diet == "chow", 3]
length(controlPopulation)
hfPopulation <- dat[dat$Sex == "F" & dat$Diet == "hf", 3]
length(hfPopulation)
mean(controlPopulation)
mean(hfPopulation)
mean(hfPopulation) - mean(controlPopulation)
mean(control)
mean(treatment)
mean(treatment);12,3
x=dat[1:12,3]
babies <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/babies.txt", header=TRUE, quote="\"")
View(babies)
bwt.nonsmoke = babies$bwt[babies$smoke==0]
bwt.smoke = babies$bwt[babies$smoke==1]
mean(bwt.nonsmoke)-mean(bwt.smoke)
X.ns = mean(dat.ns)
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- tempfile()
download(url,destfile=filename)
dat <- read.csv(filename)
head(dat)
hfPopulation <- dat[dat$Sex=="F" & dat$Diet=="hf",3]
controlPopulation <- dat[dat$Sex=="F" & dat$Diet=="chow",3]
mu_hf <- mean(hfPopulation)
mu_control <- mean(controlPopulation)
print(mu_hf - mu_control)
x<-controlPopulation
N<-length(x)
popvar <- mean((x-mean(x))^2)
identical(var(x),popvar)
identical(var(x)*(N-1)/N, popvar)
popvar <- function(x) mean( (x-mean(x))^2)
popsd <- function(x) sqrt(popvar(x))
N <- 12
hf <- sample(hfPopulation,N)
control <- sample(controlPopulation,N)
t.test(hf,control)
t.test(hf,control)
N <- 50
hf <- sample(hfPopulation,N)
control <- sample(controlPopulation,N)
t.test(hf,control)
t.test(hf,control)$p.vslur
t.test(hf,control)$p.value
Ns <- c(3,12,25,50)
B <- 10000 #number of simulations
res <-  sapply(Ns,function(n){
replicate(B,mean(sample(hfPopulation,n))-mean(sample(controlPopulation,n)))
})
n=6
B <- 10000 #number of simulations
res <-  sapply(Ns,function(n){
replicate(B,mean(sample(hfPopulation,n))-mean(sample(controlPopulation,n)))
})#random variable: takes a random sample diff between means
hist(null)
qnorm(.75)
qqnorm(rnorm(1000))
qqline(rnorm(1000))
null = replicate(B,mean(sample(hfPopulation,n))-mean(sample(controlPopulation,n)))
null = replicate(B,t.test(sample(hfPopulation,n))-mean(sample(controlPopulation,n))$stat)
null = replicate(B,t.test(sample(hfPopulation,n))-sample(controlPopulation,n)$stat)
babies <- read.table("~/Desktop/CSHL Course Statistics Course/Day3/babies.txt", header=TRUE, quote="\"")
View(babies)
baby <- babies[,order(smoke)]
baby <- babies[,order('smoke')]
baby
view(baby)
View(baby)
baby <- babies[order('smoke'),]
View(baby)
baby <- babies[,order('smoke')]
View(baby)
setorder(babies, smoke, na.last=FALSE)
baby = setorder(babies, smoke, na.last=FALSE)
temp = babies[order('smoke'),]
head(temp)
temp = babies[order(babies$smoke),]
head(temp)
tail(temp)
order(babies$smoke)
order(-babies$smoke)
temp = babies[c(3, 5, 6, 7, 19),]
temp
temp = babies[(babies$smoke==0),]
table(temp$smoke)
View(temp)
temp = babies[(babies$smoke==0),]#remove the rows that have a smoke value of "0"
temp = babies[order(babies$smoke),]#order the table according to the column smoke
order(-babies$smoke)#reverse order the table by smore
order(-babies$smoke)#reverse order the table by smoke
