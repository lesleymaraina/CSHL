table(greaterThan1)
r10K = (r10k - mean(r10k))/sd(r10k)
r10K.m1 = (r10K.m1 - mean(r10K.m1))/sd(r10F.m1)
r10K.m1 = (r10K.m1 - mean(r10K.m1))/sd(r10K.m1)
r10K.s2 = (r10K.s2 - mean(r10K.s2))/sd(r10K.s2)
mean(z10K)
mean(r10K)
z10K = (r10k - mean(r10k))/sd(r10k)
z10K.s2 = (r10K.s2 - mean(r10K.s2))/sd(r10K.s2)
z10K.m1 = (r10K.m1 - mean(r10K.m1))/sd(r10K.m1)
mean(z10K.s2)
p = pnorm(z10k)
p = pnorm(z10K)
hist(p)
load("derisi.rda")
load.package('rda')
pnorm(c(-3,-2,-1,0,1,2,3),lower.tail = T)
pnorm(c(-3,-2,-1,0,1,2,3),lower.tail = F)#find the probability that a variable is greater than a given number
pnorm(abs(-3:3),lower.tail = F)
pnorm(abs(-3:3),lower.tail = F)#lower.tail=F only gives values that range from 0 to 0.5
plot(seq(-4,4,.1), pnorm(abs(seq(-4,4,.1)),lower.tail=F), pch = 16,ylab="p-value",xlab="x",type="l")
plot(seq(-4,4,.1), pnorm(abs(seq(-4,4,.1)),lower.tail=F), pch = 16,ylab="p-value",xlab="x",type="l")#define?
set.seed(1)
r10K = rnorm(10000)
hist(abs(r10K),xlim=c(-4,4),col='red',breaks=20,main="Histogram",xlab="")
hist(r10K,col="blue",add=T,breaks=40)
legend(x = -4, y = 1500,legend = c("r10K","abs(r10K)"),fill=c("blue","red"))
hist(pnorm(abs(r10K),lower.tail=F),breaks = 10,xlim=c(0,1),main="",xlab="P-value")
hist(2*pnorm(abs(r10K),lower.tail=F),breaks = 20,add=T,col="grey")
legend("topright",legend = c("pnorm(abs(r10K), lower.tail = F)","2 * pnorm(abs(r10K), lower.tail = F)"),fill=c("white","grey"))
load("/Users/lesleychapman/Downloads/derisi.rda")
head(derisi.rda)
apple = load("/Users/lesleychapman/Downloads/derisi.rda")
head(appl)
head(apple)
head(derisi)
colnames(derisi)
total = derisi[,c(3:9,17:23)]
head (total)
bkg = derisi[,c(10:16,24:30)]
head(bkg)
intesity = total-bkg
head(intesity)
ratios = intensity[,1:7]/intensity[,8:14] #experimental/control = RX/GX
head(ratios)
intensity
intensity = total-bkg
ratios = intensity[,1:7]/intensity[,8:14] #experimental/control = RX/GX
head(ratios)
ratios = intensity[,1:7]/intensity[,8:14] #gives intesity of expression of each gene experimental/control = RX/GX
head(ratios)
colnames(intesity)
plot.ecdf(intensity,R1)
logratios = log2(ratios)
rownames(logratios) = derisi$ORF
head(logratios)
cols = rainbow(7)
logratios = log2(ratios)
rownames(logratios) = derisi$ORF#rownames(logratios)->gets the orf gene naems from the derisi data set and are used as a new row name for the logratios matrix; the rownames(X) adds column names to a matrix and must be the same the length as the matrix
head(logratios)
cols = rainbow(7)
plot.ecdf(logratios[,1], col=cols[1],lwd=3)
plot.ecdf(logratios[,2], col=cols[2], add=T,lwd=3)
plot.ecdf(logratios[,3], col=cols[3], add=T,lwd=3)
plot.ecdf(logratios[,4], col=cols[4], add=T,lwd=3)
plot.ecdf(logratios[,5], col=cols[5], add=T,lwd=3)
plot.ecdf(logratios[,6], col=cols[6], add=T,lwd=3)
plot.ecdf(logratios[,7], col=cols[7], add=T,lwd=3)
legend("topleft",legend = paste0(seq(9,21,2),"hrs"),lwd=3,col=cols)
sd.null = sd(logratios[,1])
summary(sd.null)
sd.null2 = sd(logratios[,6])
summary(sd.null2)
z1 = (logratios[,1]-mean(logrations[,1]))/sd.null
z1 = (logratios[,1]-mean(logratios[,1]))/sd.null
mean(z1)
variance(z1)
var(z1)
z2 = (logratios[,7]-mean(logratios[,1]))/sd.null
mean(z2)
sd(z2)
z2 = (logratios[,7]-mean(logratios[,7]))/sd.null
mean(z2)
sd(z2)
plot.ecdf(z7,lwd=3,col="red")
plot.ecdf(z1,add=T,lwd=3)
legend("topleft",legend = c("z1","z7"),lwd=3,col=c("black","red"))
plot.ecdf(z2,lwd=3,col="red")
plot.ecdf(z1,add=T,lwd=3)
legend("topleft",legend = c("z1","z2"),lwd=3,col=c("black","red"))
p1 = pnorm(z1,lower.tail=F)
hist(p1)
plot.ecdf(p1)
p7 = pnorm(z2,lower.tail=F)
p7
hist(p7)
plot.ecdf(p7)
rownames(logratios)[which.min(p7)]
rownames(logratios)[which.max(abs(logratios[,7]))]
obs = sum(p7<0.001)
obs
exp = .001*nrow(derisi)
exp
FDR = exp/obs
FDR
obs2 = sum(p7<0.00001)
obs2
exp2 = .00001*nrow(derisi)
exp2
FDR2 = exp2/obs2
FDR2
z=apply(logratios, MARGIN = 2, FUN = function(x) (x-mean(x))/sd.null)
save(z,file="z.rda")
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
temp = babies[(babies$smoke==0),]#rearrange the rows in a table accordin to the column value 'smoke'
clear
load("/Users/lesleychapman/Downloads/z.rda")
source("http://bioconductor.org/biocLite.R")
biocLite("org.Sc.sgd.db")
load("/Users/lesleychapman/Downloads/derisi-2.rda")
derisi <- load("/Users/lesleychapman/Downloads/derisi-2.rda")
View(derisi)
load("derisi.rda")
head(derisi)
load("derisi.rda")
head(derisi)
derisi.orfs <- derisi$ORF
load.packages(GO.db)
load.packages('GO.db')
library(GO.db)
library(org.Sc.sgd.db)
install.packages(GO.db)
install.packages('GO.db')
source("http://bioconductor.org/biocLite.R")
biocLite("GO.db")
source("http://bioconductor.org/biocLite.R")
biocLite("GO.db")
hclust(d)
source("http://www.bioconductor.org/biocLite.R")
biocLite("genefilter")
install.packages("devtools")
library(devtools)
install_github("ririzarr/rafalib")
load("GSE5859Subset.rda")
load("/Users/lesleychapman/Downloads/tissuesGeneExpression-2.rda")
load("/Users/lesleychapman/Downloads/tissuesGeneExpression.rda")
load("/Users/lesleychapman/Downloads/GSE5859Subset.rda")
dim(geneExpression)
dim(sampleInfo)
head(sampleInfo)
sampleInfo$group
match(sampleInfo$filename,colnames(geneExpression))
heatmap(e[1:100,])
heatmap(e[1:100,])#heatmap of the first 100 elements in matrix e
library(rafalib); mypar2(1,1) ##optional
install.packages(RColorBrewer)
install.packages('RColorBrewer')
install.packages("RColorBrewer")
library(rafalib); mypar2(1,1) ##optional
library(rafalib); mypar2(1, 1) ##optional
myplclust(hclust(d), labels=as.factor(tissue),
lab.col = as.numeric(as.factor(tissue)))
d <- dist(t(e)) #transpose every element in the matrix -> turns every row into a column before comparing
x=cutree(h, h=9)
library(rafalib); mypar2(1, 1) ##optional
myplclust(hclust(d), labels=as.factor(tissue),
lab.col = as.numeric(as.factor(tissue)))
suppressPackageStartupMessages(library(MASS))
n = 100
mypar2(1,1)
set.seed(1)
y=t(mvrnorm(n,c(0,0), matrix(c(1,0.95,0.95,1),2,2)))
plot(y[1,], y[2,], xlab="Twin 1 (standardized height)",
ylab="Twin 2 (standardized height)", xlim=c(-3,3), ylim=c(-3,3))
points(y[1,1:2], y[2,1:2], col=2, pch=16)
library(rafalib)
mypar2(1,1)
ftissue = as.factor(tissue)
plot(pc$rotation[,1], pc$rotation[,2], bg=as.numeric(ftissue), pch=21,
xlab="First dimension", ylab="Second dimension")
legend("topleft", levels(ftissue), col=seq(along=levels(ftissue)),
pch=15, cex=.5)
pc <- prcomp( e - rowMeans(e))
plot(pc$rotation[,1], pc$rotation[,2], bg=as.numeric(ftissue), pch=21,
xlab="First dimension", ylab="Second dimension")
legend("topleft", levels(ftissue), col=seq(along=levels(ftissue)),
pch=15, cex=.5)
ve = pc$sdev
plot(ve^2 / sum(ve^2)*100, xlab="PC", ylab="Variance explained")#proportion of variance compared to the total
source("http://bioconductor.org/biocLite.R")
biocLite(c("AnnotationHub", "Homo.sapiens",
"TxDb.Hsapiens.UCSC.hg19.knownGene",
"BSgenome.Hsapiens.UCSC.hg19", "biomaRt",
"TxDb.Athaliana.BioMart.plantsmart22"))
install.packages(‘BSgenome.Hsapiens.UCSC.hg19’ )
install.packages(‘BSgenome.Hsapiens.UCSC.hg19’)
source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
load("AnnotationHub")
load(AnnotationHub)
library('AnnotationHub')
library('Homo.sapiens')
library('biomaRt')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("TxDb.Athaliana.BioMart.plantsmart2")
ah <- AnnotationHub()
library('AnnotationHub')
search()
LS(2)
ls(2)
ah = AnnotationHub()#makes an annotation hub object
ah = AnnotationHub()#makes an annotation hub object
search()
ls("tools:rstudio")
colnames(mcols(ah))
hubCache()
length(res)
ah <- AnnotationHub()
hubCache()
ah <- AnnotationHub()
ah <- AnnotationHub()
