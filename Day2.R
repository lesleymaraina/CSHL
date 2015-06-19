install.packages("RSQLite")
install.packages("RSQLite")
source('http://bioconductor.org/biocLite.R')
library('org.Hs.eg.db')
install.packages("DBI")
source('http://bioconductor.org/biocLite.R')
biocLite('org.Hs.eg.db')
library('org.Hs.eg.db')
human = org.Hs.eg.db
columns(human)
keytypes(human)
head(keys(human,keytype='ENTREZID'))
head(keys(human,keytype='PFAM'))
head(keys(human,keytype='SYMBOL'))
head(keys(human,keytype='SYMBOL',pattern='AIC'))#pulls out all of the keys that have the pattern AIC
head(keys(human,keytype='SYMBOL',pattern='AIC',column='ENTREZID'))#search symbol key for pattern AIC, and return the gene name for each row
select(human,keys=c('672','1'),keytype='ENTRselect(human,keys=c('672','1'),keytype='ENTREZID',columns=c('ONTOLOGY'))
#CC: cellular compartment; BP: biological process; MF: molecular function
select(human,keys=c('672','1'),keytype='ENTREZID',columns=c('GO','SYMBOL'))
a = select(human,keys=c('672','1'),keytype='ENTREZID',columns=c('GO','SYMBOL'))
class(a)
a = select(human,keys=keys(human,'ENTREZID'),keytype='ENTREZID',columns=c('GO','SYMBOL'))#retrieves every key from every row in the human database
dim(a)
b = mapIds(human,keys=c('672','1'),keytype='ENTREZID',column='GO',multiVals='list')#specify keytype and output types[columns], if there are multiple values put them in a list
class(b)
D = mapIds(human,keys=b[[1]],keytype='GO',column='SYMBOLS',multiVals='list')
lapply(d,length)
a = sample(1:24,12)#sample any number between 1 and 24, and give 12 numbers
a
b = sample(1:24,12)#sample any number between 1 and 24, and give 12 numbers
intersect(a,b)
union(a,b)
length(union(a,b))
a %in% b #determine if a is in b
a[a %in% b] #determine if a is in b
union(a,b) %in% 1:24
1:24 %in% union(a,b)
table(1:24 %in% union(a,b))#takes a set of values (T/F, factors, etc) and counts these values
EZID',columns=c('SYMBOL','REFSEQ'))#select from human
select(human,keys=c('672','1'),keytype='ENTREZID',columns=c('SYMBOL','REFSEQ'))#select from human#672=BRCA1
HUMAN
human
human #shows version information
