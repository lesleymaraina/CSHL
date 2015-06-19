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
select(human,keys=c('672','1'),keytype='ENTREZID',columns=c('SYMBOL','REFSEQ'))#select from human
select(human,keys=c('672','1'),keytype='ENTREZID',columns=c('SYMBOL','REFSEQ'))#select from human#672=BRCA1
HUMAN
human
human #shows version information
