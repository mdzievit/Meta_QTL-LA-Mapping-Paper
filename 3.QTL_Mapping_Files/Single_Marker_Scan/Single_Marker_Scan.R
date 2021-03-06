library(tidyverse)
library(qqman)

print(sessionInfo())

#######Single marker scan for B73 populations ############

##Read coded SNP data in before it sliding window
data <- read_tsv("B73-Genotypes.txt",
                 col_names = FALSE,
                 na = "N") %>%
  t()
colnames(data) <- data[1,]
data <- data[-1,]
##These are the snp positions formatted for manhattan plot
snps <- read_tsv("B73-SNPs.txt")

##Read phenotype data in. Contains the F2 and F2:3 data
pheno <- read_tsv("B73-Phenotypes.txt")

##Bind the genotype and phenotype data together
data <- cbind(pheno[,2:3],data)

##Set up the f2 and f3 pvalue vector for loop
len <- dim(data)[2]
f2_pvalues = c(1:(len - 3))
f3_pvalues = c(1:(len - 3))

##Run loop that performs an f-test on each SNP
for (i in 4:ncol(data) )
{
  test = aov(data$F2_MLA ~ data[,i], )
  holder = summary(test)[[1]][["Pr(>F)"]]
  f2_pvalues[i - 3] = holder[1]
  
  test = aov(data$F3_MLA ~ data[,i])
  holder = summary(test)[[1]][["Pr(>F)"]]
  f3_pvalues[i - 3] = holder[1]
  
}

##Binds SNP information to the pvalues
snps_f2_b73 = cbind(snps,f2_pvalues)
colnames(snps_f2_b73)[4] = "P"

snps_f3_b73 = cbind(snps,f3_pvalues)
colnames(snps_f3_b73)[4] = "P"

manhattan(snps_f2_b73)

#######Single marker scan for Mo17 population ############

##Read coded SNP data in before it sliding window
data <- read_tsv("Mo17-Genotypes.txt",
                 col_names = FALSE,
                 na = "N") %>%
  t()
colnames(data) <- data[1,]
data <- data[-1,]

##These are the snp positions formatted for manhattan plot
snps = read_tsv("Mo17-SNPs.txt")

##Read phenotype data in. Contains the F2 and F2:3 data
pheno <- read_tsv("Mo17-Phenotypes.txt")

##Bind the genotype and phenotype data together
data <- cbind(pheno[,2:3],data)

##Set up the f2 and f3 pvalue vector for loop
len <- dim(data)[2]
f2_pvalues = c(1:(len - 3))
f3_pvalues = c(1:(len - 3))

##Run loop that performs an f-test on each SNP
for (i in 4:ncol(data) )
{
  test = aov(data$F2_MLA ~ data[,i], )
  holder = summary(test)[[1]][["Pr(>F)"]]
  f2_pvalues[i - 3] = holder[1]
  
  test = aov(data$F3_MLA ~ data[,i])
  holder = summary(test)[[1]][["Pr(>F)"]]
  f3_pvalues[i-3] = holder[1]
  
}

##Binds SNP information to the pvalues
snps_f2_mo17 = cbind(snps,f2_pvalues)
colnames(snps_f2_mo17)[4] = "P"

snps_f3_mo17 = cbind(snps,f3_pvalues)
colnames(snps_f3_mo17)[4] = "P"

manhattan(snps_f2_mo17)
manhattan(snps_f3_mo17)
##Output the datafiles

write_tsv(snps_f2_b73 %>%
            rename(F2 = P) %>%
            left_join(snps_f3_b73 %>%
                        rename(F3 = P),
                      by = c("SNP","CHR","BP")),
          path = "pvalues_B73.txt")
write_tsv(snps_f2_mo17 %>%
            rename(F2 = P) %>%
            left_join(snps_f3_mo17 %>%
                        rename(F3 = P),
                      by = c("SNP","CHR","BP")),
          path = "pvalues_Mo17.txt")

