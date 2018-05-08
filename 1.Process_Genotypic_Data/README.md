##**Dissection of leaf angle variation in maize through genetic mapping and meta-analysis**

Matthew J. Dzievit, Xianran Li, Jianming Yu


###Folder organization
1. Filter - it contains the java scripts for filtering the data
2. Imputation - it contains the sliding window script
3. Other - it contains a text file of SNP types that came from the metadata from Cornell

###Genotypic data information

The raw genotypic data is located here:

Descriptions of genotypic file

1. Rows are samples, columns are SNPs

2. SNPs are numbered by position only from chr1 - chr10. So pos 1....pos 1,000,000 then starts over to pos 2 ..... 2,000,000. This is first chr1 and then chr2, etc. The true raw data includes SNPs that did not map to a chromosome, but those were removed immediately.

###Pipeline for generating SNPs
Please note: java scripts were run inside eclipse and you should change the file location of the files (and possible the name) inside the code. I did not optimize the code to run with inputting a path and filename, etc.


1. Used Tassel (http://www.maizegenetics.net/tassel) to open genotypic data file and separate into three data sets: 

	A. Parent file (tassel\_parents.txt) - pulled the parent replicates to call SNPs

	B. Progeny file for B73 population ("progeny\_B73.txt") - pulled the progeny related to the two reciprocal populations using B73 (see paper for more details)

	C. Progeny file for Mo17 population ("progeny\_Mo17.txt") - pulled the progeny related to the two reciprocal populations using Mo17 (see paper for more details)

2. Ran the "parent\_snps.java" code to call SNPs for each parent replicate. It reads in a file called "tassel\_parents.txt" and a file in a folder called "Other\snp_type.txt" It outputs a file called "parents\_final\_SNPs.txt" that contains SNPs in rows and parents in columns. This program condenses the replicates from the parents and calls a SNP for each position.

3. Next the java file "parent\_poly\_snp\_subset.java" was run with the output file from the previous step. This file separates the parents and compares the two to call SNPs (PHW30 vs B73 and PHW30 vs Mo17). It outputs a file into a folder for "B73" or "Mo17" depending on what parent it is calling for. The file name is "parents\_final\_SNPs\_filtered\_" with population name (B73 of Mo17) and ".txt" in its respective folder. 

4.  "line\_freq.java" runs the progeny filterer and genotypic coder (A/AB/B) script. It reads in the parent SNPs from previous step, the progeny SNPs for the respective population (script has a loop), and filters the SNPs based on the described criteria in the paper. It outputs a few summary files about the SNPs (line frequency, SNP frequency) and information about those files is noted in the script. The output file that goes to the next step is named "coded\_snps\_" with respective population name ("B73" or "Mo17") ".txt". Key thing is to note PHW30 is always the "A" allele and the other parent is always "B" allele.

5. "window.java" runs the sliding window algorithm as described in the paper with the output file previously listed. It returns two files, one is the genotyping error corrected file, and the other is information about the SNPs (Pos, Chr, name).