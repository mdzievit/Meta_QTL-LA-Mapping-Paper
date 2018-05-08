##**Dissection of leaf angle variation in maize through genetic mapping and meta-analysis**

Matthew J. Dzievit, Xianran Li, Jianming Yu

##QTL Mapping Files

###Organization of folders

1. Consensus Map - it contains the linkage map files from the IcIM software output and R script to create the consensus linkage map
2. Rest of the contents in this folder were used for genetic mapping using the IcIM software


###Outline of methods
1. Imported the PopName"\_Bin\_File.xlsx" into the software and followed manual and options depicted in paper
2. Took the map file from QTL IcIM file and combined the two linkage maps into a file "Consensus\_Map\Genetic\_Maps.txt" and  then ran the "consensus_map.r" script to generate the "ConsensusMap.txt" file. 
3. Formatted excel file to reimport into IcIM software to run the genetic mapping
4. Output the LOD score files and combined them into "All_Results.txt" file
5. Ran the "Single\_Marker\_Scan\Single\_Marker\_Scan.R" script to generate the pvalue results
6. Plotted the results in "QTL\_Maps\_Code.r" using the single marker scan results and the LOD results to generate figure 2 and supplemental figure 2 from thepaper





