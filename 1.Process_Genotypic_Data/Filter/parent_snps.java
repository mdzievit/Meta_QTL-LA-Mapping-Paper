package Filter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class parent_snps {

	public static void main(String[] args) throws IOException {
			
		// This is the IUPAC SNPs that we are going to be looking for to create a count
		String snpFile = "C:\\Users\\mdzievit.IASTATE\\Desktop\\F3_SNPs\\Other\\snp_type.txt";
		 
		//This is going to load the types of SNPs into a an arary
		ArrayList<String> snpType = new ArrayList();
		String str = "";
		BufferedReader read = new BufferedReader(new FileReader(snpFile));
		while ((str = read.readLine()) != null)
	 	{
			snpType.add(str);
	 	}
		
		//This is going to setup the matrix that holds the parent SNP data, it does not include SNP headers
		String mainFile = "C:\\Users\\mdzievit.IASTATE\\Desktop\\F3_SNPs\\tassel_parents.txt";
		BufferedReader read2 = new BufferedReader( new FileReader(mainFile));
		String header = read2.readLine();
		String [] snpHeader = header.split("\t");
		int snpLength = snpHeader.length;
		int snpTypeLength = snpType.size();
		int mb_pass = 4;
		int ph_pass = 6;
	
		//This is where the counts of the SNPs go
		int [][] B73snp = new int [snpTypeLength][snpLength];
		int [][] Mo17snp = new int [snpTypeLength][snpLength];
		int [][] PHW30snp = new int [snpTypeLength][snpLength];
		String [] holder = new String [snpLength];

		String finder = ":";
		String line = "";
		int end = 0;
		while ((str = read2.readLine()) != null)
 		{
			holder = str.split("\t");
			end = str.indexOf(finder) + 1;
			line = str.substring(0,end);
			if (line.equals("B73:"))
			{
				for (int j = 1; j <snpLength; j++)
				{
					for (int i = 0; i < snpTypeLength; i++)
					{
						if(snpType.get(i).equals(holder[j]))
						{
							B73snp[i][j] = B73snp[i][j] + 1;
						}
					}
				}
			}
			if (line.equals("Mo17:"))
			{
				for (int j = 1; j <snpLength; j++)
				{
					for (int i = 0; i < snpTypeLength; i++)
					{
						if(snpType.get(i).equals(holder[j]))
						{
							Mo17snp[i][j] = Mo17snp[i][j] + 1;
						}
					}
				}
			}
			if (line.equals("PHW30:"))
			{
				for (int j = 1; j <snpLength; j++)
				{
					for (int i = 0; i < snpTypeLength; i++)
					{
						if(snpType.get(i).equals(holder[j]))
						{
							PHW30snp[i][j] = PHW30snp[i][j] + 1;
						}
					}
				}
			}
 		}	
		
		//This is where the final SNP for the parent is going to go
		String [] snpCallB73 = new String[snpLength];
		String [] snpCallMo17 = new String[snpLength];
		String [] snpCallPHW30 = new String[snpLength];
		

	
		for (int j = 1; j <snpLength; j++)
		{
			int finalSNP1 = 0;
			int finalSNP2 = 0;
			int finalSNP3 = 0;
			int max1 = 0;
			int max2 = 0;
			int max3 = 0;
			for (int i = 0; i < snpTypeLength; i++)
			{
				if (B73snp[i][j] > max1)
				{
					max1 = B73snp[i][j];
					finalSNP1 = i;
				}
				if (Mo17snp[i][j] > max2)
				{
					max2 = Mo17snp[i][j];
					finalSNP2 = i;
				}
				if (PHW30snp[i][j] > max3)
				{
					max3 = PHW30snp[i][j];
					finalSNP3 = i;
				}
			}
			if (max1 >= mb_pass)
			{
				snpCallB73[j] = snpType.get(finalSNP1);
			}
			else
			{
				snpCallB73[j] = "NA";
			}
			if (max2 >= mb_pass)
			{
				snpCallMo17[j] = snpType.get(finalSNP2);
			}
			else
			{
				snpCallMo17[j] = "NA";
			}
			if (max3 >= ph_pass)
			{
				snpCallPHW30[j] = snpType.get(finalSNP3);
			}
			else
			{
				snpCallPHW30[j] = "NA";
			}
		}
		

		//This is going output the SNPs into a text file. Parents will be in a column, SNPs by rows
		
		String outputFilename = "C:\\Users\\mdzievit.IASTATE\\Desktop\\F3_SNPs\\Test\\parents_final.txt";
		BufferedWriter out = new BufferedWriter (new FileWriter(outputFilename,true));
		out.write("Position" + "\t" + "B73" + "\t" + "Mo17" + "\t" + "PWH30" + "\t" + "\n");
		
		for (int i = 1; i <snpLength ; i++)
		{
			out.write( snpHeader[i] + "\t" + snpCallB73[i] + "\t" + snpCallMo17[i] + "\t" + snpCallPHW30[i] +"\t" +"\n");
		}
		out.close();
		System.out.println(snpLength-1);
		System.out.println("All Done I Hope!!");	
		
		
	}

}
