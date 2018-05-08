package Filter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class line_freq {

	public static void main(String[] args) throws IOException {
		
		/* 
		 
		 This filters the progeny SNPs, for missing data and HWE, and then codes them, A, B, AB. PHW30 is always A, and B73 or Mo17 is always B
		
		*/
		
		String pops [] = new String []{"B73","Mo17"};
		for (int pop = 0; pop<pops.length ; pop++) 
		{	 		
			System.out.println("Start");
			//***********These are all of the input files:
			
			//These are the progeny snps that we are going to process
			String snpFile = "C:\\Users\\mdzievit.IASTATE\\Box Sync\\Desktop\\F3_SNPs\\" + pops[pop] + "\\progeny_" + pops[pop] + ".txt";
			BufferedReader read = new BufferedReader(new FileReader(snpFile));
			String str;
			str = read.readLine();
			String header [] = str.split("\t");
			ArrayList<String>finalParentSNPs = new ArrayList<String>();
		
			//This is the parental file that we are processing.
			String parentFile = "C:\\Users\\mdzievit.IASTATE\\Box Sync\\Desktop\\F3_SNPs\\" + pops[pop] + "\\parents_final_SNPs_filtered_" + pops[pop] + ".txt";
			BufferedReader parentRead = new BufferedReader (new FileReader(parentFile));
			String str2;
			str2  = parentRead.readLine();
			
			//This is going to output the coded file that is filtered
			String outputFilename = "C:\\Users\\mdzievit.IASTATE\\Box Sync\\Desktop\\F3_SNPs\\" + pops[pop] + "\\Less\\" + "coded_snps_" + pops[pop] + ".txt";
			BufferedWriter out = new BufferedWriter (new FileWriter(outputFilename,true));
			
			//This is going to output all the lines and a summary for each one.
			String lineFilename = "C:\\Users\\mdzievit.IASTATE\\Box Sync\\Desktop\\F3_SNPs\\" + pops[pop] + "\\Less\\" + "line_frequency_" + pops[pop] + ".txt";
			BufferedWriter lineOut = new BufferedWriter (new FileWriter(lineFilename,true));
			
			//This outputs all the SNPs and a summary of the calls
		 	String snpFilename = "C:\\Users\\mdzievit.IASTATE\\Box Sync\\Desktop\\F3_SNPs\\" + pops[pop] + "\\Less\\" + "snp_frequency_" + pops[pop] + ".txt";
			BufferedWriter snpOut = new BufferedWriter (new FileWriter(snpFilename,true));
			
		 	//This outputs the filtered SNP positions
		 	String finalSNPsfilename = "C:\\Users\\mdzievit.IASTATE\\Box Sync\\Desktop\\F3_SNPs\\" + pops[pop] + "\\Less\\" + "finalSNP_positions_" + pops[pop] + ".txt";
			BufferedWriter finalSNPs = new BufferedWriter (new FileWriter(finalSNPsfilename,true));
			
		 	String transSNPsfilename = "C:\\Users\\mdzievit.IASTATE\\Box Sync\\Desktop\\F3_SNPs\\" + pops[pop] + "\\Less\\" + "transposed_snps_" + pops[pop] + ".txt";
			BufferedWriter transSNPs = new BufferedWriter (new FileWriter(transSNPsfilename,true));
					
			//Matrix of the parental SNPs
			ArrayList<ArrayList<String>> parentSNPs = new ArrayList<ArrayList<String>>();

			
			while(( str2 = parentRead.readLine()) != null)
			{
				//This holds the line for process
				String holder [] = str2.split("\t");
				//This holds the line of parental SNPs to put into the matrix
				ArrayList<String> parents = new ArrayList<String>();
				if(!pops[pop].equals("Both"))
				{
					parents.add(holder[0]);
					parents.add(holder[pop+1]);
					parents.add(holder[3]);
					parentSNPs.add(parents);
				}
				else
				{
					parents.add(holder[0]);
					parents.add(holder[1]);
					parents.add(holder[3]);
					parentSNPs.add(parents);
				}
			}
			ArrayList<ArrayList<String>> progenySNPs = new ArrayList<ArrayList<String>>();
			
			while ((str = read.readLine()) != null)
			{
				String[] snps = str.split("\t");
				ArrayList<String> holder = new ArrayList<String>();
				double afterN = 0;
				String finder = ":";
				String line = "";
				int end = 0;
				end = snps[0].indexOf(finder) - 4;
			 	line = snps[0].substring(0,end);
			 	ArrayList<String> finalHolder = new ArrayList<String>();
			 	finalHolder.add("Positions");
			 	
			 	//This checks the line and subsets it depending on the population it belongs to
			 	//B73 or Mo17 is always "B", while PHW30 is always "A"
			 	if (line.equals(pops[pop] + "_PHW30") || line.equals("PHW30_" + pops[pop]) || pops[pop].equals("Both"))
			 	{
			 		holder.add(snps[0]);
			 		int tracker = 1;
					for (int i = 0; i<parentSNPs.size(); i++)
					{
						boolean quit = true;
						while (quit)
						{
							if (tracker < snps.length && header[tracker].equals(parentSNPs.get(i).get(0)))
							{	
								
								finalHolder.add(header[tracker]);
								if(snps[tracker].equals(parentSNPs.get(i).get(1)) && !snps[tracker].equals("N"))
								{
									holder.add( "B" );
									quit = false;
								}
								else if (snps[tracker].equals(parentSNPs.get(i).get(2)))
								{
									holder.add ( "A" );
									quit = false;
								}
								//Converts all het calls to AB
								else if (snps[tracker].equals("R") || snps[tracker].equals("Y") || snps[tracker].equals("S") || snps[tracker].equals("W") || snps[tracker].equals("K") || snps[tracker].equals("M"))
								{
									holder.add ("AB");
									quit = false;
								}
								else
								{
									holder.add("N");
									afterN = afterN + 1;
									quit = false;
								}
								tracker = tracker + 1;
							}
							else
							{
								tracker = tracker + 1;
							}
						}
						
					}
					progenySNPs.add(holder);
			 	}
			 	finalParentSNPs = finalHolder;
			}
			
			int progRow = progenySNPs.size();
			int progCol = progenySNPs.get(0).size();
			
			System.out.println( progRow + "\t" + progCol);
			ArrayList<String> snpType = new ArrayList<>(Arrays.asList("A","B","AB","N"));
			
			//This is going to output all the lines and a summary for each one.
			lineOut.write("Lines" + "\t");
			for (int x = 0; x < snpType.size(); x++)
			{
				lineOut.write(snpType.get(x) + "\t");
			}
			lineOut.write("\n");
			
			//This will go through each line of the progeny file. It will count the SNP type for each line and each SNP			
		 	for(int i = 0; i < progRow ; i++)
			{
		 		//This is the array that contains the SNP type summary for each line
		 		double countType [] = new double [snpType.size()];
			 	
								 	
			 	//This is going to go through each type of SNP, and then go through the entire line and see if it equal to that type, if it is, the countType gets +1
				for (int j = 1; j < progCol; j++)
			 	{
					int x = 0;
					Boolean quit = true;
					while (quit && x <snpType.size())
					{	
						if (progenySNPs.get(i).get(j).equals(snpType.get(x)))
			 			{
			 				countType [x] = countType [x] + 1.0;
			 				quit = false;
			 			}
						else
						{
							x = x+1;
						}
			 		}			 		
			 	}
				lineOut.write(progenySNPs.get(i).get(0) + "\t");
				for (int x = 0; x<snpType.size(); x++)
				{
					lineOut.write(countType[x] + "\t");
				}
				lineOut.write("\n");
				
			}
		 	lineOut.close();
		 	
		 	//This outputs all the SNPs and a summary of the calls
			snpOut.write("Positions" + "\t");
			
			for (int x = 0; x < snpType.size(); x++)
			{
				snpOut.write(snpType.get(x) + "\t");
			}
			snpOut.write("\n");
			
			int numLines = progenySNPs.size();
			int numSNPs = progenySNPs.get(0).size();
			System.out.println("This is the summary:" + numLines + "\t" + numSNPs);
			
			for (int i = 0; i < numSNPs; i++)
			{	
				transSNPs.write(i + "\t");
				for (int j = 0; j < numLines ; j++)
				{
					transSNPs.write(progenySNPs.get(j).get(i) + "\t" );
				}
				transSNPs.write("\n");
			}
			transSNPs.close();
			ArrayList<Integer> snpTracker = new ArrayList<Integer>();
			
		 	for (int j = 1; j < progCol; j++)
		 	{
		 		double countType [] = new double [snpType.size()];
		 		for (int i = 0; i < progRow; i++)
		 		{
		 			int x = 0;
					Boolean quit = true;
					while (quit && x <snpType.size())
					{	
						if (progenySNPs.get(i).get(j).equals(snpType.get(x)))
			 			{
			 				countType [x] = countType [x] + 1.0;
			 				quit = false;
			 			}
						else
						{
							x = x+1;
						}
			 		}
		 		}
		 		double percentMiss = countType[3]/(progRow);
		 		
		 		//testStat1 tests for the 3 class method (A vs AB vs B)
		 		//testStat2 tests for the 2 class method (A & B vs AB)
		 		//testStat3 tests for alt 2 class method (A vs B)
		 		double missData = 0.05;
		 		double testStat1 = 3.841 ; //this is 3 classes, 1 df, .05 = 3.841, .01 = 6.6348 , .001 = 10.8275, .0001 = 15.1366, .00001 = 19.5114, .000001 = 23.9284
		 									 // .000000001 = 35.9999
		 		double testStat2 = 6.635; //this is for 2 classes, so 1 df, .5% = 7.879,  1% = 6.635, 5% = 3.841
		 		double testStat3 = 3.841; //this is for 2 classes, so 1 df, .5% = 7.879,  1% = 6.635, 5% = 3.841
		 		
		 		 		
		 		//True chi square filtering, 3 classes
		 		
		 		double obsHet = countType[2];
		 		double obsHom1 = countType[0];
		 		double obsHom2 = countType[1];
		 		double total = obsHet + obsHom1 + obsHom2;
		 		double total2 = obsHom1 + obsHom2;
		 		double expectedHet = total/2;
		 		double expectedHom = total/4;
		 		double expectedHom2 = total2/2;
		 		double obsTestStat1 = ((Math.pow(obsHet - expectedHet, 2))/expectedHet) + ((Math.pow(obsHom1 - expectedHom, 2))/expectedHom) + 
		 				((Math.pow(obsHom2 - expectedHom, 2))/expectedHom);
		 		double obsTestStat3 = ((Math.pow(obsHom1 - expectedHom2,2))/expectedHom2) + ((Math.pow(obsHom2 - expectedHom2,2))/expectedHom2);
		 		
		 		//alternative filtering. Only 2 chi square classes for filtering.
		 		
		 		obsHet = countType[2];
		 		double obsHom = countType[0] + countType[1];
		 		total = obsHet + obsHom;
		 		expectedHet = total/2;
		 		expectedHom = total/2;
		 		double obsTestStat = ((Math.pow(obsHet - expectedHet, 2))/expectedHet) + ((Math.pow(obsHom - expectedHom, 2))/expectedHom);
		 		
		 		
		 		//This filters out missing data, Less than 5%, and checks if the expected ratio at a 5% significance is greater than the observed ratio
		 		//CHANGE the TESTSTAT Variables depending on which class you are using!!
		 		
		 		if (percentMiss <= missData &&  obsTestStat1 <= testStat1)
		 		{
		 			snpTracker.add(j);
		 		}
		 		snpOut.write(finalParentSNPs.get(j) + "\t");
				for (int x = 0; x<snpType.size(); x++)
				{
					snpOut.write(countType[x] + "\t");
				}
				snpOut.write("\n");
				
		 		
		 		/*
		 		snpTracker.add(j);
		 		snpOut.write(finalParentSNPs.get(j) + "\t");
				for (int x = 0; x<snpType.size(); x++)
				{
					snpOut.write(countType[x] + "\t");
				}
				snpOut.write("\n");		 	
				*/
			}
		 	snpOut.close();
		 	


		 	
		 	//This outputs the filtered SNP positions
			finalSNPs.write("SNPs" + "\n");
			
		 	//This outputs the filtered SNP matrix SNPS by lines
			out.write("Lines" + "\t");

		 	for(int j =0; j < snpTracker.size(); j++)
			{
				out.write(finalParentSNPs.get(snpTracker.get(j)) + "\t");
				finalSNPs.write(finalParentSNPs.get(snpTracker.get(j)) + "\n");
			}
			
		 	out.write("\n");
			finalSNPs.close();
			
			for (int i = 0; i < progRow; i++)
			{
				out.write(progenySNPs.get(i).get(0) + "\t");
				for (int j = 0; j < snpTracker.size() ; j++)
				{
					out.write(progenySNPs.get(i).get(snpTracker.get(j)) + "\t" );
				}
				out.write("\n");
			}
			
			System.out.println(progenySNPs.get(0).size());
			System.out.println(progenySNPs.size());
			System.out.println(snpTracker.size());
			out.close();
			System.out.println("Remaining SNPs: " + snpTracker.size());
			System.out.println( "Done with " + pops[pop]);
		}
	}
}



