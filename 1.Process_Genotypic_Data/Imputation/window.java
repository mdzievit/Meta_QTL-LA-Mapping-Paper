package Imputation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class window {

	public static void main(String[] args) throws IOException {
		String pops [] = new String []{"B73","Mo17"};
		for (int pop = 0; pop< pops.length; pop++) 
		{	 				
			String snpFile = "C:\\Users\\mdzievit.IASTATE\\Box Sync\\Desktop\\F3_SNPs\\" + pops[pop] + "\\Less\\" + "coded_snps_" + pops[pop] + ".txt";
			BufferedReader read = new BufferedReader(new FileReader(snpFile));
			String str;
			str = read.readLine();
			String header [] = str.split("\t");
			
			//***********These are all of the output files:
			//This is the imputed SNPs output file
			String outputFilename = "C:\\Users\\mdzievit.IASTATE\\Box Sync\\Desktop\\F3_SNPs\\" + pops[pop] + "\\Less\\Imputed\\" + "imputed_" + pops[pop] + ".txt";
			BufferedWriter out = new BufferedWriter (new FileWriter(outputFilename,true));
			
			//This is the imputed SNP file with CHR, overall SNP number, and position of the second imputed SNP.
			String snpOutputFilename = "C:\\Users\\mdzievit.IASTATE\\Box Sync\\Desktop\\F3_SNPs\\" + pops[pop] + "\\Less\\Imputed\\" + "imputed_snps_" + pops[pop] + ".txt";
			BufferedWriter snpOut = new BufferedWriter (new FileWriter(snpOutputFilename,true));
			snpOut.write("SNP" + "\t" + "Coded" + "\t" + "CHR" + "\t" + "BP" + "\t" + "\n");
			
			//Creates the matrix for storing imputed calls
			ArrayList<ArrayList<String>> allWindows = new ArrayList<ArrayList<String>>();
			
			//This is an array that holds all the SNP names (positions)
			ArrayList<String> names = new ArrayList<String>();
			
			//These are arrays that hold the CHR, BP, and SNP number to print later.
			ArrayList<String> snpWindow = new ArrayList<String>();
			ArrayList<Integer> snpTrackerHolder = new ArrayList<Integer>();
			ArrayList<Integer> chrTrackerHolder = new ArrayList<Integer>();
			
			while(( str = read.readLine()) != null)
			{
				int step = 1;
				String snps [] = str.split("\t");
				names.add(snps[0]);
				//Holds the lines imputed SNP call 
				ArrayList<String> lineWindow = new ArrayList<String>();
				
				//Holds the SNP position that is being imputed, the second SNP, since the window is moving by 2
				ArrayList<String> snpHolder = new ArrayList<String>();
				
				//Holds the overall SNP number. This will be helpful for comparisons between raw and different imputing methods
				ArrayList<Integer> snpTracker = new ArrayList<Integer>();
				
				//This trackers the chr and adds it to an array.
				ArrayList<Integer> chrTracker = new ArrayList<Integer>();
				//The following is for doing window at the end
				
				int windowSize = 18;
				int snpCriteria = 15;
				
				for (int i = 1; i<snps.length - (windowSize - 1); i=i + step)
				{
					
					int snpStart = Integer.parseInt(header[i]);
					int snpEnd = Integer.parseInt(header[i + (windowSize - 1)]);
					if ( snpStart < snpEnd)
					{
						int Acounter = 0;
						int ABcounter = 0;
						int Bcounter = 0;
						int Ncounter = 0;
						snpHolder.add(header[i + (windowSize - 1)]);
						snpTracker.add(i + (windowSize - 1));
						
						for (int j = 0; j < windowSize; j++)
						{
							if (snps[j+i].equals("A"))
							{
								Acounter = Acounter + 1;
							}
							else if (snps[j+i].equals("AB"))
							{
								ABcounter = ABcounter + 1;			
							}
							else if (snps[j+i].equals("B"))
							{
								Bcounter = Bcounter + 1;
							}
							else if (snps[j+i].equals("N"))
							{
								Ncounter = Ncounter + 1;
							}
						}
						if (Acounter >= snpCriteria)
						{
							lineWindow.add("A");
						}
						else if (Bcounter >= snpCriteria)
						{
							lineWindow.add("B");
						}
						else if (ABcounter > 2)
						{
							lineWindow.add("AB");
						}
												else if (Ncounter > (snpCriteria - 5))
						{
							lineWindow.add("-1");
						}
						else
						{
							//lineWindow.add("-1");
							lineWindow.add("AB");
						}
					}
				}
				
				int chr  = 1;
				chrTracker.add(chr);
				for (int bb = 1 ; bb < snpHolder.size(); bb++)
				{
					int start = Integer.parseInt(snpHolder.get(bb-1));
					int end  = Integer.parseInt(snpHolder.get(bb));
					if(end>start)
					{
						chrTracker.add(chr);
					}
					else
					{
						chr = chr + 1;
						chrTracker.add(chr);
					}
				}
				allWindows.add(lineWindow);
				snpWindow = snpHolder;
				snpTrackerHolder = snpTracker;
				chrTrackerHolder = chrTracker;
			}
			int lines = allWindows.size();
			int snp = allWindows.get(0).size();
			
			out.write("Lines" + "\t");
			for (int i = 0; i <lines ; i++)
			{
				out.write(names.get(i) + "\t");
			}
			out.write("\n");
			
			for (int i = 0; i < snp; i++)
			{
				out.write( snpWindow.get(i) + "\t" );
				for (int j = 0; j < lines; j ++)
				{
					out.write(allWindows.get(j).get(i) + "\t" );
				}
				out.write("\n");
			}
			out.close();
			for (int i = 0; i < snpWindow.size(); i++)
			{
				String coded = snpTrackerHolder.get(i) + "_" + chrTrackerHolder.get(i) + "_" + snpWindow.get(i);
				snpOut.write(snpTrackerHolder.get(i) + "\t" + coded + "\t" + chrTrackerHolder.get(i) + "\t" + snpWindow.get(i) + "\t" + "\n");
			}
			snpOut.close();
			System.out.println("All done with " + pops[pop]);
		}
	}

}
