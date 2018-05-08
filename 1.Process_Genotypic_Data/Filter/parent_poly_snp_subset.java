package Filter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class parent_poly_snp_subset {

	public static void main(String[] args) throws IOException {

		
		String pops [] = new String []{"B73","Mo17","Both"};
		
		for (int pop = 0; pop < pops.length  ; pop++) 
		{	 		
			// This is the master file of the parent snps, it is column (parents) rows (SNPs)
			String parentFile = "C:\\Users\\mdzievit.IASTATE\\Desktop\\F3_SNPs\\parents_final_SNPs.txt";
			 
			BufferedReader read = new BufferedReader(new FileReader(parentFile));
			String str = read.readLine();;
			String holder [] = new String [4];
			ArrayList<ArrayList<String>> snps = new ArrayList<ArrayList<String>>();
			int counter = 0;
			while ((str = read.readLine()) != null)
		 	{
				holder = str.split("\t");
				if(holder[3].equals("A") || holder[3].equals("C") || holder[3].equals("G") || holder[3].equals("T"))
				{
					if ( pop == 0 || pop == 1)
					{
						if(holder[pop+1].equals("A") || holder[pop+1].equals("C") || holder[pop + 1].equals("G") || holder[pop + 1].equals("T"))
						{
							if (!holder[pop+1].equals(holder[3]))
							{
								snps.add( new ArrayList<String>());
								for (int j = 0; j <4; j++)
								{
									snps.get(counter).add(holder[j]);							 
								}
								counter = counter + 1;
							}
						}
					}
					
					else if (pop==2)
					{
						if(holder[pop].equals("A") || holder[pop].equals("C") || holder[pop].equals("G") || holder[pop].equals("T"))
						{
							if(holder[pop - 1].equals("A") || holder[pop - 1].equals("C") || holder[pop - 1].equals("G") || holder[pop - 1].equals("T"))
							{
								if( !holder[1].equals(holder[3]) && !holder[2].equals(holder[3]))
								{
									snps.add( new ArrayList<String>());
									for (int j = 0; j <4; j++)
									{
										snps.get(counter).add(holder[j]);							 
									}
									counter = counter + 1;
								}
							}
						}	
					}	
				}
		 	}

			String outputFilename = "C:\\Users\\mdzievit.IASTATE\\Desktop\\F3_SNPs\\" + pops[pop] + "\\parents_final_SNPs_filtered_" + pops[pop] + ".txt";
			BufferedWriter out = new BufferedWriter (new FileWriter(outputFilename,true));
			out.write("Position" + "\t" + "B73" + "\t" + "Mo17" + "\t" + "PWH30" + "\t" + "\n");
			
			for (int i = 0; i <snps.size() ; i++)
			{
				for (int j = 0; j <snps.get(0).size(); j++)
				{
					out.write( snps.get(i).get(j) + "\t");
				}
				out.write("\n");
			}
			out.close();
			System.out.println(pops[pop] + "\t" + snps.size());
			System.out.println("All Done I Hope!!");		

		}
	}

}
