/*
 * Gets information from Jasmine output based on its FORMAT fields
 * For now just checks type uniqueness
 */

import java.io.File;
import java.io.FileInputStream;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.TreeSet;

public class FormatFieldStats
{
	public static void main(String[] args) throws Exception
	{
		//String mergedGtVcf = "/home/mkirsche/eichler/survmerged.vcf";
		String mergedGtVcf = "/home/mkirsche/eichler/mergedgt.vcf";
		
		if(args.length == 1)
		{
			mergedGtVcf = args[0];
		}
		else
		{
			System.out.println("Usage: java FormatFieldStats <mergedVcfWithGenotypes>");
		}
		int oneType = 0, multipleTypes = 0;
		int oneGt = 0, multipleGt = 0;
		int uniqueSpecific = 0, uniqueSensitive = 0;
		Scanner input = new Scanner(new FileInputStream(new File(mergedGtVcf)));
		TreeMap<String, Integer> multiTypeSets = new TreeMap<String, Integer>();
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				continue;
			}
			VcfEntry entry = new VcfEntry(line);
			TreeSet<String> types = new TreeSet<String>();
			TreeSet<String> gts = new TreeSet<String>();
			String[] tokens = line.split("\t");
			
			String suppVec = entry.getInfo("SUPP_VEC");
			int numSamples = countOnes(suppVec);
			int firstOneIndex = suppVec.indexOf('1');
			
			if(numSamples == 1)
			{
				char specificCharacter = tokens[firstOneIndex + 9].split(":")[1].charAt(0);
				if(specificCharacter == '0')
				{
					uniqueSensitive++;
				}
				else
				{
					uniqueSpecific++;
				}
			}
			for(int i = 9; i<tokens.length; i++)
			{
				if(entry.getInfo("SVMETHOD").equals("JASMINE"))
				{
					String[] formatTokens = tokens[i].split(":");
					String type = formatTokens[2];
					if(!type.equals("NA") && !type.equals("."))
					{
						types.add(type);
					}
					
					String gt = formatTokens[0];
					if(!gt.equals("./."))
					{
						gts.add(gt);
					}
				}
				else if(entry.getInfo("SVMETHOD").startsWith("SURVIVOR"))
				{
					String[] formatTokens = tokens[i].split(":");
					String type = formatTokens[6];
					if(!type.equalsIgnoreCase("NAN") && !type.equals("."))
					{
						String[] subTypes = type.split(",");
						for(String subType : subTypes)
						{
							types.add(subType);
						}
					}
					String gt = formatTokens[0];
					if(!gt.equals("./."))
					{
						gts.add(gt);
					}
				}
			}
			if(types.size() == 1)
			{
				oneType++;
			}
			else if(types.size() > 1)
			{
				multipleTypes++;
				int newCount = multiTypeSets.containsKey(types.toString()) ? (1 + multiTypeSets.get(types.toString())) : 1;
				multiTypeSets.put(types.toString(), newCount);
			}
			if(gts.size() == 1)
			{
				oneGt++;
			}
			else if(gts.size() > 1)
			{
				multipleGt++;
			}
		}
		
		System.out.println("Sample-unique specific calls: " + uniqueSpecific);
		System.out.println("Sample-unique sensitive calls: " + uniqueSensitive);
		System.out.println("Merged variants of all one genotype: " + oneGt);
		System.out.println("Merged variants with multiple genotypes: " + multipleGt);
		System.out.println("Merged variants of all one type: " + oneType);
		System.out.println("Merged variants with multiple types: " + multipleTypes);
		for(String s : multiTypeSets.keySet())
		{
			System.out.println(s+" "+multiTypeSets.get(s));
		}
		
		input.close();
	}
	
	/*
	 * The number of times the character '1' occurs in a string
	 */
	static int countOnes(String s)
	{
		int res = 0;
		for(int i = 0; i<s.length(); i++)
		{
			if(s.charAt(i) == '1')
			{
				res++;
			}
		}
		return res;
	}
}
