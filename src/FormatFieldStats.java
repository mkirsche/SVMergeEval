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
		int oneType = 0, multipleTypes = 0;
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
			String[] tokens = line.split("\t");
			for(int i = 9; i<tokens.length; i++)
			{
				if(entry.getInfo("SVMETHOD").equals("JASMINE"))
				{
					String type = tokens[i].split(":")[2];
					if(type.equals("NA") || type.equals("."))
					{
						continue;
					}
					types.add(type);
				}
				else if(entry.getInfo("SVMETHOD").startsWith("SURVIVOR"))
				{
					String type = tokens[i].split(":")[6];
					if(type.equalsIgnoreCase("NAN") || type.equals("."))
					{
						continue;
					}
					String[] subTypes = type.split(",");
					for(String subType : subTypes)
					{
						types.add(subType);
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
		}
		
		System.out.println("Merged variants of all one type: " + oneType);
		System.out.println("Merged variants with multiple types: " + multipleTypes);
		for(String s : multiTypeSets.keySet())
		{
			System.out.println(s+" "+multiTypeSets.get(s));
		}
		
		input.close();
	}
}
