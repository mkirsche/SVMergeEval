/*
 * Script for testing the order-independence of variant mergers
 * It runs a given merger with a list of samples, plus that list reversed in different ways 
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;

public class ReverseComparison {
	static String merger = "JASMINE";
	//static String merger = "SURVIVOR";
	static String revMode = "FILELIST";
	//static String revMode = "VCF";
	//static String revMode = "SHUFFLE";
	static boolean overwrite = false;
public static void main(String[] args) throws Exception
{
	String fileList = "", revFileList = "", forwardVcf = "", reverseVcf = "";
	if(args.length == 4)
	{
		fileList = args[0];
		revFileList = args[1];
		forwardVcf = args[2];
		reverseVcf = args[3];
	}
	else
	{
		System.out.println("Usage: java ReverseComparison <filelist> <revfilelist> <forwardvcf> <reversevcf>");
		if(merger.equals("JASMINE"))
		{
			if(revMode.equals("FILELIST"))
			{
				fileList = "/home/mkirsche/eichler/filelist.txt";
				revFileList = "/home/mkirsche/eichler/jasmine_revfilelist.txt";
				forwardVcf = "/home/mkirsche/eichler/jasmine_merged.vcf";
				reverseVcf = "/home/mkirsche/eichler/jasmine_revmerged.vcf";
				compareReverseFilelist(fileList, revFileList, forwardVcf, reverseVcf);
			}
			else if(revMode.equals("VCF"))
			{
				fileList = "/home/mkirsche/eichler/filelist.txt";
				revFileList = "/home/mkirsche/eichler/jasmine_eachrevfilelist.txt";
				forwardVcf = "/home/mkirsche/eichler/jasmine_merged.vcf";
				reverseVcf = "/home/mkirsche/eichler/jasmine_eachrevmerged.vcf";
				compareReverseFilelist(fileList, revFileList, forwardVcf, reverseVcf);
			}
			else if(revMode.equals("SHUFFLE"))
			{
				fileList = "/home/mkirsche/eichler/filelist.txt";
				revFileList = "/home/mkirsche/eichler/jasmine_shufflefilelist.txt";
				forwardVcf = "/home/mkirsche/eichler/jasmine_merged.vcf";
				reverseVcf = "/home/mkirsche/eichler/jasmine_shufflemerged.vcf";
				compareReverseFilelist(fileList, revFileList, forwardVcf, reverseVcf);
			}
		}
		else if(merger.equals("SURVIVOR"))
		{
			if(revMode.equals("FILELIST"))
			{
				fileList = "/home/mkirsche/eichler/filelist.txt";
				revFileList = "/home/mkirsche/eichler/survivor_revfilelist.txt";
				forwardVcf = "/home/mkirsche/eichler/survivor_merged.vcf";
				reverseVcf = "/home/mkirsche/eichler/survivor_revmerged.vcf";
				compareReverseFilelist(fileList, revFileList, forwardVcf, reverseVcf);
			}
			else if(revMode.equals("VCF"))
			{
				fileList = "/home/mkirsche/eichler/filelist.txt";
				revFileList = "/home/mkirsche/eichler/survivor_eachrevfilelist.txt";
				forwardVcf = "/home/mkirsche/eichler/survivor_merged.vcf";
				reverseVcf = "/home/mkirsche/eichler/survivor_eachrevmerged.vcf";
				compareReverseFilelist(fileList, revFileList, forwardVcf, reverseVcf);
			}
			else if(revMode.equals("SHUFFLE"))
			{
				fileList = "/home/mkirsche/eichler/filelist.txt";
				revFileList = "/home/mkirsche/eichler/survivor_shufflefilelist.txt";
				forwardVcf = "/home/mkirsche/eichler/survivor_merged.vcf";
				reverseVcf = "/home/mkirsche/eichler/survivor_shufflemerged.vcf";
				compareReverseFilelist(fileList, revFileList, forwardVcf, reverseVcf);
			}
		}
	}
}

/*
 * Compare a merged variant set on a given filelist to the callset on the same filelist in reverse order
 */
static void compareReverseFilelist(String fileList, String revFileList, String forwardVcf, String reverseVcf) throws Exception
{
	int[] order = new int[] {};
	if(revMode.equals("FILELIST"))
	{
		order = revFileList(fileList, revFileList);
	}
	else if(revMode.equals("VCF"))
	{
		revAllFiles(fileList, revFileList);
	}
	else if(revMode.equals("SHUFFLE"))
	{
		order = shuffleAll(fileList, revFileList);
	}
	if(merger.equals("JASMINE"))
	{
		if(overwrite || !(new File(forwardVcf)).exists())
		{
			RunMergingSoftware.runJasmine(fileList, forwardVcf);
		}
		else
		{
			System.out.println("Jasmine output already exists: " + forwardVcf);
		}
		if(overwrite || !(new File(reverseVcf)).exists())
		{
			RunMergingSoftware.runJasmine(revFileList, reverseVcf);
		}
		else
		{
			System.out.println("Jasmine output already exists: " + reverseVcf);
		}
	}
	else if(merger.equals("SURVIVOR"))
	{
		if(overwrite || !(new File(forwardVcf)).exists())
		{
			RunMergingSoftware.runSurvivor(fileList, forwardVcf);
		}
		if(overwrite || !(new File(reverseVcf)).exists())
		{
			RunMergingSoftware.runSurvivor(revFileList, reverseVcf);
		}
	}
	MergeComparison.compareMergedSets(forwardVcf, reverseVcf, order);
}

/*
 * Reverses the order of files in a list and outputs the resulting list to a new file
 */
static int[] revFileList(String inputFile, String outputFile) throws Exception
{
	Scanner input = new Scanner(new FileInputStream(new File(inputFile)));
	PrintWriter out = new PrintWriter(new File(outputFile));
	ArrayList<String> lines = new ArrayList<String>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0)
		{
			continue;
		}
		lines.add(line);
	}
	Collections.reverse(lines);
	for(String line : lines)
	{
		out.println(line);
	}
	input.close();
	out.close();
	
	int n = lines.size();
	int[] order = new int[n];
	for(int i = 0; i<n; i++)
	{
		order[i] = n - 1 - i;
	}
	return order;
}

/*
 * Reverse the order of variants within each VCF file
 */
static void revAllFiles(String inputFile, String outputFile) throws Exception
{
	Scanner fileListInput = new Scanner(new FileInputStream(new File(inputFile)));
	PrintWriter fileListOut = new PrintWriter(new File(outputFile));
	while(fileListInput.hasNext())
	{
		String fileName = fileListInput.nextLine();
		if(fileName.length() == 0)
		{
			continue;
		}
		String newFileName = fileName + ".rev.vcf";
		Scanner input = new Scanner(new FileInputStream(new File(fileName)));
		PrintWriter out = new PrintWriter(new File(newFileName));
		ArrayList<String> entries = new ArrayList<String>();
		while(input.hasNext())
		{
			String vcfLine = input.nextLine();
			if(vcfLine.length() == 0 || vcfLine.startsWith("#"))
			{
				out.println(vcfLine);
			}
			else
			{
				entries.add(vcfLine);
			}
		}
		Collections.reverse(entries);
		for(String vcfLine : entries)
		{
			out.println(vcfLine);
		}
		
		input.close();
		out.close();
		
		fileListOut.println(newFileName);
	}
	fileListInput.close();
	fileListOut.close();
}

/*
 * Shuffles both the order of VCF files and the order of variants within each VCF file
 */
static int[] shuffleAll(String inputFile, String outputFile) throws Exception
{
	Scanner fileListInput = new Scanner(new FileInputStream(new File(inputFile)));
	PrintWriter fileListOut = new PrintWriter(new File(outputFile));
	ArrayList<String> fileNames = new ArrayList<String>();
	while(fileListInput.hasNext())
	{
		String fileName = fileListInput.nextLine();
		if(fileName.length() == 0)
		{
			continue;
		}
		String newFileName = fileName + ".shuffle.vcf";
		Scanner input = new Scanner(new FileInputStream(new File(fileName)));
		PrintWriter out = new PrintWriter(new File(newFileName));
		ArrayList<String> entries = new ArrayList<String>();
		while(input.hasNext())
		{
			String vcfLine = input.nextLine();
			if(vcfLine.length() == 0 || vcfLine.startsWith("#"))
			{
				out.println(vcfLine);
			}
			else
			{
				entries.add(vcfLine);
			}
		}
		Collections.shuffle(entries, new Random(41));
		for(String vcfLine : entries)
		{
			out.println(vcfLine);
		}
		
		input.close();
		out.close();
		
		fileNames.add(newFileName);
	}
	ArrayList<Integer> orderList = new ArrayList<Integer>();
	for(int i = 0; i<fileNames.size(); i++)
	{
		orderList.add(i);
	}
	ArrayList<String> newFileNames = new ArrayList<String>();
	Collections.shuffle(orderList, new Random(67));
	for(int i = 0; i<orderList.size(); i++)
	{
		newFileNames.add(fileNames.get(orderList.get(i)));
	}
	fileNames = newFileNames;
	for(String filename : fileNames)
	{
		fileListOut.println(filename);
	}
	fileListInput.close();
	fileListOut.close();
	
	int n = orderList.size();
	int[] order = new int[n];
	for(int i = 0; i<n; i++)
	{
		order[i] = orderList.get(i);
		//order[orderList.get(i)] = i;
	}
	return order;
}
}
