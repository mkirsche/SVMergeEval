/*
 * Script for testing the order-independence of variant mergers
 * It runs a given merger with a list of samples, plus that list reversed in different ways 
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

public class ReverseComparison {
	//static String merger = "JASMINE";
	static String merger = "SURVIVOR";
	//static String revMode = "FILELIST";
	static String revMode = "VCF";
	static boolean overwrite = false;
public static void main(String[] args) throws Exception
{
	if(merger.equals("JASMINE"))
	{
		if(revMode.equals("FILELIST"))
		{
			String fileList = "/home/mkirsche/eichler/filelist.txt";
			String revFileList = "/home/mkirsche/eichler/jasmine_revfilelist.txt";
			String forwardVcf = "/home/mkirsche/eichler/jasmine_merged.vcf";
			String reverseVcf = "/home/mkirsche/eichler/jasmine_revmerged.vcf";
			compareReverseFilelist(fileList, revFileList, forwardVcf, reverseVcf);
		}
		else if(revMode.equals("VCF"))
		{
			String fileList = "/home/mkirsche/eichler/filelist.txt";
			String revFileList = "/home/mkirsche/eichler/jasmine_eachrevfilelist.txt";
			String forwardVcf = "/home/mkirsche/eichler/jasmine_merged.vcf";
			String reverseVcf = "/home/mkirsche/eichler/jasmine_eachrevmerged.vcf";
			compareReverseFilelist(fileList, revFileList, forwardVcf, reverseVcf);
		}
	}
	else if(merger.equals("SURVIVOR"))
	{
		if(revMode.equals("FILELIST"))
		{
			String fileList = "/home/mkirsche/eichler/filelist.txt";
			String revFileList = "/home/mkirsche/eichler/survivor_revfilelist.txt";
			String forwardVcf = "/home/mkirsche/eichler/survivor_merged.vcf";
			String reverseVcf = "/home/mkirsche/eichler/survivor_revmerged.vcf";
			compareReverseFilelist(fileList, revFileList, forwardVcf, reverseVcf);
		}
		else if(revMode.equals("VCF"))
		{
			String fileList = "/home/mkirsche/eichler/filelist.txt";
			String revFileList = "/home/mkirsche/eichler/survivor_eachrevfilelist.txt";
			String forwardVcf = "/home/mkirsche/eichler/survivor_merged.vcf";
			String reverseVcf = "/home/mkirsche/eichler/survivor_eachrevmerged.vcf";
			compareReverseFilelist(fileList, revFileList, forwardVcf, reverseVcf);
		}
	}
}

/*
 * Compare a merged variant set on a given filelist to the callset on the same filelist in reverse order
 */
static void compareReverseFilelist(String fileList, String revFileList, String forwardVcf, String reverseVcf) throws Exception
{
	if(revMode.equals("FILELIST"))
	{
		if(overwrite || !(new File(revFileList)).exists())
		{
			revFileList(fileList, revFileList);
		}
	}
	else if(revMode.equals("VCF"))
	{
		if(overwrite || !(new File(revFileList)).exists())
		{
			revFileList(fileList, revFileList);
		}
	}
	if(merger.equals("JASMINE"))
	{
		if(overwrite || !(new File(forwardVcf)).exists())
		{
			runJasmine(fileList, forwardVcf);
		}
		else
		{
			System.out.println("Jasmine output already exists: " + forwardVcf);
		}
		if(overwrite || !(new File(reverseVcf)).exists())
		{
			runJasmine(revFileList, reverseVcf);
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
			runSurvivor(fileList, forwardVcf);
		}
		if(overwrite || !(new File(reverseVcf)).exists())
		{
			runSurvivor(revFileList, reverseVcf);
		}
	}
	if(revMode.equals("FILELIST"))
	{
		MergeComparison.compareMergedSets(forwardVcf, reverseVcf, true);
	}
	else if(revMode.equals("VCF"))
	{
		MergeComparison.compareMergedSets(forwardVcf, reverseVcf, false);
	}
}

/*
 * Reverses the order of files in a list and outputs the resulting list to a new file
 */
static void revFileList(String inputFile, String outputFile) throws Exception
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
 * Runs Jasmine with min_support=1 on given input and output files
 */
static void runJasmine(String fileList, String outFile) throws Exception
{
	System.out.println("Running Jasmine on " + fileList);
	String command = "/usr/lib/jvm/jdk-12.0.2/bin/java -cp /home/mkirsche/eclipse-workspace/Thriver/src:/home/mkirsche/eclipse-workspace/Thriver/Iris/src Main "
	 + "file_list=" + fileList + " out_file=" + outFile + " min_support=1";
	Process child = Runtime.getRuntime().exec(command);
    int exit = child.waitFor();
    if(exit != 0)
    {
    	System.out.println("Command failed: " + command);
    	System.out.println("Exit code: " + exit);
    	System.exit(1);
    }
    else
    {
    	System.out.println("Jasmine ran successfully");
    }
}

/*
 * Runs SURVIVOR with max_dist=1000 on given input and output files
 */
static void runSurvivor(String fileList, String outFile) throws  Exception
{
	System.out.println("Running SURVIVOR on " + fileList);
	String command = "/home/mkirsche/git/SURVIVOR/Debug/SURVIVOR merge " + fileList + " 1000 1 1 1 1 1 " + outFile;
	Process child = Runtime.getRuntime().exec(command);
    int exit = child.waitFor();
    if(exit != 0)
    {
    	System.out.println("Command failed: " + command);
    	System.out.println("Exit code: " + exit);
    	System.exit(1);
    }
    else
    {
    	System.out.println("SURVIVOR ran successfully");
    }
}
}
