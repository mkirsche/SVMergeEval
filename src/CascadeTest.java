import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Scanner;

public class CascadeTest
{
	
	static boolean RERUN_MERGING = true;
	
	public static void main(String[] args) throws Exception
	{
		String fileListPath = "/home/mkirsche/eichler/filelist.txt";
		if(args.length == 1)
		{
			fileListPath = args[0];
		}
		else
		{
			System.out.println("Usage: java CascadeTest <fileListPath>");
		}
		runCascadeTest(fileListPath);
	}
	
	/*
	 * Tests the effects of running a sequence of pairwise merges rather than merging everything at once
	 */
	static void runCascadeTest(String fileListPath) throws Exception
	{
		Scanner fileListInput = new Scanner(new FileInputStream(new File(fileListPath)));
		ArrayList<String> allFileNames = new ArrayList<String>();
		while(fileListInput.hasNext())
		{
			String line = fileListInput.nextLine();
			if(line.length() > 0)
			{
				allFileNames.add(line);
			}
		}
		fileListInput.close();
		
		Path currentRelativePath = Paths.get("");
		String outDir = currentRelativePath.toAbsolutePath().toString() + "/" + "merged";
		
		String cascadeMergeOut = "";
		
		// Run the cascading merging
		File f = new File(outDir);
		f.mkdir();
					
		for(int i = 1; i < allFileNames.size(); i++)
		{
			String oldFile = (i == 1) ? allFileNames.get(0) : outFileName(outDir, i-1);
			String newFile = allFileNames.get(i);
			System.out.println("Merging " + oldFile + " and " + newFile);
			String fileList = outDir + "/" + "filelist" + i + ".txt";
			PrintWriter fileListOut = new PrintWriter(new File(fileList));
			fileListOut.println(oldFile);
			fileListOut.println(newFile);
			fileListOut.close();
			String outFileName = outFileName(outDir, i);
			cascadeMergeOut = outFileName;
			if(new File(outFileName).exists() && !RERUN_MERGING)
			{
				System.out.println("Already merged");
				continue;
			}
			RunMergingSoftware.runJasmine(fileList, outFileName);
		}
		
		String fullOut = outDir + "/allMerged.vcf";
		if(!new File(fullOut).exists() || RERUN_MERGING)
		{
			System.out.println("Merging everything at once");
			RunMergingSoftware.runJasmine(fileListPath, fullOut);
		}
		else
		{
			System.out.println("Full merge file already exists, so reusing that");
		}
		
		MergeComparison.USE_EXT_IDLISTS[1] = true;
		MergeComparison.compareMergedSets(fullOut, cascadeMergeOut, new int[] {});
	}
	
	static String outFileName(String outDir, int mergeNum)
	{
		return outDir + "/" + "merged" + mergeNum + ".vcf";
	}

}
