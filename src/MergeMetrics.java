import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Scanner;
import java.util.TreeSet;

public class MergeMetrics {
public static void main(String[] args) throws Exception
{
	String mergedVcf = "/home/mkirsche/eichlersim/surv.vcf";
	String filelist = "/home/mkirsche/eichler/filelist.txt";
	
	System.err.println("Reading merged VCF");
	
	int intraSampleMerges = 0;
	Scanner input = new Scanner(new FileInputStream(new File(mergedVcf)));
	TreeSet<MergeComparison.MergedVariant> merges = new TreeSet<MergeComparison.MergedVariant>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0)
		{
			continue;
		}
		if(line.startsWith("#"))
		{
			continue;
		}
		VcfEntry entry = new VcfEntry(line);
		merges.add(new MergeComparison.MergedVariant(line));
		if(entry.getInfo("SVMETHOD").startsWith("SURVIVOR"))
		{
			for(int i = 9; i<entry.tabTokens.length; i++)
			{
				String[] cur = entry.tabTokens[i].split(":");
				intraSampleMerges += cur[cur.length-1].split(",").length - 1;
			}
		}
	}
	input.close();
	System.err.println("Found " + merges.size() + " total merged variants");
	
	System.err.println("Reading original files");
	ArrayList<HashMap<String, VcfEntry>> idToEntry = new ArrayList<HashMap<String, VcfEntry>>();
	Scanner filelistInput = new Scanner(new FileInputStream(new File(filelist)));
	while(filelistInput.hasNext())
	{
		String filename = filelistInput.nextLine();
		if(filename.length() == 0)
		{
			continue;
		}
		System.err.println("Reading " + filename);
		HashMap<String, VcfEntry> map = new HashMap<String, VcfEntry>();
		input = new Scanner(new FileInputStream(new File(filename)));
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				continue;
			}
			VcfEntry entry = new VcfEntry(line);
			entry.originalLine = null;
			map.put(entry.getId(), entry);
		}
		System.err.println("Found " + map.size() + " variants");
		input.close();
		idToEntry.add(map);
	}
	filelistInput.close();
	
	int[] suppFreq = new int[16];
	int[][] pairwiseMergingCounts = new int[15][15];
	int interSampleMerges = 0;
	
	System.err.println("Processing merged variants");
	int extremeMerges = 0;
	
	// Now we are iterating over all merged variants and can look at the entries that make them up
	for(MergeComparison.MergedVariant merge : merges)
	{
		int[] samples = merge.samples;
		interSampleMerges += samples.length - 1;
		for(int i = 0; i<samples.length; i++)
			for(int j = 0; j<samples.length; j++)
			{
				if(i != j)
				{
					pairwiseMergingCounts[samples[i]][samples[j]]++;
				}
			}
		String[] ids = merge.ids;
		suppFreq[samples.length]++;
		VcfEntry[] entries = new VcfEntry[samples.length];
		for(int i = 0; i<entries.length; i++)
		{
			entries[i] = idToEntry.get(samples[i]).get(ids[i]);
		}
		long minStart = entries[0].getPos(), maxStart = entries[0].getPos();
		for(VcfEntry entry : entries)
		{
			minStart = Math.min(minStart,  entry.getPos());
			maxStart = Math.max(maxStart, entry.getPos());
		}
		if(maxStart - minStart > 3000)
		{
			extremeMerges++;
		}
	}
	
	System.out.println();
	System.out.println("Intersample merges: " + interSampleMerges);
	System.out.println("Intrasample merges: " + intraSampleMerges);
	System.out.println("Merges spanning >3k bases: " + extremeMerges);
	System.out.println("Component size histogram: " + Arrays.toString(suppFreq));
	System.out.println();
	System.out.println("Merges per sample pair (table[i][j] is percentage of SVs in sample j which get merged with something in sample i)");
	for(int i = 0; i<pairwiseMergingCounts.length; i++)
	{
		for(int j = 0; j<pairwiseMergingCounts[i].length; j++)
		{
			System.out.printf("%05.2f%% ", pairwiseMergingCounts[i][j] * 100.0 / idToEntry.get(j).size());
		}
		System.out.println();
	}
}
}
