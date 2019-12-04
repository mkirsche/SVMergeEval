/*
 * Compares the sets of merge variants produced by two different runs of SV merging software
 */
import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.TreeSet;

public class MergeComparison {
	
	// Option to run in different comparison modes based on the desired similarity metric
	static String[] MODE_DESCRIPTIONS = new String[] {
			"full merged variants",
			"consecutive variants in merged sets",
			"all pairs of merged variants"
	};
	public static void main(String[] args) throws Exception
	{
		String fn1 = "/home/mkirsche/eichlersim/surv.vcf";
		String fn2 = "/home/mkirsche/eichlersim/merged.vcf";
		if(args.length == 2)
		{
			fn1 = args[0];
			fn2 = args[1];
		}
		compareMergedSets(fn1, fn2, false);
	}
	
	public static void compareMergedSets(String fn1, String fn2, boolean fn2Reversed) throws Exception
	{
		for(int comparisonMode = 0; comparisonMode <= 2; comparisonMode++)
		{
			TreeSet<MergedVariant> first = getAllMerged(fn1, comparisonMode, false), second = getAllMerged(fn2, comparisonMode, fn2Reversed);
			int[] counts = subsetCounts(first, second);
			System.out.println("Comparing " + MODE_DESCRIPTIONS[comparisonMode]);
			System.out.printf("First only: %d (%.2f%% of first callset)\n",counts[1], 100.0 * counts[1] / first.size());
			System.out.printf("Second only: %d (%.2f%% of second callset)\n",counts[2], 100.0 * counts[2] / second.size());
			System.out.printf("Both: %d (%.4f%% of larger callset)\n",counts[3], 100.0 * counts[3] / Math.max(first.size(), second.size()));
			System.out.println();
		}
	}
	
	/*
	 * Given two sets of merged variants, outputs four values:
	 *   1. The number of variants in the union of both sets
	 *   2. The number of variants in the first set but not the second
	 *   3. The number of variant in the second set but not the first
	 *   4. The number of variants in both sets
	 */
	static int[] subsetCounts(TreeSet<MergedVariant> a, TreeSet<MergedVariant> b)
	{
		TreeSet<MergedVariant> merged = union(a, b);
		int[] res = new int[4];
		for(MergedVariant mv : merged)
		{
			res[0]++;
			int cur = 0;
			if(a.contains(mv)) cur |= 1;
			if(b.contains(mv)) cur |= 2;
			res[cur]++;
		}
		return res;
	}
	
	/*
	 * Compute the union of two sets of merged variants
	 */
	static TreeSet<MergedVariant> union(TreeSet<MergedVariant> a, TreeSet<MergedVariant> b)
	{
		TreeSet<MergedVariant> res = new TreeSet<MergedVariant>();
		for(MergedVariant mva : a) res.add(mva);
		for(MergedVariant mvb : b) res.add(mvb);
		return res;
	}
	
	/*
	 * Get a set of all merged variants from a VCF file
	 */
	static TreeSet<MergedVariant> getAllMerged(String fn, int comparisonMode, boolean reversed) throws Exception
	{
		TreeSet<MergedVariant> res = new TreeSet<MergedVariant>();
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				continue;
			}
			MergedVariant mv = new MergedVariant(line, reversed);
			if(comparisonMode == 0)
			{
				res.add(mv);
			}
			else if(comparisonMode == 1)
			{
				int length = mv.ids.length;
				for(int i = 0; i<length-1; i++)
				{
					int[] samples = new int[] {mv.samples[i], mv.samples[i+1]};
					String[] ids = new String[] {mv.ids[i], mv.ids[i+1]};
					res.add(new MergedVariant(ids, samples));
				}
			}
			else if(comparisonMode == 2)
			{
				int length = mv.ids.length;
				for(int i = 0; i<length-1; i++)
				{
					for(int j = i+1; j<length; j++)
					{
						int[] samples = new int[] {mv.samples[i], mv.samples[j]};
						String[] ids = new String[] {mv.ids[i], mv.ids[j]};
						res.add(new MergedVariant(ids, samples));
					}
				}
			}
		}
		input.close();
		return res;
	}
	
	/*
	 * A merged variant represented as the samples it drew from and the IDs
	 * of the variants within those samples
	 */
	static class MergedVariant implements Comparable<MergedVariant>
	{
		String[] ids;
		int[] samples;
		MergedVariant(String[] ids, int[] samples)
		{
			this.ids = ids;
			this.samples = samples;
		}
		MergedVariant(String line, boolean reversed) throws Exception
		{
			VcfEntry entry = new VcfEntry(line);
			
			// Check SUPP field for number of samples involved
			int supp = Integer.parseInt(entry.getInfo("SUPP"));
			
			// Get the indices of the samples from SUPP_VEC
			String suppVec = entry.getInfo("SUPP_VEC");
			int[] sampleIndices = new int[supp];
			int idx = 0;
			for(int i = 0; i<suppVec.length(); i++)
			{
				if(suppVec.charAt(i) == '1')
				{
					sampleIndices[idx++] = reversed ? (suppVec.length() - i - 1) : i;
				}
			}
			
			// Get the variant IDs
			String method = entry.getInfo("SVMETHOD");
			String[] idList = new String[supp];
			if(method.equals("JASMINE"))
			{
				idList = entry.getInfo("IDLIST").split(",");
			}
			else if(method.startsWith("SURVIVOR"))
			{
				ArrayList<String> ids = new ArrayList<String>();
				for(int i = 9; i<entry.tabTokens.length; i++)
				{
					String val = entry.tabTokens[i].split(":")[7];
					if(!val.equalsIgnoreCase("nan"))
					{
						ids.add(val);
					}
				}
				for(int i = 0; i<ids.size(); i++)
				{
					idList[i] = ids.get(i);
				}
			}
			
			// Initialize all values
			ids = new String[supp];
			samples = new int[supp];
			for(int i = 0; i<supp; i++)
			{
				ids[i] = idList[i];
				samples[i] = sampleIndices[i];
			}
			for(int i = 0; i<ids.length; i++)
			{
				for(int j = 0; j<ids.length-1; j++)
				{
					if(samples[j] > samples[j+1])
					{
						int tmp = samples[j];
						samples[j] = samples[j+1];
						samples[j+1] = tmp;
						String tmpId = ids[j];
						ids[j] = ids[j+1];
						ids[j+1] = tmpId;
					}
				}
			}
		}
		
		// Order by increasing value of (sample, id) pairs
		public int compareTo(MergedVariant o) {
			for(int i = 0; i < ids.length || i < o.ids.length; i++)
			{
				if(i == ids.length) return -1;
				if(i == o.ids.length) return 1;
				if(samples[i] != o.samples[i])
				{
					return samples[i] - o.samples[i];
				}
				if(!ids[i].equals(o.ids[i]))
				{
					return ids[i].compareTo(o.ids[i]);
				}
			}
			return 0;
		}
	}
}
