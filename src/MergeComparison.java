/*
 * Compares the sets of merge variants produced by two different runs of SV merging software
 */
import java.io.File;
import java.io.FileInputStream;
import java.util.Scanner;
import java.util.TreeSet;

public class MergeComparison {
	public static void main(String[] args) throws Exception
	{
		String fn1 = "/home/mkirsche/eichler/merged.vcf";
		String fn2 = "/home/mkirsche/eclipse-workspace/SVMergeEval/newmerged.vcf";
		TreeSet<MergedVariant> first = getAllMerged(fn1), second = getAllMerged(fn2);
		int[] counts = subsetCounts(first, second);
		System.out.printf("First only: %d (%.2f%% of first callset)\n",counts[1], 100.0 * counts[1] / first.size());
		System.out.printf("Second only: %d (%.2f%% of second callset)\n",counts[2], 100.0 * counts[2] / second.size());
		System.out.printf("Both: %d (%.4f%% of larger callset)\n",counts[3], 100.0 * counts[3] / Math.max(first.size(), second.size()));
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
	static TreeSet<MergedVariant> getAllMerged(String fn) throws Exception
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
			res.add(new MergedVariant(line));
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
		MergedVariant(String line) throws Exception
		{
			VcfEntry entry = new VcfEntry(line);
			
			// Check SUPP field for number of samples involved
			int supp = Integer.parseInt(entry.getInfo("SUPP"));
			ids = new String[supp];
			samples = new int[supp];
			
			// Get the indices of the samples from SUPP_VEC
			String suppVec = entry.getInfo("SUPP_VEC");
			int[] sampleIndices = new int[supp];
			int idx = 0;
			for(int i = 0; i<suppVec.length(); i++)
			{
				if(suppVec.charAt(i) == '1')
				{
					sampleIndices[idx++] = i;
				}
			}
			
			// Get the variant IDs
			String[] idlist = entry.getInfo("IDLIST").split(",");
			
			// Initialize all values
			for(int i = 0; i<supp; i++)
			{
				ids[i] = idlist[i];
				samples[i] = sampleIndices[i];
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
