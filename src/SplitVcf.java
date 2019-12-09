/*
 * Program for splitting a merged VCF - introduces noise to see if
 * merging software can reconstruct the original merged variants
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;
 
public class SplitVcf
{
	// The probability of each base randomly mutating
	static double SEQ_ERROR = 0.1;
	
	// The average amount the breakpoints move from their merged value when splitting
	static int AVERAGE_POS_ERROR = 50;
	
	// Random number generator
	static Random rand;
	public static void main(String[] args) throws Exception
	{
		rand = new Random(31);
		String inputFile = "/home/mkirsche/eichlersim/surv.vcf";
		String outDir = "eichler_split_surv";
		if(args.length == 2)
		{
			inputFile = args[0];
			outDir = args[1];
		}
		else
		{
			System.out.println("Usage: java SplitVcf <mergedVcf> <outputDirectory>");
		}
		splitVcf(inputFile, outDir);
	}
	
	/* 
	 * Splits a merged VCF into individual files with noise
	 */
	static void splitVcf(String inputFile, String outDir) throws Exception
	{
		int n = countInputs(inputFile);	
		if(n == 0)
		{
			return;
		}
		
		// Make output directory for the split VCF files
		Path currentRelativePath = Paths.get("");
		outDir = currentRelativePath.toAbsolutePath().toString() + "/" + outDir;
		File f = new File(outDir);
		f.mkdir();
		
		String[] fns = new String[n];
		for(int i = 0; i<n; i++) fns[i] = outDir + "/" + "sample" + (i+1) + ".vcf";
		
		PrintWriter[] outs = new PrintWriter[n];
		for(int i = 0; i<n; i++)
		{
			outs[i] = new PrintWriter(new File(fns[i]));
		}
		
		// Write split VCF filenames to filelist.txt
		PrintWriter out = new PrintWriter(new File(outDir + "/filelist.txt"));
		for(int i = 0; i<n; i++)
		{
			out.println(fns[i]);
		}
		out.close();
		
		Scanner input = new Scanner(new FileInputStream(new File(inputFile)));
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0)
			{
				continue;
			}
			if(line.startsWith("#"))
			{
				for(int i = 0; i<n; i++)
				{
					outs[i].println(line);
				}
			}
			else
			{
				VcfEntry entry = new VcfEntry(line);
				String suppVec = entry.getInfo("SUPP_VEC");
				String method = entry.getInfo("SVMETHOD");
				int supp = Integer.parseInt(entry.getInfo("SUPP"));
				
				// Get the list of IDs - the way of doing that depends on the method used to merge them
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
						if(!val.equalsIgnoreCase("nan")) ids.add(val);
					}
					for(int i = 0; i<ids.size(); i++)
					{
						idList[i] = ids.get(i);
					}
				}
				
				int idIdx = 0;
				for(int i = 0; i<n; i++)
				{
					boolean inSample = suppVec.charAt(i) == '1';
					if(inSample)
					{
						VcfEntry newEntry = new VcfEntry(entry.originalLine);
						
						newEntry.setId(idList[idIdx]);
						
						long oldStart = newEntry.getPos();
						
						int curNoise = Math.min(AVERAGE_POS_ERROR, (int)(Math.abs(newEntry.getLength()) * 0.25));
						
						// Add some noise to different parts of it
						int dPos = randomDelta(curNoise);
						int dLen = randomDelta(curNoise);
						
						long newStart = oldStart + dPos;
						int newLen = newEntry.getLength() + dLen;
						long newEnd = newStart + newLen;
						
						// Set start
						newEntry.setPos(newEntry.getPos() + dPos);
						
						// Set SVLEN INFO field
						if(newEntry.hasInfoField("SVLEN"))
						{
							newEntry.setInfo("SVLEN", newLen + "");
						}
						
						String type = newEntry.getType();
						int refLength = newEntry.getRef().length();
						int altLength = newEntry.getAlt().length();
						if(!newEntry.getAlt().contains("<") && !newEntry.getRef().contains("N"))
						{
							if(refLength > altLength)
							{
								type = "DEL";
							}
							if(refLength < altLength)
							{
								type = "INS";
							}
						}
						
						// Set END INFO field
						if(newEntry.hasInfoField("END"))
						{
							if(type.equals("INS"))
							{
								newEntry.setInfo("END", newStart + "");
							}
							else
							{
								newEntry.setInfo("END", newEnd + "");
							}
						}
						
						// Set REF and ALT fields depending on SV type
						if(!newEntry.getAlt().contains("<") && !newEntry.getRef().contains("N"))
						{
							if(type.equals("INS"))
							{
								if(entry.getLength() < 0)
								{
									System.out.println(entry);
								}
								entry.setAlt(adjustStringLength(entry.getAlt(), newLen + 1));
								if(entry.hasInfoField("SEQ"))
								{
									entry.setInfo("SEQ", entry.getAlt().substring(1));
								}
							}
							else if(type.equals("DEL"))
							{
								entry.setRef(adjustStringLength(entry.getRef(), Math.abs(newLen) + 1));
								if(entry.hasInfoField("SEQ"))
								{
									entry.setInfo("SEQ", entry.getRef().substring(1));
								}
							}
							else if(type.equals("DEL/INV"))
							{
								entry.setRef(adjustStringLength(entry.getRef(), Math.abs(newLen) + 1));
								if(entry.hasInfoField("SEQ"))
								{
									entry.setInfo("SEQ", entry.getRef().substring(1));
								}
							}
							else
							{
								System.out.println(type);
								System.out.println(entry.getRef()+" "+entry.getAlt());
							}
						}
						
						outs[i].println(newEntry);
						idIdx++;
					}
				}
			}
		}
		for(int i = 0; i<n; i++)
		{
			outs[i].close();
		}
		input.close();
	}
	
	/*
	 * Adjusts a string to incorporate errors and changes in its length
	 */
	static String adjustStringLength(String s, int newLength)
	{
		char[] res = new char[newLength];
		res[0] = s.charAt(0);
		
		if(newLength == s.length())
		{
			// In this case no insertions or deletions are needed
			for(int i = 1; i<res.length; i++)
			{
				res[i] = s.charAt(i);
			}
		}
		else if(newLength < s.length())
		{
			// Deletions needed so pick random subsequence of indices to keep
			int[] subseq = getRandomSubsequence(1, s.length() - 1, newLength - 1);
			for(int i = 1; i<res.length; i++)
			{
				res[i] = s.charAt(subseq[i-1]);
			}
		}
		else
		{
			// Insert random bases in gaps
			for(int i = 1; i<res.length; i++)
			{
				res[i] = intToBase(rand.nextInt(4));
			}
			
			// Randomly select how the bases from the original sequence get placed
			int[] subseq = getRandomSubsequence(1, newLength - 1, s.length() - 1);
			for(int i = 0; i<subseq.length; i++)
			{
				res[subseq[i]] = s.charAt(i+1);
			}
		}
		
		// Mutate each base with a fixed probability
		for(int i = 1; i<res.length; i++)
		{
			res[i] = mutateBaseWithProb(res[i], SEQ_ERROR);
		}
		
		return new String(res);
	}
	
	/*
	 * Mutates a given base to a random other base with a fixed probability
	 */
	static char mutateBaseWithProb(char base, double p)
	{
		if(rand.nextDouble() < p)
		{
			int baseAsInt = baseToInt(base);
			char newBase = intToBase((baseAsInt + rand.nextInt(3) + 1)%4);
			if(Character.isLowerCase(base))
			{
				newBase += 'a' - 'A';
			}
			return newBase;
		}
		else
		{
			return base;
		}
	}
	
	/*
	 * Converts integer in [0, 3] to the corresponding element of [A, C, G, T]
	 */
	static char intToBase(int val)
	{
		if(val == 0)
		{
			return 'A';
		}
		if(val == 1)
		{
			return 'C';
		}
		if(val == 2)
		{
			return 'G';
		}
		return 'T';
	}
	
	/*
	 * Converts basepair characters to an index (A = 0, C = 1, G = 2, T = 3)
	 */
	static int baseToInt(char c)
	{
		if(c == 'a' || c == 'A')
		{
			return 0;
		}
		if(c == 'c' || c == 'C')
		{
			return 1;
		}
		if(c == 'g' || c == 'G')
		{
			return 2;
		}
		return 3;
	}
	
	/*
	 * Gets a random subsequence of a given length of integers in the range [minVal, maxVal]
	 */
	static int[] getRandomSubsequence(int minVal, int maxVal, int subseqLength)
	{
		ArrayList<Integer> vals = new ArrayList<Integer>();
		for(int i = minVal; i<=maxVal; i++)
		{
			vals.add(i);
		}
		Collections.shuffle(vals);
		int[] res = new int[subseqLength];
		
		for(int i = 0; i<res.length; i++)
		{
			res[i] = vals.get(i);
		}
		
		Arrays.sort(res);
		return res;
	}
	
	/*
	 * Gets a uniform random number between -2*average and 2*average, inclusive
	 */
	static int randomDelta(int average)
	{
		return rand.nextInt(4*average + 1) - 2*average;
	}
	
	/*
	 * Counts how many input files were used to produced a given merged VCF File based on length of SUPP_VEC field
	 */
	static int countInputs(String mergedVcf) throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(mergedVcf)));
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				continue;
			}
			
			VcfEntry entry = new VcfEntry(line);
			input.close();
			return entry.getInfo("SUPP_VEC").length();
		}
		input.close();
		return 0;
	}
	
}
