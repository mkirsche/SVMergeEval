/*
 * Methods for storing VCF v4.2 entries for structural variants
 * It allows easy access to main variant fields, plus some methods
 * to do more parsing and error-checking for things like type and length.
 */

import java.util.Arrays;

public class VcfEntry {

	String originalLine;
	String[] tabTokens;
	
	public VcfEntry(String line) throws Exception
	{
		originalLine = line;
		tabTokens = line.split("\t");
		if(tabTokens.length < 8)
		{
			throw new Exception("VCF line had too few entries: "
					+ Arrays.toString(tabTokens));
		}
	}
	
	/*
	 * Reconstruct the VCF line by concatenating and tab-separating the fields
	 */
	public String toString()
	{
		StringBuilder sb = new StringBuilder("");
		for(int i = 0; i<tabTokens.length; i++)
		{
			sb.append(tabTokens[i]);
			if(i < tabTokens.length - 1)
			{
				sb.append("\t");
			}
		}
		return sb.toString();
	}
	
	/*
	 * Get the chromosome field
	 */
	public String getChromosome()
	{
		return tabTokens[0];
	}
	
	/*
	 * Set the chromosome field
	 */
	public void setChromosome(String s)
	{
		tabTokens[0] = s;
	}
	
	/*
	 * Get hte POS field
	 */
	public long getPos() throws Exception
	{
		return Long.parseLong(tabTokens[1]);
	}
	
	/*
	 * Set the POS field
	 */
	public void setPos(long val)
	{
		tabTokens[1] = val+"";
	}
	
	/*
	 * Get the variant ID field
	 */
	public String getId()
	{
		return tabTokens[2];
	}
	
	/*
	 * Set the variant ID field
	 */
	public void setId(String s)
	{
		tabTokens[2] = s;
	}
	
	/*
	 * Get the REF sequence field
	 */
	public String getRef()
	{
		return tabTokens[3];
	}
	
	/*
	 * Set the REF sequence field
	 */
	public void setRef(String s)
	{
		tabTokens[3] = s;
	}
	
	/*
	 * Get the ALT sequence field
	 */
	public String getAlt()
	{
		return tabTokens[4];
	}
	
	/*
	 * Set the ALT sequence field
	 */
	public void setAlt(String s)
	{
		tabTokens[4] = s;
	}
	
	/*
	 * The length of a variant
	 */
	public int getLength() throws Exception
	{
		try {
			String s = getInfo("SVLEN");
			return (int)(.5 + Double.parseDouble(s));
		} catch(Exception e) {
            String seq = getSeq();
            String type = getType();
            if(type.equals("INS")) return seq.length();
            else return -seq.length();
		}
	}
	
	/*
	 * The end position of a variant
	 */
	public long getEnd() throws Exception
	{
		if(hasInfoField("END")) return Long.parseLong(getInfo("END"));
		String type = getType();
		if(type.equals("INS"))
		{
			return getPos();
		}
		else
		{
			return getPos() + Math.abs(getLength());
		}
	}
	
	/*
	 * Get the variant type
	 */
	public String getType() throws Exception
	{
		String res = getInfo("SVTYPE");
		if(res.length() == 0)
		{
			int refLength = getRef().length();
			int altLength = getAlt().length();
			if(getAlt().startsWith("<") && getAlt().endsWith(">"))
			{
				return getAlt().substring(1, getAlt().length() - 1);
			}
			if(refLength > altLength)
			{
				return "INS";
			}
			else if(refLength < altLength)
			{
				return "DEL";
			}
			else
			{
				return "";
			}
		}
		else return res;
	}
	
	/*
	 * Set the SVTYPE INFO field to a certain value
	 */
	public void setType(String s) throws Exception
	{
		setInfo("SVTYPE", s);
	}
	
	/*
	 * The strands on which the variant occurs
	 */
	public String getStrand() throws Exception
	{
		return getInfo("STRANDS");
	}
	
	/*
	 * The genomic sequence being affected/added/deleted by the variant
	 */
	public String getSeq() throws Exception
	{
		if(hasInfoField("SEQ"))
		{
			return getInfo("SEQ");
		}
		String ref = getRef(), alt = getAlt();
		if(alt.startsWith("<"))
		{
			return "";
		}
		String type = getType();
		
		// If the SV is a deletion, swap REF and ALT and treat as an insertion
		if(type.equals("DEL"))
		{
			String tmp = ref;
			ref = alt;
			alt = tmp;
			type = "INS";
		}
		if(ref.equals("X") || ref.equals("N"))
		{
			return alt;
		}
		else if(alt.equals("X") || ref.equals("N"))
		{
			return ref;
		}
		else if(type.equals("INS"))
		{
			int startPad = 0, endPad = 0;
			int totalPad = ref.length();
			while(startPad + endPad < totalPad && startPad < alt.length())
			{
				if(ref.charAt(startPad) == alt.charAt(startPad))
				{
					startPad++;
				}
				else if(ref.charAt(ref.length() - 1 - endPad) == alt.charAt(alt.length() - 1 - endPad))
				{
					endPad++;
				}
				else
				{
					break;
				}
			}
			if(startPad + endPad == totalPad)
			{
				return alt.substring(startPad, alt.length() - endPad);
			}
			else
			{
				return alt;
			}
		}
		else
		{
			return alt;
		}
	}
	
	/*
	 * Set a particular VCF INFO field, adding the field if it doesn't already exist
	 */
	public void setInfo(String field, String val) throws Exception
	{
		String[] infoFields = tabTokens[7].split(";");
		for(String semitoken : infoFields)
		{
			int equalIndex = semitoken.indexOf('=');
			if(equalIndex == -1)
			{
				continue;
			}
			String key = semitoken.substring(0, equalIndex);
			if(key.equals(field))
			{
				String updatedToken = key + "=" + val;
				
				// Special case if this is the first INFO field
				if(tabTokens[7].startsWith(semitoken))
				{
					tabTokens[7] = tabTokens[7].replaceFirst(semitoken, updatedToken);
				}
				else
				{
					tabTokens[7] = tabTokens[7].replaceAll(";" + semitoken, ";" + updatedToken);
				}
				return;
			}
		}
		
		// Field not found, so add it!
		tabTokens[7] += ";" + field + "=" + val;
	}
	
	/*
	 * Get list of supporting read names - first check for RNAMES field, and then anything containing RNAMES
	 */
	public String[] getRnames() throws Exception
	{
		if(hasInfoField("RNAMES"))
		{
			return getInfo("RNAMES").split(",");
		}
		String infoToken = tabTokens[7];
		String[] semicolonSplit = infoToken.split(";");
		for(String semitoken : semicolonSplit)
		{
			int equalIndex = semitoken.indexOf('=');
			if(equalIndex == -1)
			{
				continue;
			}
			String key = semitoken.substring(0, equalIndex);
			if(key.toUpperCase().contains("RNAMES"))
			{
				return semitoken.substring(1 + equalIndex).split(",");
			}
		}
		return new String[] {};
	}
	
	/*
	 * Get the value of a particular INFO field
	 */
	public String getInfo(String field) throws Exception
	{
		String infoToken = tabTokens[7];
		String[] semicolonSplit = infoToken.split(";");
		for(String semitoken : semicolonSplit)
		{
			int equalIndex = semitoken.indexOf('=');
			if(equalIndex == -1)
			{
				continue;
			}
			String key = semitoken.substring(0, equalIndex);
			if(key.equals(field))
			{
				return semitoken.substring(1 + equalIndex);
			}
		}
		return "";
	}
	
	/*
	 * Whether or not the line has a particular INFO field
	 */
	public boolean hasInfoField(String fieldName)
	{
		String infoToken = tabTokens[7];
		String[] semicolonSplit = infoToken.split(";");
		for(String semitoken : semicolonSplit)
		{
			int equalIndex = semitoken.indexOf('=');
			if(equalIndex == -1)
			{
				continue;
			}
			String key = semitoken.substring(0, equalIndex);
			if(key.equals(fieldName))
			{
				return true;
			}
		}
		return false;
	}

}
