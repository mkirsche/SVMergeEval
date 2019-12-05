
public class RunMergingSoftware 
{
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
