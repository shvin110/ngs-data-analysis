package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ArrayBlockingQueue;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import structures.ListNum;
import structures.StringPair;
import var2.CallVariants2.Sample;

/**
 * Merges VCF files from multiple samples into a unified variant call set.
 * Creates a comprehensive variant list where every position that has a variant
 * in ANY input sample is represented in the output, with sample-specific
 * information preserved.
 * 
 * This is particularly useful for multi-sample variant calling workflows where
 * you want to ensure that interesting variants from any individual sample are
 * evaluated across all samples, even if they weren't initially called in every sample.
 * 
 * The merging process:
 * 1. Synchronously reads corresponding lines from all input VCF files
 * 2. Aggregates statistical evidence across samples for each position
 * 3. Preserves individual sample genotype information
 * 4. Combines header metadata appropriately
 * 
 * @author Brian Bushnell
 * @contributor Isla Winglet
 * @date December 18, 2016
 */
public class MergeSamples {
	
	/**
	 * Main method for command-line execution.
	 * 
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		MergeSamples x=new MergeSamples(args);
		//x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	/**
	 * Default constructor for programmatic use.
	 */
	public MergeSamples(){
		threads=Shared.threads();
		inq=new ArrayBlockingQueue<ListNum<VCFLine[]>>(threads+1);
	}
	
	/**
	 * Constructor that parses command-line arguments.
	 * 
	 * @param args Command line arguments array
	 */
	public MergeSamples(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("invalid")){
				outInvalid=b;
			}else if(a.equals("lines")){
				maxLines=Long.parseLong(b);
				if(maxLines<0){maxLines=Long.MAX_VALUE;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;
			
			in1=parser.in1;
			out1=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		threads=Shared.threads();
		inq=new ArrayBlockingQueue<ListNum<VCFLine[]>>(threads+1);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Core Methods          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Merges VCF files from a list of Sample objects (used by CallVariants2).
	 * This is the primary entry point for multi-sample variant calling workflows.
	 * 
	 * @param list List of Sample objects containing VCF file paths
	 * @param scafMap Scaffold mapping for coordinate resolution
	 * @param outVcf Output VCF filename
	 * @param scoreHistFile Optional score histogram output filename
	 */
	public void mergeSamples(ArrayList<Sample> list, ScafMap scafMap, String outVcf, String scoreHistFile){
		map=scafMap;
		ArrayList<StringPair> vcfList=new ArrayList<StringPair>(list.size());
		for(Sample s : list){vcfList.add(new StringPair(s.name, s.vcfName));}
		mergeFiles(vcfList, outVcf, scoreHistFile);
	}
	
	/**
	 * Core merging method that processes multiple VCF files synchronously.
	 * 
	 * The merging algorithm assumes that all input VCF files contain the same
	 * positions in the same order (typically generated by CallVariants2 with
	 * the same reference and coordinate system).
	 * 
	 * @param list List of (sample_name, vcf_filename) pairs
	 * @param outVcf Output merged VCF filename
	 * @param scoreHistFile Optional score histogram output filename
	 */
	public void mergeFiles(ArrayList<StringPair> list, String outVcf, String scoreHistFile){
		System.err.println("Merging "+list);
		final int ways=list.size();
		ByteFile[] bfa=new ByteFile[ways];
		final boolean allowSubprocess=(ways<=4);
		for(int i=0; i<ways; i++){
			StringPair pair=list.get(i);
			FileFormat ff=FileFormat.testInput(pair.b, FileFormat.VCF, null, allowSubprocess, false);
			bfa[i]=ByteFile.makeByteFile(ff);
		}
		
		mergeMT(outVcf, bfa);

		if(scoreHistFile!=null){
			CVOutputWriter.writeScoreHist(scoreHistFile, scoreArray);
		}
	}
	
	/**
	 * Single-threaded merging implementation (legacy).
	 * Processes files sequentially for simpler debugging but slower performance.
	 * 
	 * @param outVcf Output VCF filename
	 * @param bfa Array of ByteFile objects for input VCF files
	 */
	private void mergeST(String outVcf, ByteFile[] bfa){
		ByteStreamWriter bswVcf=null;
		if(outVcf!=null){
			bswVcf=new ByteStreamWriter(outVcf, true, false, true, FileFormat.VCF);
			bswVcf.start();
		}
		
		ByteBuilder bb=new ByteBuilder(34000);
		VCFLine[] row=processRow(bfa, bb);
		while(row!=null){
			if(row[0]!=null){
				VCFLine merged=merge(row);
				merged.toText(bb);
				bb.nl();
				if(bb.length>32000){
					if(bswVcf!=null){bswVcf.print(bb);}
					bb=new ByteBuilder(34000);
				}
			}
			row=processRow(bfa, bb);
		}
		
		if(bswVcf!=null){
			if(bb.length>0){bswVcf.print(bb);}
			bswVcf.poisonAndWait();
		}
	}
	
	/**
	 * Multithreaded merging implementation (current).
	 * Uses producer-consumer pattern for better performance on large datasets.
	 * 
	 * @param outVcf Output VCF filename
	 * @param bfa Array of ByteFile objects for input VCF files
	 */
	private void mergeMT(String outVcf, ByteFile[] bfa){
		ByteStreamWriter bswVcf=null;
		if(outVcf!=null){
			FileFormat ff=FileFormat.testOutput(outVcf, FileFormat.VCF, null, true, true, append, true);
			bswVcf=new ByteStreamWriter(ff);
			bswVcf.start();
		}
		
		ArrayList<MergeThread> alpt=spawnThreads(bswVcf);
		
		long nextID=0;
		ByteBuilder header=new ByteBuilder(34000);
		
		// Process header lines first
		VCFLine[] row=processRow(bfa, header);
		while(row!=null && row[0]==null){//Header
			row=processRow(bfa, header);
		}
		if(bswVcf!=null){
			bswVcf.add(header, nextID);
			nextID++;
		}
		
		// Process data lines in batches
		ListNum<VCFLine[]> list=new ListNum<VCFLine[]>(new ArrayList<VCFLine[]>(200), nextID);
		while(row!=null){
			if(row[0]!=null){
				list.add(row);
				if(list.size()>=200){
					putList(list);
					nextID++;
					list=new ListNum<VCFLine[]>(new ArrayList<VCFLine[]>(200), nextID);
				}
			}
			row=processRow(bfa, header);
		}
		if(list.size()>0){
			putList(list);
			nextID++;
		}
		
		putList(POISON_LIST);
		
		waitForFinish(alpt);
		
		if(bswVcf!=null){bswVcf.poisonAndWait();}
	}
	
	/**
	 * Reads one line from each input file and processes it.
	 * This method ensures that all files stay synchronized - each call processes
	 * the corresponding line from all input VCF files.
	 * 
	 * @param bfa Array of input ByteFile objects
	 * @param bb ByteBuilder for header accumulation
	 * @return Array of VCFLine objects (one per input file), or null if EOF
	 */
	VCFLine[] processRow(ByteFile[] bfa, ByteBuilder bb){
		byte[][] lines=new byte[bfa.length][];
		for(int i=0; i<bfa.length; i++){
			byte[] line=bfa[i].nextLine();
			if(line==null){return null;}
			lines[i]=line;
		}
		
		VCFLine[] row=new VCFLine[bfa.length];
		if(lines[0][0]=='#'){
			processHeader(lines, bb);
			return row;
		}
		for(int i=0; i<lines.length; i++){
			byte[] line=lines[i];
			row[i]=new VCFLine(line);
			if(i>0){assert(row[i].pos==row[0].pos) : "\n"+row[0]+"\n"+row[i];}
		}
		return row;
	}
	
	/**
	 * Processes header lines from all input files and merges them appropriately.
	 * Combines statistical metadata (read counts, quality averages) across samples
	 * while preserving format definitions and other header information.
	 * 
	 * @param lines Array of header lines (one from each input file)
	 * @param bb ByteBuilder for output accumulation
	 */
	void processHeader(byte[][] lines, ByteBuilder bb){
		String[][] matrix=new String[lines.length][];
		for(int i=0; i<lines.length; i++){
			matrix[i]=new String(lines[i]).split("=");
		}
		
		if(matrix[0][0].equals("##ploidy")){
			ploidy=Integer.parseInt(matrix[0][1]);
			bb.append("##ploidy="+ploidy+"\n");
		}else if(matrix[0][0].equals("##reads")){
			for(String[] split : matrix){
				reads+=Long.parseLong(split[1]);
			}
			bb.append("##reads="+reads+"\n");
		}else if(matrix[0][0].equals("##pairedReads")){
			for(String[] split : matrix){
				pairedReads+=Long.parseLong(split[1]);
			}
			bb.append("##pairedReads="+pairedReads+"\n");
		}else if(matrix[0][0].equals("##properlyPairedReads")){
			for(String[] split : matrix){
				properlyPairedReads+=Long.parseLong(split[1]);
			}
			properPairRate=properlyPairedReads*1.0/(Tools.max(1, reads));
			bb.append("##properlyPairedReads="+properlyPairedReads+"\n");
			bb.append("##properPairRate="+Tools.format("%.4f\n", properPairRate));
		}else if(matrix[0][0].equals("##properPairRate")){
			//do nothing - recalculated above
		}else if(matrix[0][0].equals("##totalQualityAvg")){
			totalQualityAvg=0;
			for(String[] split : matrix){
				totalQualityAvg+=Float.parseFloat(split[1]);
			}
			totalQualityAvg/=lines.length;
			bb.append("##totalQualityAvg="+Tools.format("%.4f\n", totalQualityAvg));
		}else if(matrix[0][0].equals("##mapqAvg")){
			mapqAvg=0;
			for(String[] split : matrix){
				mapqAvg+=Float.parseFloat(split[1]);
			}
			mapqAvg/=lines.length;
			bb.append("##mapqAvg="+Tools.format("%.2f\n", mapqAvg));
		}else if(matrix[0][0].equals("##readLengthAvg")){
			readLengthAvg=0;
			for(String[] split : matrix){
				readLengthAvg+=Float.parseFloat(split[1]);
			}
			readLengthAvg/=lines.length;
			bb.append("##readLengthAvg="+Tools.format("%.2f\n", readLengthAvg));
		}else if(matrix[0][0].startsWith("#CHROM\tPOS\t")){
			// Combine sample columns from all input files
			bb.append(lines[0]);
			for(int i=1; i<lines.length; i++){
				String[] split=new String(lines[i]).split("\t");
				bb.tab().append(split[split.length-1]);
			}
			bb.nl();
		}else{
			// Copy other header lines as-is from first file
			bb.append(lines[0]);
			bb.nl();
		}
	}
	
	/**
	 * Merges variant calls from multiple samples at the same genomic position.
	 * 
	 * The merging strategy:
	 * 1. Converts each VCFLine to a Var object for statistical aggregation
	 * 2. Sums all statistical evidence (read counts, quality scores, etc.)
	 * 3. Preserves individual sample genotype information
	 * 4. Uses the highest individual quality score as the merged quality
	 * 5. Regenerates VCF format with combined statistics
	 * 
	 * @param row Array of VCFLine objects (one per sample) at the same position
	 * @return Merged VCFLine representing the combined evidence
	 */
	VCFLine merge(VCFLine[] row){
		Var sum=null;
		VCFLine best=null;
		
		// Find best individual call and aggregate statistics
		for(VCFLine line : row){
			if(best==null || line.qual>best.qual){best=line;}
			Var v=line.toVar();
			assert(v!=null);
			if(sum==null){sum=v;}
			else{
				sum.add(v);        // Aggregate statistical evidence
				sum.addCoverage(v); // Combine coverage information
			}
		}
		assert(best!=null);
		assert(sum!=null) : row.length+", "+row[0];
		
		// Generate merged VCF line with combined statistics
		ByteBuilder bb=sum.toVCF(new ByteBuilder(), properPairRate, totalQualityAvg, mapqAvg, readLengthAvg, ploidy, map, filter, trimWhitespace);
		VCFLine merged=new VCFLine(bb.toBytes());
		
		// Preserve individual sample information
		merged.samples.clear();
		for(VCFLine line : row){
			merged.samples.addAll(line.samples);
		}
		
		// Use best individual quality if better than merged
		if(merged.qual<best.qual){
			merged.qual=best.qual;
			merged.filter=best.filter;
		}
		
		scoreArray[(int)merged.qual]++;
		return merged;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------     Threading Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Takes a batch of variant rows from the processing queue
	 * @return Batch to process */
	final ListNum<VCFLine[]> takeList(){
		ListNum<VCFLine[]> list=null;
		while(list==null){
			try {
				list=inq.take();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		return list;
	}
	
	/** Adds a batch of variant rows to the processing queue
	 * @param list Batch to queue */
	final void putList(ListNum<VCFLine[]> list){
		while(list!=null){
			try {
				inq.put(list);
				list=null;
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
	}
	
	/** Spawns worker threads for parallel merging
	 * @param bsw ByteStreamWriter for output
	 * @return List of spawned threads */
	private ArrayList<MergeThread> spawnThreads(ByteStreamWriter bsw){
		ArrayList<MergeThread> alpt=new ArrayList<MergeThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new MergeThread(bsw));
		}
		if(verbose){outstream.println("Spawned threads.");}
		
		for(MergeThread pt : alpt){
			pt.start();
		}
		if(verbose){outstream.println("Started threads.");}
		
		return alpt;
	}
	
	/** Waits for all worker threads to complete
	 * @param alpt List of threads to wait for */
	private void waitForFinish(ArrayList<MergeThread> alpt){
		boolean allSuccess=true;
		for(MergeThread pt : alpt){
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					pt.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	/**
	 * Worker thread that processes batches of variant rows.
	 * Each thread merges variants independently and sends output to ByteStreamWriter.
	 */
	private class MergeThread extends Thread {

		/** Creates a merge thread
		 * @param bsw_ ByteStreamWriter for output */
		MergeThread(ByteStreamWriter bsw_){
			bsw=bsw_;
		}

		/** Main thread execution loop */
		@Override
		public void run(){
			ListNum<VCFLine[]> list=takeList();
			while(list!=null && list!=POISON_LIST){
				processList(list);
				list=takeList();
			}
			putList(POISON_LIST);
		}

		/** Processes a batch of variant rows
		 * @param list Batch of rows to merge */
		private void processList(ListNum<VCFLine[]> list){
			ByteBuilder bb=new ByteBuilder(4096);
			for(VCFLine[] row : list){
				mergeRow(row, bb);
			}
			if(bsw!=null){bsw.add(bb, list.id);}
		}
		
		/** Merges a single row of variants across samples
		 * @param row Array of VCFLines at same position
		 * @param bb ByteBuilder for output accumulation */
		private void mergeRow(VCFLine[] row, ByteBuilder bb){
			if(row[0]!=null){
				VCFLine merged=merge(row);
				merged.toText(bb);
				bb.nl();
			}
		}
		
		/** ByteStreamWriter for ordered output */
		private final ByteStreamWriter bsw;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Poison pill for ending thread processing */
	final ListNum<VCFLine[]> POISON_LIST=new ListNum<VCFLine[]>(null, -1);
	/** Queue for passing variant batches between threads */
	private final ArrayBlockingQueue<ListNum<VCFLine[]>> inq;
	/** Number of worker threads */
	private final int threads;
	
	/** Sum of reads across all samples */
	long readsSum;
	/** Sum of pairs across all samples */
	long pairsSum;
	/** Sample ploidy */
	int ploidy=1;
	
	/** Combined proper pair rate */
	double properPairRate;
	/** Average total quality across samples */
	double totalQualityAvg;
	/** Average mapping quality across samples */
	double mapqAvg;
	/** Average read length across samples */
	double readLengthAvg;
	
	/** Total reads processed */
	long reads;
	/** Total paired reads */
	long pairedReads;
	/** Total properly paired reads */
	long properlyPairedReads;
	
	/** Filtering criteria */
	VarFilter filter;
	/** Scaffold mapping */
	ScafMap map;
	/** Whether to trim whitespace from scaffold names */
	boolean trimWhitespace=true;
	
	/** Primary input filename */
	private String in1=null;
	/** Primary output filename */
	private String out1=null;
	/** Invalid variants output filename */
	private String outInvalid=null;
	
	/** Score histogram for quality distribution analysis */
	long[] scoreArray=new long[200];
	
	/** Lines processed counter */
	private long linesProcessed=0;
	/** Valid lines counter */
	private long linesValid=0;
	/** Bytes processed counter */
	private long bytesProcessed=0;
	/** Maximum lines to process */
	private long maxLines=Long.MAX_VALUE;
	
	/** Output stream for messages */
	private PrintStream outstream=System.err;
	/** Verbose output flag */
	public static boolean verbose=false;
	/** Error state flag */
	public boolean errorState=false;
	/** Overwrite output files flag */
	private boolean overwrite=true;
	/** Append to output files flag */
	private boolean append=false;
}