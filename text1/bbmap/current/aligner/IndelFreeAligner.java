package aligner;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.SIMDAlignByte;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import stream.SamLine;
import structures.ByteBuilder;
import structures.IntHashMap;
import structures.IntList;
import structures.IntListHashMap;
import structures.ListNum;
import template.Accumulator;
import template.ThreadWaiter;
import tracker.ReadStats;

/**
 * Performs indel-free alignments of queries.
 * Designed for a small query set with stays in memory,
 * and a large reference set which is streamed.
 * Allows an arbitrary number of substitutions.
 * Uses SIMD.
 * 
 * @author Brian Bushnell
 * @contributor Isla
 * @date June 2, 2025
 *
 */
public class IndelFreeAligner implements Accumulator<IndelFreeAligner.ProcessThread> {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();

		//Create an instance of this class
		IndelFreeAligner x=new IndelFreeAligner(args);

		//Run the object
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public IndelFreeAligner(String[] args){

		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}

		//Set shared static variables prior to parsing
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());

		{//Parse the arguments
			final Parser parser=parse(args);
			Parser.processQuality();

			maxReads=parser.maxReads;
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			in1=parser.in1;
			in2=parser.in2;
			extin=parser.extin;

			out1=parser.out1;
			extout=parser.extout;
		}

		validateParams();
		doPoundReplacement(); //Replace # with 1 and 2
		fixExtensions(); //Add or remove .gz or .bz2 as needed
		checkFileExistence(); //Ensure files can be read and written
		checkStatics(); //Adjust file-related static fields as needed for this program 

		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);

		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	/** Parse arguments from the command line */
	private Parser parse(String[] args){

		//Create a parser object
		Parser parser=new Parser();

		//Set any necessary Parser defaults here
		//parser.foo=bar;

		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];

			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("ref")){
				refFile=b;
			}else if(a.equals("subs") || a.equals("maxsubs")){
				maxSubs=Integer.parseInt(b);
			}else if(a.equals("hits") || a.equals("minhits") || a.equals("seedhits")){
				minSeedHits=Math.max(1, Integer.parseInt(b));
			}else if(a.equals("minprob") || a.equals("minhitsprob")){
				minHitsProb=Float.parseFloat(b);
			}else if(a.equals("maxclip") || a.equals("clip")){
				Query.maxClip=Tools.max(0, Float.parseFloat(b));
				assert(Query.maxClip<1 || Query.maxClip==(int)Query.maxClip);
			}else if(a.equals("index")){
				indexQueries=Parse.parseBoolean(b);
			}else if(a.equals("prescan")){
				prescan=Parse.parseBoolean(b);
			}else if(a.equals("seedmap") || a.equals("map")){
				useSeedMap=Parse.parseBoolean(b);
			}else if(a.equals("seedlist") || a.equals("list")){
				useSeedMap=!Parse.parseBoolean(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
				assert(k<16) : "0<=k<16 : "+k;
				indexQueries=(k>0);
			}else if(a.equals("qstep") || a.equals("step") || a.equals("qskip")){
				qStep=Integer.parseInt(b);
				assert(qStep>0);
			}else if(a.equals("mm")){
				midMaskLen=(Tools.isNumeric(b) ? Integer.parseInt(b) :
					Parse.parseBoolean(b) ? 1 : 0);
			}else if(a.equals("blacklist") || a.equals("banhomopolymers")){
				Query.blacklistRepeatLength=(Tools.isNumeric(b) ? Integer.parseInt(b) :
					Parse.parseBoolean(b) ? 1 : 0);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Parse.parseKMG(b);
				//Set a variable here
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		return parser;
	}

	/** Replace # with 1 and 2 in headers */
	private void doPoundReplacement(){
		//Do input file # replacement
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}

		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
	}

	/** Add or remove .gz or .bz2 as needed */
	private void fixExtensions(){
		in1=Tools.fixExtension(in1);
		in2=Tools.fixExtension(in2);
	}

	/** Ensure files can be read and written */
	private void checkFileExistence(){
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1, in2)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, in2, out1)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
	}

	/** Adjust file-related static fields as needed for this program */
	private static void checkStatics(){
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		assert(FastaReadInputStream.settingsOK());
	}

	/** Ensure parameter ranges are within bounds and required parameters are set */
	private boolean validateParams(){
		assert((k>=1 && k<=15) || !indexQueries);
		assert(minHitsProb<=1);
		assert(midMaskLen<k-1);
		assert(maxSubs>=0);
		return ((k>=1 && k<=15) || !indexQueries) && 
				(minHitsProb<=1) && (midMaskLen<k-1) && (maxSubs>=0);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){

		//Reset counters
		readsProcessed=readsOut=0;
		basesProcessed=basesOut=0;

		SamLine.RNAME_AS_BYTES=false;

		Query.setMode(k, midMaskLen, indexQueries);
		final ArrayList<Query> queries=fetchQueries(ffin1, ffin2);

		//Turn off read validation in the input threads to increase speed
		final boolean vic=Read.VALIDATE_IN_CONSTRUCTOR;
		Read.VALIDATE_IN_CONSTRUCTOR=Shared.threads()<4;

		//Create a read input stream
		final ConcurrentReadInputStream cris=makeCris(refFile);

		//Optionally create a read output stream
		final ByteStreamWriter bsw=ByteStreamWriter.makeBSW(ffout1);

		//Process the reads in separate threads
		spawnThreads(cris, bsw, queries);

		if(verbose){outstream.println("Finished; closing streams.");}

		//Write anything that was accumulated by ReadStats
		errorState|=ReadStats.writeAll();
		//Close the read streams
		errorState|=ReadWrite.closeStreams(cris);

		if(bsw!=null) {bsw.poisonAndWait();}

		//Reset read validation
		Read.VALIDATE_IN_CONSTRUCTOR=vic;

		//Report timing and results
		t.stop();
		outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
		outstream.println(Tools.readsBasesOut(readsProcessed, basesProcessed, readsOut, basesOut, 8, false));
		outstream.println(Tools.things("Alignments", alignmentCount, 8));
		outstream.println(Tools.things("Seed Hits", seedHitCount, 8));

		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	/** Spawn process threads */
	private void spawnThreads(final ConcurrentReadInputStream cris, 
			final ByteStreamWriter bsw, final ArrayList<Query> qList){

		//Do anything necessary prior to processing

		//Determine how many threads may be used
		final int threads=Shared.threads();

		//Fill a list with ProcessThreads
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(threads);
		for(int i=0; i<threads; i++){
			alpt.add(new ProcessThread(cris, bsw, qList, maxSubs, i));
		}

		//Start the threads and wait for them to finish
		boolean success=ThreadWaiter.startAndWait(alpt, this);
		errorState&=!success;

		//Do anything necessary after processing

	}

	@Override
	public final void accumulate(ProcessThread pt){
		synchronized(pt) {
			readsProcessed+=pt.readsProcessedT;
			basesProcessed+=pt.basesProcessedT;
			alignmentCount+=pt.alignmentsT;
			seedHitCount+=pt.seedHitsT;

			readsOut+=pt.readsOutT;
			basesOut+=pt.basesOutT;
			errorState|=(!pt.success);
		}
	}

	@Override
	public final boolean success(){return !errorState;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Creates a ConcurrentReadInputStream for the reference file.
	 * @param fname Path to the reference file
	 * @return Initialized and started input stream
	 */
	private ConcurrentReadInputStream makeCris(String fname){
		FileFormat ff=FileFormat.testInput(fname, null, true);
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff, null);
		cris.start(); // Start the stream
		if(verbose){outstream.println("Started cris");}
		return cris;
	}

	/**
	 * Loads query sequences from input files and converts them to Query objects.
	 * Pre-computes k-mer indices and alignment parameters for each query.
	 * @param ff1 Primary input file format
	 * @param ff2 Secondary input file format (may be null)
	 * @return List of Query objects ready for alignment
	 */
	public ArrayList<Query> fetchQueries(FileFormat ff1, FileFormat ff2){
		Timer t=new Timer(outstream, false);
		ArrayList<Read> reads=ConcurrentReadInputStream.getReads(maxReads, false, ff1, ff2, null, null);
		ArrayList<Query> queries=new ArrayList<Query>(reads.size());
		if(indexQueries) {
			Query.mhc=new MinHitsCalculator(k, maxSubs, midMaskLen, minHitsProb, Query.maxClip); // Initialize hit calculator
		}
		for(Read r : reads) {//TODO: Could be multithreaded.
			readsProcessed+=r.pairCount();
			basesProcessed+=r.pairLength();
			Query query=new Query(r.id, queries.size(), r.bases, r.quality);
			queries.add(query);
			if(r.mate==null) {continue;}
			r=r.mate; // Process mate if present
			query=new Query(r.id, queries.size(), r.bases, r.quality);
			queries.add(query);
		}
		t.stop("Loaded "+queries.size()+" queries in ");
		return queries;
	}

	/**
	 * Performs sparse alignment using seed hits to guide alignment positions.
	 * More efficient than brute force when seed hits are available and selective.
	 * @param query Query sequence bases
	 * @param ref Reference sequence bases  
	 * @param maxSubs Maximum allowed substitutions
	 * @param maxClips Maximum allowed clipped bases
	 * @param seedHits List of potential alignment start positions from seed matching
	 * @return List of alignment start positions with ≤maxSubs substitutions, or null if none found
	 */
	public static IntList alignSparse(byte[] query, byte[] ref, int maxSubs, int maxClips, IntList seedHits) {
		if(seedHits==null || seedHits.isEmpty()) {
			return null; // No seed guidance available
		}

		IntList results=null;

		for(int i=0; i<seedHits.size; i++) {
			int rStart=seedHits.array[i];
			int subs;

			// Choose appropriate alignment method based on position
			if(rStart<0) {
				subs=alignClipped(query, ref, maxSubs, maxClips, rStart); // Left overhang
			} else if(rStart>ref.length-query.length) {
				subs=alignClipped(query, ref, maxSubs, maxClips, rStart); // Right overhang
			} else {
				subs=align(query, ref, maxSubs, rStart); // Perfect fit within reference
			}

			if(subs<=maxSubs) {
				if(results==null) {results=new IntList(4);}
				results.add(rStart);
			}
		}

		return results;
	}

	/**
	 * Performs comprehensive alignment testing all possible positions.
	 * Uses SIMD optimization when available and appropriate.
	 * @param query Query sequence bases
	 * @param ref Reference sequence bases
	 * @param maxSubs Maximum allowed substitutions
	 * @param maxClips Maximum allowed clipped bases
	 * @return List of all alignment positions with ≤maxSubs substitutions, or null if none found
	 */
	public static IntList alignAllPositions(byte[] query, byte[] ref, int maxSubs, int maxClips) {
		if(Shared.SIMD && (query.length<256 || maxSubs<256)) {
//			return SIMDAlignByte.alignDiagonal(query, ref, maxSubs); // Use vectorized alignment
			return SIMDAlignByte.alignDiagonal(query, ref, maxSubs, maxClips); // Use vectorized alignment
			//TODO: Pad short contigs to avoid scalar mode
		}
		IntList list=null;
		int rStart=-maxSubs;
		
		// Left overhang region (negative start positions)
		for(; rStart<0; rStart++) {
			int subs=alignClipped(query, ref, maxSubs, maxClips, rStart);
			if(subs<=maxSubs) {
				if(list==null) {list=new IntList(4);}
				list.add(rStart);
			}
		}
		
		// Perfect fit region (query completely within reference)
		for(final int limit=ref.length-query.length; rStart<=limit; rStart++) {
			int subs=align(query, ref, maxSubs, rStart);
			if(subs<=maxSubs) {
				if(list==null) {list=new IntList(4);}
				list.add(rStart);
			}
		}
		
		// Right overhang region (query extends past reference end)
		for(final int limit=ref.length-query.length+maxSubs; rStart<=limit; rStart++) {
			int subs=alignClipped(query, ref, maxSubs, maxClips, rStart);
			if(subs<=maxSubs) {
				if(list==null) {list=new IntList(4);}
				list.add(rStart);
			}
		}
		return list;
	}

	/**
	 * Aligns query to reference starting at specified position with no clipping.
	 * Optimized for cases where query fits completely within reference bounds.
	 * @param query Query sequence bases
	 * @param ref Reference sequence bases
	 * @param maxSubs Maximum allowed substitutions (for early termination)
	 * @param rStart Starting position in reference (0-based)
	 * @return Number of substitutions found
	 */
	static int align(byte[] query, byte[] ref, final int maxSubs, final int rStart) {
		int subs=0;
		for(int i=0, j=rStart; i<query.length && subs<=maxSubs; i++, j++) {
			final byte q=query[i], r=ref[j];
			final int incr=(q!=r || AminoAcid.baseToNumber[q]<0 ? 1 : 0); // Count mismatches and N's
			subs+=incr;
		}
		return subs;
	}

	/**
	 * Aligns query to reference with clipping support for overhangs.
	 * Handles cases where query extends beyond reference boundaries.
	 * @param query Query sequence bases
	 * @param ref Reference sequence bases
	 * @param maxSubs Maximum allowed substitutions
	 * @param maxClips Maximum allowed clipped bases
	 * @param rStart Starting position in reference (may be negative)
	 * @return Number of substitutions found (including excess clipping as substitutions)
	 */
	static int alignClipped(byte[] query, byte[] ref, int maxSubs, final int maxClips, 
			final int rStart) {
		final int rStop1=rStart+query.length; // Position after final base
		final int leftClip=Math.max(0, -rStart), rightClip=Math.max(0, rStop1-ref.length);
		int clips=leftClip+rightClip;
		if(clips>=query.length) {return query.length;} // Entirely clipped
		int subs=Math.max(0, clips-maxClips); // Excess clipping counts as substitutions
		int i=leftClip, j=rStart+leftClip; // Skip clipped bases
		
		// Align overlapping region
		for(final int limit=Math.min(rStop1, ref.length); j<limit && subs<=maxSubs; i++, j++) {
			final byte q=query[i], r=ref[j];
			final int incr=(q!=r || AminoAcid.baseToNumber[q]<0 ? 1 : 0);
			subs+=incr;
		}
		return subs;
	}

	/**
	 * Builds a k-mer index for a reference sequence.
	 * Maps masked k-mers to their positions for efficient seed finding.
	 * @param ref Reference sequence bases
	 * @return Hash map from masked k-mers to lists of positions, or null if indexing disabled
	 */
	IntListHashMap buildReferenceIndex(byte[] ref) {
		if(!indexQueries || k<=0) {return null;}

		final int defined=Math.max(k-maxSubs, 2);
		final int kSpace=(1<<(2*defined));
		final long maxKmers=Math.min(kSpace, (ref.length-k+1)*2L);
		final int initialSize=(int)Math.min(4000000, ((maxKmers*3)/2));
		final IntListHashMap index=new IntListHashMap(initialSize);

		final int shift=2*k, shift2=shift-2, mask=~((-1)<<shift); // Bit manipulation constants
		int kmer=0, rkmer=0, len=0; // Rolling k-mer state

		for(int i=0; i<ref.length; i++) {
			final byte b=ref[i];
			final int x=AminoAcid.baseToNumber[b], x2=AminoAcid.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask; // Roll forward k-mer
			rkmer=((rkmer>>>2)|(x2<<shift2))&mask; // Roll reverse k-mer

			if(x<0) {len=0; rkmer=0;} else {len++;} // Reset on ambiguous base
			if(len>=k) {
				// Apply wildcard mask and store both orientations
				int maskedKmer=(kmer&Query.midMask);
				int maskedRkmer=(rkmer&Query.midMask);
				index.put(maskedKmer, i-k+1); // Store k-mer start position
				if(maskedKmer!=maskedRkmer) { // Avoid duplicate entries for palindromes
					index.put(maskedRkmer, i-k+1);
				}
			}
		}
		return index;
	}

	/**
	 * Processes alignment hits and generates SAM output.
	 * Creates properly formatted SAM lines with CIGAR strings, mapping quality, and flags.
	 * @param q Query sequence that was aligned
	 * @param ref Reference sequence that was aligned to
	 * @param hits List of alignment start positions
	 * @param reverseStrand True if alignments are to reverse complement
	 * @param count Running count of alignments for this query (affects primary flag)
	 * @param bsw Output stream writer for SAM data
	 * @return Number of alignments processed
	 */
	static int processHits(Query q, Read ref, IntList hits, boolean reverseStrand, int count,
			ByteStreamWriter bsw) {
		if(hits==null || hits.size()==0) {return 0;}
		ByteBuilder bb=new ByteBuilder();
		ByteBuilder match=new ByteBuilder(q.bases.length);

		for(int i=0; i<hits.size(); i++) {
			count++;
			int start=hits.get(i);

			// Use appropriate query sequence for match calculation
			byte[] querySeq = reverseStrand ? q.rbases : q.bases;
			toMatch(querySeq, ref.bases, start, match.clear());

			SamLine sl=new SamLine();
			sl.pos=Math.max(start+1, 1); // Convert to 1-based SAM coordinates
			sl.qname=q.name;
			sl.setRname(ref.id);
			sl.seq=q.bases; // Always report original query sequence
			sl.qual=q.quals;
			sl.setPrimary(count==1); // Primary if first hit overall
			sl.setMapped(true);
			if(reverseStrand) {sl.setStrand(Shared.MINUS);}
			sl.tlen=q.bases.length;
			sl.cigar=SamLine.toCigar14(match.toBytes(), start, start+q.bases.length-1, ref.length(), q.bases);
			int subs=sl.countSubs();
			sl.addOptionalTag("NM:i:"+subs); // Add edit distance
			sl.mapq=Tools.mid(0, (int)(40*(sl.length()*0.5-subs)/(sl.length()*0.5)), 40); // Calculate mapping quality
			sl.toBytes(bb).nl();
			if(bb.length()>=16384) { // Batch output for efficiency
				if(bsw!=null) {bsw.addJob(bb);}
				bb=new ByteBuilder();
			}
		}
		if(bsw!=null && !bb.isEmpty()) {bsw.addJob(bb);}
		return hits.size;
	}

	/**
	 * Generates match string showing alignment quality at each position.
	 * Used for CIGAR string generation and alignment visualization.
	 * @param query Query sequence bases
	 * @param ref Reference sequence bases
	 * @param rStart Starting position in reference
	 * @param match Output builder for match string
	 */
	static void toMatch(byte[] query, byte[] ref, int rStart, ByteBuilder match) {
		for(int i=0, j=rStart; i<query.length; i++, j++) {
			boolean inbounds=(j>=0 && j<ref.length);
			byte q=query[i];
			byte r=(inbounds ? ref[j] : (byte)'$'); // Use sentinel for out-of-bounds
			boolean good=(q==r && AminoAcid.isFullyDefined(q));
			match.append(good ? 'm' : inbounds ? 'S' : 'C'); // m=match, S=substitution, C=clip
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Worker thread for processing reference sequences against query set.
	 * This class is static to prevent accidental writing to shared variables.
	 * Each thread processes a portion of the reference stream independently.
	 */
	class ProcessThread extends Thread {

		/**
		 * Constructor for ProcessThread.
		 * @param cris_ Shared input stream for reference sequences
		 * @param bsw_ Shared output stream for alignment results
		 * @param qList List of query sequences to align against
		 * @param maxSubs_ Maximum substitutions allowed per alignment
		 * @param tid_ Thread identifier
		 */
		ProcessThread(final ConcurrentReadInputStream cris_, final ByteStreamWriter bsw_, 
				ArrayList<Query> qList, final int maxSubs_, final int tid_){
			cris=cris_;
			bsw=bsw_;
			queries=qList;
			maxSubs=maxSubs_;
			tid=tid_;
		}

		/**
		 * Main thread execution method called by start().
		 * Processes reference sequences and performs alignments.
		 */
		@Override
		public void run(){
			synchronized(this) {
				processInner(); // Process all assigned reference sequences
				success=true; // Indicate successful completion
			}
		}

		/**
		 * Core processing loop that fetches and processes reference sequence batches.
		 * Continues until input stream is exhausted.
		 */
		void processInner(){
			ListNum<Read> ln=cris.nextList(); // Grab the first batch of reference sequences

			while(ln!=null && ln.size()>0){
				processList(ln); // Process this batch
				
				cris.returnList(ln); // Notify input stream that batch was processed
				ln=cris.nextList(); // Fetch next batch
			}

			// Clean up final batch
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}

		/**
		 * Processes a batch of reference sequences.
		 * @param ln ListNum containing batch of reference sequences
		 */
		void processList(ListNum<Read> ln){
			final ArrayList<Read> refList=ln.list;

			// Process each reference sequence in the batch
			for(int idx=0; idx<refList.size(); idx++){
				final Read ref=refList.get(idx);

				if(!ref.validated()){ref.validate(true);} // Validate in worker threads for speed

				// Track statistics
				final int initialLength1=ref.length();
				final int initialLength2=ref.mateLength();
				readsProcessedT+=ref.pairCount();
				basesProcessedT+=initialLength1+initialLength2;

				// Align all queries against this reference sequence
				processRefSequence(ref);
			}
		}

		/**
		 * Finds seed hits for a query against a reference index.
		 * @param q Query sequence
		 * @param refIndex K-mer index of reference sequence
		 * @param reverseStrand True to use reverse complement of query
		 * @param hitCounts Reusable map for counting hits per position
		 * @param rname Reference name (for debugging)
		 * @return List of alignment start positions with sufficient seed hits
		 */
		private IntList getSeedHits(Query q, IntListHashMap refIndex, 
				boolean reverseStrand, IntHashMap hitCounts, String rname) {
			if(useSeedMap) {return getSeedHitsMap(q, refIndex, reverseStrand, hitCounts);}
			else {return getSeedHitsList(q, refIndex, reverseStrand);}
		}

		/**
		 * Finds seed hits using list-based approach.  Potentially slow with short kmers.
		 * @param q Query sequence to search for
		 * @param refIndex K-mer index of reference
		 * @param reverseStrand True to use reverse complement query
		 * @return List of alignment positions meeting minimum hit threshold
		 */
		private IntList getSeedHitsList(Query q, IntListHashMap refIndex, boolean reverseStrand) {
			int[] queryKmers=reverseStrand ? q.rkmers : q.kmers;
			if(queryKmers==null) {return null;}
			final int minHits=Math.max(minSeedHits, q.minHits);
			if(prescan) {
				int valid=prescan(q, refIndex, reverseStrand, minHits);
				if(valid<minHits) {return null;}
			}
			IntList seedHits=new IntList();//TODO: Test speed with and without lazy init

			// Check each query k-mer for matches in reference
			for(int i=0; i<queryKmers.length; i+=qStep) {
				if(queryKmers[i]==-1) {continue;} // Skip invalid k-mers

				IntList positions=refIndex.get(queryKmers[i]);
				if(positions!=null) {
					for(int j=0; j<positions.size; j++) {
						int refPos=positions.array[j];
						int alignStart=refPos-i; // Adjust for k-mer position in query
						seedHits.add(alignStart);
					}
				}
			}
			if(seedHits==null) {return null;}

			seedHitsT+=seedHits.size();
			if(seedHits.size<minHits) {return null;}
			
			// Remove duplicates and filter by minimum occurrence
			if(seedHits.size>1 || minHits>1) {
				seedHits.sort();
				seedHits.condenseMinCopies(minHits);
			}
			alignmentsT+=seedHits.size();
			return seedHits.isEmpty() ? null : seedHits;
		}

		/**
		 * Finds seed hits using map-based approach for efficient hit counting.
		 * More efficient when many hits are expected per position or with long kmers.
		 * @param q Query sequence to analyze
		 * @param refIndex Reference index for k-mer lookup
		 * @param reverseStrand Flag indicating reverse complement strand
		 * @param hitCounts Reusable map for counting hits per alignment position
		 * @return List of alignment start positions meeting hit threshold, or null if empty
		 */
		private IntList getSeedHitsMap(Query q, IntListHashMap refIndex, 
				boolean reverseStrand, IntHashMap hitCounts) {
			final int[] queryKmers=reverseStrand ? q.rkmers : q.kmers;
			if(queryKmers==null){return null;}
			final int minHits=Math.max(minSeedHits, q.minHits);
			//TODO: Early exit if there are not enough ref kmers.
			if(prescan) {
				int valid=prescan(q, refIndex, reverseStrand, minHits);
				if(valid<minHits) {return null;}
			}

			IntList seedHits=null;

			// Process query k-mers at specified step interval
			for(int i=0; i<queryKmers.length; i+=qStep) {
				if(queryKmers[i]==-1){continue;} // Skip invalid k-mers

				IntList positions=refIndex.get(queryKmers[i]);
				if(positions!=null) {
					seedHitsT+=positions.size;
					if(seedHits==null) { // Lazy allocation
						seedHits=new IntList();
						if(hitCounts==null) {hitCounts=new IntHashMap();}
						else {hitCounts.clear();}
					}
					for(int j=0; j<positions.size; j++) {
						int alignStart=positions.array[j]-i; // Calculate alignment start

						// Increment hit count and add to results when threshold met
						int newCount=hitCounts.increment(alignStart);
						if(newCount==minHits) {seedHits.add(alignStart);}
					}
				}
			}

			alignmentsT+=(seedHits==null ? 0 : seedHits.size());
			return seedHits==null || seedHits.isEmpty() ? null : seedHits;
		}

		/**
		 * Finds seed hits using map-based approach for efficient hit counting.
		 * More efficient when many hits are expected per position.
		 * @param q Query sequence to analyze
		 * @param refIndex Reference index for k-mer lookup
		 * @param reverseStrand Flag indicating reverse complement strand
		 * @return Number of query kmers shared with the ref
		 */
		private int prescan(Query q, IntListHashMap refIndex, boolean reverseStrand, final int minHits) {
			final int[] queryKmers=reverseStrand ? q.rkmers : q.kmers;
			final int maxMisses=queryKmers.length-minHits;//TODO: Does not account for qStep
			if(queryKmers==null || maxMisses<0){return 0;}

			int misses=0;
			// Process query k-mers at specified step interval
			for(int i=0; i<queryKmers.length && misses<=maxMisses; i+=qStep) {
				if(queryKmers[i]==-1){continue;} // Skip invalid k-mers
				boolean hit=refIndex.containsKey(queryKmers[i]);
				misses+=(hit ? 0 : 1);
			}

			return queryKmers.length-misses;
		}

		/**
		 * Dispatches reference sequence processing to indexed or brute force method.
		 * @param ref Reference sequence to align queries against
		 * @return Number of alignments found
		 */
		long processRefSequence(final Read ref) {
			return indexQueries ? processRefSequenceIndexed(ref) : processRefSequenceBrute(ref);
		}

		/**
		 * Processes reference sequence using indexed seed-and-extend approach.
		 * More efficient for longer references with selective seed hits.
		 * @param ref Reference sequence to process
		 * @return Total number of alignments found across all queries
		 */
		long processRefSequenceIndexed(final Read ref){
//			Timer t=new Timer();
//			t.start("Indexing ref.");
			IntListHashMap refIndex=buildReferenceIndex(ref.bases); // Build k-mer index for this reference
			IntHashMap seedMap=(useSeedMap ? new IntHashMap() : null);
//			t.stopAndStart("Time:");
//			System.err.println("refIndex: "+refIndex.size());
//			System.err.println("seedMap: "+(seedMap==null ? 0 : seedMap.size()));
//			Shared.printMemory();

			long sum=0;
			for(Query q : queries) {
				int count=0;

				// Forward strand alignment
				IntList seedHits=getSeedHits(q, refIndex, false, seedMap, ref.name());
				IntList hits=alignSparse(q.bases, ref.bases, maxSubs, q.maxClips, seedHits);
				count+=processHits(q, ref, hits, false, 0, bsw);
//				System.err.println("seedHits+:"+(seedHits==null ? 0 : seedHits.size()));

				// Reverse strand alignment
				seedHits=getSeedHits(q, refIndex, true, seedMap, ref.name());
				hits=alignSparse(q.rbases, ref.bases, maxSubs, q.maxClips, seedHits);
				count+=processHits(q, ref, hits, true, count, bsw);
//				System.err.println("seedHits-:"+(seedHits==null ? 0 : seedHits.size()));

				readsOutT+=count;
				basesOutT+=count*q.bases.length;
				sum+=count;
			}
//			t.stop("Time:");
//			Shared.printMemory();
			return sum;
		}

		/**
		 * Processes reference sequence using brute force alignment at all positions.
		 * Used when indexing is disabled or for short references.
		 * @param ref Reference sequence to process
		 * @return Total number of alignments found across all queries
		 */
		long processRefSequenceBrute(final Read ref){
			long sum=0;
			for(Query q : queries) {
				int count=0;

				// Forward strand alignment - test all positions
				IntList hits=alignAllPositions(q.bases, ref.bases, maxSubs, q.maxClips);
				count+=processHits(q, ref, hits, false, 0, bsw);

				// Reverse strand alignment using pre-computed reverse complement query
				hits=alignAllPositions(q.rbases, ref.bases, maxSubs, q.maxClips);
				count+=processHits(q, ref, hits, true, count, bsw);

				readsOutT+=count;
				basesOutT+=count*q.bases.length;
				sum+=count;
			}
			alignmentsT+=(queries.size()*(long)ref.length()); // All positions tested
			return sum;
		}

		/** Number of reads processed by this thread */
		protected long readsProcessedT=0;
		/** Number of bases processed by this thread */
		protected long basesProcessedT=0;
		/** Number of alignments performed by this thread */
		protected long alignmentsT=0;
		/** Number of seed hits found by this thread */
		protected long seedHitsT=0;

		/** Number of reads retained by this thread */
		protected long readsOutT=0;
		/** Number of bases retained by this thread */
		protected long basesOutT=0;

		/** True only if this thread has completed successfully */
		boolean success=false;

		/** Shared input stream for reference sequences */
		private final ConcurrentReadInputStream cris;
		/** Shared output stream for alignment results */
		private final ByteStreamWriter bsw;
		/** List of query sequences to align */
		private final ArrayList<Query> queries;
		/** Maximum substitutions allowed per alignment */
		final int maxSubs;
		/** Thread identifier */
		final int tid;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file path for queries */
	private String in1=null;
	/** Secondary input file path for queries */
	private String in2=null;

	/** Primary output file path for alignments */
	private String out1=null;

	/** Override input file extension */
	private String extin=null;
	/** Override output file extension */
	private String extout=null;

	/** Reference file path for alignment targets */
	String refFile=null;
	/** Maximum substitutions allowed per alignment */
	int maxSubs=5;
	/** K-mer length for indexing */
	int k=13;
	/** Length of middle region to mask in k-mers for fuzzy matching */
	int midMaskLen=1;
	/** Whether to build k-mer indices for queries */
	boolean indexQueries=true;
	/** Count seed hits before filling list */
	boolean prescan=false;
	/** Step size for sampling query k-mers (1=all k-mers, 2=every other, etc.) */
	int qStep=1;

	/** Minimum seed hits required for alignment consideration */
	int minSeedHits=1;
	/** Probability threshold for calculating minimum hits based on query length */
	private float minHitsProb=0.9999f;
	/** Whether to use map-based seed hit counting (vs list-based) */
	boolean useSeedMap=false;

	/*--------------------------------------------------------------*/

	/** Number of reference sequences processed */
	protected long readsProcessed=0;
	/** Number of reference bases processed */
	protected long basesProcessed=0;

	/** Total number of alignment operations performed */
	protected long alignmentCount=0;
	/** Total number of seed hits found */
	protected long seedHitCount=0;

	/** Number of alignment results output */
	protected long readsOut=0;
	/** Number of alignment result bases output */
	protected long basesOut=0;

	/** Quit after processing this many input reads; -1 means no limit */
	private long maxReads=-1;

	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary input file format for queries */
	private final FileFormat ffin1;
	/** Secondary input file format for queries */
	private final FileFormat ffin2;

	/** Primary output file format for alignments */
	private final FileFormat ffout1;

	@Override
	public final ReadWriteLock rwlock() {return rwlock;}
	/** Read-write lock for thread-safe operations */
	private final ReadWriteLock rwlock=new ReentrantReadWriteLock();

	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;

}
