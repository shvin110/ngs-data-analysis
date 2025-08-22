package aligner;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Implements a WaveFront alignment algorithm for global alignment.
 * This version is actually fast.
 * 
 *@author Brian Bushnell
 *@contributor Isla
 *@date May 27, 2025
 */
public class WaveFrontAligner2 implements IDAligner {

	/** Main() passes the args and class to Test to avoid redundant code */
	public static <C extends IDAligner> void main(String[] args) throws Exception {
	    StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
		@SuppressWarnings("unchecked")
		Class<C> c=(Class<C>)Class.forName(stackTrace[(stackTrace.length<3 ? 1 : 2)].getClassName());
		Test.testAndPrint(c, args);
	}

    public WaveFrontAligner2() {}

	@Override
	public final String name() {return "WaveFront";}
    @Override
    public final float align(byte[] a, byte[] b) {return alignStatic(a, b, null);}
    @Override
    public final float align(byte[] a, byte[] b, int[] pos) {return alignStatic(a, b, pos);}
    @Override
    public final float align(byte[] a, byte[] b, int[] pos, int minScore) {return alignStatic(a, b, pos);}
    @Override
    public final float align(byte[] a, byte[] b, int[] pos, int rStart, int rStop) {
        return alignStatic(a, b, pos, rStart, rStop);
    }
    
	private static AtomicLong loops=new AtomicLong(0);
	public long loops() {return loops.get();}
	public void setLoops(long x) {loops.set(x);}
	public static String output=null;
	
	public static float alignStatic(byte[] query, byte[] ref, int[] posVector){
	    final int qLen=query.length;
	    final int rLen=ref.length;
	    
	    // Handle special cases
	    if(qLen==0 || rLen==0){
	        if(posVector!=null){
	            posVector[0]=0;
	            posVector[1]=Math.max(0, rLen-1);
	        }
	        return qLen==0 && rLen==0?1.0f:0.0f;
	    }
	    
	    // The maximum reasonable edit distance for alignment
	    final int maxEditDist=(int)Math.max(qLen, rLen)+8;
	    
	    // Number of diagonals (k ranges from -qLen to rLen)
	    final int numDiagonals=qLen+rLen+1;
	    
	    // Diagonal offset for indexing (k=0 at index qLen)
	    final int diagOffset=qLen;
	    
	    // We'll use a rolling buffer approach with arrays tracking both position and matches
	    int[] currentWF=new int[numDiagonals];
	    int[] nextWF=new int[numDiagonals];
	    int[] currentMatches=new int[numDiagonals];  // Track matches for each wavefront position
	    int[] nextMatches=new int[numDiagonals];
	    
	    // Initialize all positions to -1 (not reached)
	    Arrays.fill(currentWF, -1);
	    Arrays.fill(nextWF, -1);
	    
	    // Start with edit distance 0, only diagonal 0 is active
	    currentWF[diagOffset]=0;
	    currentMatches[diagOffset]=0;
	    
	    // Track the result
	    int finalEditDist=-1;
	    int finalMatches=0;
	    
	    // Counter for wavefront size estimation (total cells explored)
	    int loopCounter=0, loopCounter2=0;
	    
	    // Main loop - iterate through edit distances
	    for(int d=0; d<=maxEditDist; d++){
	        boolean anyActive=false;
	        
	        // Process current wavefront
	        for(int k=-qLen; k<=rLen; k++){
	            int diagIdx=k+diagOffset;
	            
	            // Skip if this diagonal hasn't been reached yet
	            if(diagIdx<0 || diagIdx>=numDiagonals) continue;
	            
	            // Current furthest reach on this diagonal
	            int reach=currentWF[diagIdx];
	            if(reach<0) continue;
	            
	            loopCounter++;
	            
	            // Convert to matrix coordinates
	            int i=reach;
	            int j=reach+k;
	            int matchCount=currentMatches[diagIdx];
	            
	            // Extend matches as far as possible
	            while(i<qLen && j<rLen && query[i]==ref[j]){
	                i++;
	                j++;
	                matchCount++;
	                loopCounter2++;
	            }
	            
	            // Update furthest reach and match count for current position
	            currentWF[diagIdx]=i;
	            currentMatches[diagIdx]=matchCount;
	            
	            // Check if we've reached the end
	            if(i>=qLen && j>=rLen){
	                finalEditDist=d;
	                finalMatches=matchCount;
	                break;
	            }
	            
	            // Apply edit operations for the next wavefront
	            
	            // Insertion: extend to diagonal k-1
	            if(diagIdx>0){
	                int newReach = i + 1;
	                if(nextWF[diagIdx-1] < newReach || 
	                   (nextWF[diagIdx-1] == newReach && nextMatches[diagIdx-1] < matchCount)){
	                    nextWF[diagIdx-1] = newReach;
	                    nextMatches[diagIdx-1] = matchCount;
	                }
	                anyActive=true;
	            }
	            
	            // Substitution: extend to diagonal k
	            int newReach = i + 1;
	            if(nextWF[diagIdx] < newReach || 
	               (nextWF[diagIdx] == newReach && nextMatches[diagIdx] < matchCount)){
	                nextWF[diagIdx] = newReach;
	                nextMatches[diagIdx] = matchCount;
	            }
	            anyActive=true;
	            
	            // Deletion: extend to diagonal k+1
	            if(diagIdx<numDiagonals-1){
	                if(nextWF[diagIdx+1] < i || 
	                   (nextWF[diagIdx+1] == i && nextMatches[diagIdx+1] < matchCount)){
	                    nextWF[diagIdx+1] = i;
	                    nextMatches[diagIdx+1] = matchCount;
	                }
	                anyActive=true;
	            }
	        }
	        
	        // If we've found a complete alignment, we're done
	        if(finalEditDist>=0){
	            break;
	        }
	        
	        // If no active diagonals for next wavefront, we can't find an alignment
	        if(!anyActive){
	            // Didn't find complete path - use best progress
	            for(int k=-qLen; k<=rLen; k++){
	                int diagIdx=k+diagOffset;
	                if(diagIdx>=0 && diagIdx<numDiagonals && currentWF[diagIdx]>=0){
	                    int i=currentWF[diagIdx];
	                    int j=i+k;
	                    if(i==qLen || j==rLen){
	                        finalEditDist=d;
	                        finalMatches=currentMatches[diagIdx];
	                        break;
	                    }
	                }
	            }
	            if(finalEditDist<0){
	                finalEditDist=maxEditDist;
	                finalMatches=0;
	            }
	            break;
	        }
	        
	        // Swap buffers for next iteration
	        int[] temp=currentWF;
	        currentWF=nextWF;
	        nextWF=temp;
	        
	        temp=currentMatches;
	        currentMatches=nextMatches;
	        nextMatches=temp;
	        
	        // Reset the next wavefront buffer
	        Arrays.fill(nextWF, -1);
	        // Don't need to reset nextMatches - they'll be overwritten
	    }
	    
	    // If we hit max edit distance without finding alignment
	    if(finalEditDist<0){
	        finalEditDist=maxEditDist;
	        finalMatches=0;
	    }
	    
	    // Set loop count for metrics
	    loops.addAndGet(Math.max(loopCounter, loopCounter2));
	    
	    // Calculate exact identity
	    // alignment_length = matches + edits
	    // identity = matches / alignment_length
	    int alignmentLength = finalMatches + finalEditDist;
	    float identity = alignmentLength > 0 ? (float)finalMatches / alignmentLength : 0.0f;
	    
	    // Fill position vector
	    if(posVector!=null){
	        posVector[0]=0;
	        posVector[1]=Math.max(0, rLen-1);
	    }
	    
	    return identity;
	}
    
    /**
     * Wrapper for aligning within a reference window.
     */
    public static final float alignStatic(final byte[] query, final byte[] ref, 
            final int[] posVector, int refStart, int refEnd) {
        refStart = Math.max(refStart, 0);
        refEnd = Math.min(refEnd, ref.length-1);
        final int rlen = refEnd - refStart + 1;
        final byte[] region = (rlen == ref.length ? ref : Arrays.copyOfRange(ref, refStart, refEnd+1));
        final float id = alignStatic(query, region, posVector);
        if(posVector != null) {
            posVector[0] += refStart;
            posVector[1] += refStart;
        }
        return id;
    }
}