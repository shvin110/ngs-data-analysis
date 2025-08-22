package bin;

import java.util.ArrayList;

import shared.Tools;

/**
 * Abstract base class for bin refinement strategies.
 * Accepts a potentially impure bin and returns null for no change,
 * or a list of refined bins if beneficial splits are found.
 * 
 * @author Brian Bushnell & UMP45
 * @date June 20, 2025
 */
abstract class AbstractRefiner extends BinObject {

	/**
	 * Attempts to refine/split the given bin.
	 * @param input Potentially impure cluster to analyze
	 * @return null if no refinement recommended, or ArrayList of 2+ bins if split beneficial
	 */
	abstract ArrayList<Bin> refine(Bin input);

	/**
	 * Validates that a proposed split actually improves cluster quality.
	 * Override this for custom validation logic.
	 */
	protected boolean isSplitBeneficial(Bin original, ArrayList<Bin> splits) {
		if (splits == null || splits.size() < 2) return false;

		// Basic sanity checks
		long totalSize = 0;
		for (Bin bin : splits) {
			if (bin.numContigs() == 0) return false;
			totalSize += bin.size();
		}

		// Conservation of mass
		if (totalSize != original.size()) return false;

		// Each split should be reasonably sized
		for (Bin bin : splits) {
			if (bin.size() < original.size() * 0.1f) return false; // No tiny fragments
		}

		return true;
	}

	public static AbstractRefiner makeRefiner(Oracle oracle) {
		return makeRefiner(oracle, DEFAULT_TYPE);
	}

	public static AbstractRefiner makeRefiner(Oracle oracle, int type) {
		if(type==CRYSTAL) {return new CrystalChamber(oracle);}
		// Future methods for other refiner types:
		// public static AbstractRefiner createGraphCutter(Oracle oracle) { ... }
		// public static AbstractRefiner createEvidenceClusterer(Oracle oracle) { ... }
		throw new RuntimeException();
	}
	
	public static int findType(String s) {
		int idx=Tools.find(s, types);
		assert(idx>=0) : "Can't find type "+s;
		return idx;
	}

	public static final int CRYSTAL=0, GRAPH=1, EVIDENCE=2;
	public static final String types[]= {"CRYSTAL", "GRAPH", "EVIDENCE"};
	public static int DEFAULT_TYPE=CRYSTAL;
	
}
