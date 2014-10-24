package nl.systemsgenetics.qtlsnpbedannotator;

import java.util.List;
import umcg.genetica.collections.intervaltree.PerChrIntervalTree;
import umcg.genetica.io.bed.BedEntry;

/**
 *
 * @author Patrick Deelen
 */
public class BedAnnotationSource {
	
	private final String name;
	private final PerChrIntervalTree<BedEntry> intervalTree;

	public BedAnnotationSource(String name, PerChrIntervalTree<BedEntry> intervalTree) {
		this.name = name;
		this.intervalTree = intervalTree;
	}
	
	public List<BedEntry> searchPosition(String chr, int pos) {
		return intervalTree.searchPosition(chr, pos);
	}

	public String getName() {
		return name;
	}
	
}
