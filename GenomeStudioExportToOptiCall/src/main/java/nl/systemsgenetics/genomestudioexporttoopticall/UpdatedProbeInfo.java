package nl.systemsgenetics.genomestudioexporttoopticall;

/**
 *
 * @author Patrick Deelen
 */
public class UpdatedProbeInfo {
	
	private final String id;
	private final String chr;
	private final int pos;
	private final String refAllele;
	private final String altAlelle;
	private final char aAllele;
	private final char bAllele;

	public UpdatedProbeInfo(String id, String chr, String pos, String refAllele, String altAlelle, String aAllele, String bAllele) {
		this.id = id;
		this.chr = chr;
		this.pos = Integer.parseInt(pos);
		this.refAllele = refAllele;
		this.altAlelle = altAlelle;
		
		if(aAllele.length() > 1 || bAllele.length() > 1){
			throw new RuntimeException();
		}
		
		this.aAllele = aAllele.charAt(0);
		this.bAllele = bAllele.charAt(0);
	}

	public String getId() {
		return id;
	}

	public String getChr() {
		return chr;
	}

	public int getPos() {
		return pos;
	}

	public String getRefAllele() {
		return refAllele;
	}

	public String getAltAlelle() {
		return altAlelle;
	}

	public char getaAllele() {
		return aAllele;
	}

	public char getbAllele() {
		return bAllele;
	}
	
	
	
}