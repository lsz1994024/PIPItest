package proteomics.Types;

public class ResultEntry implements Comparable<ResultEntry> {

    public final double score;
    public final String peptide;
    private final int hashCode;
    private final boolean isDecoy;
    private final Peptide0 peptide0;

    public ResultEntry(Peptide0 peptide0, double score, String peptide, boolean isDecoy) {
        this.score = score;
        this.peptide = peptide;
        this.peptide0 = peptide0;
        String toString = peptide + "-" + score;
        hashCode = toString.hashCode();
        this.isDecoy = isDecoy;
    }

    public int compareTo(ResultEntry other) {
        return Double.compare(score, other.score);
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof ResultEntry) {
            ResultEntry temp = (ResultEntry) other;
            return ((temp.peptide.contentEquals(peptide)) && (temp.score == score));
        } else {
            return false;
        }
    }

    public Peptide0 getPeptide0() {
        return this.peptide0;
    }
}