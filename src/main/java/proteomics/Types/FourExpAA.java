package proteomics.Types;

public class FourExpAA implements Comparable<FourExpAA> {

    private final ExpAA[] fourExpAa;
    private int hashCode;
    private final double totalIntensity;
    private final String ptmFreeAAString;
    private int regionIdx;

    public FourExpAA(ExpAA aa1, ExpAA aa2, ExpAA aa3, ExpAA aa4, Double noise) {
        fourExpAa = new ExpAA[]{aa1, aa2, aa3, aa4};
        String toString = fourExpAa[0].toString() + "-" + fourExpAa[1].toString() + "-" + fourExpAa[2].toString() + "-" + fourExpAa[3].toString();
        hashCode = toString.hashCode();

        StringBuilder sb = new StringBuilder(5);
        for (ExpAA aa : fourExpAa) {
            sb.append(aa.getPtmFreeAA());
        }
        ptmFreeAAString = sb.toString();

        double intensity = fourExpAa[0].getHeadIntensity();
        for (ExpAA aa : fourExpAa) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity * intensity / noise;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        return (other instanceof FourExpAA) && (this.hashCode() == other.hashCode());
    }

    public boolean approximateEquals(FourExpAA other, double tolerance) {
        for (int i = 0; i < this.size(); ++i) {
            if (!this.get(i).approximateEquals(other.get(i), tolerance)) {
                return false;
            }
        }
        return true;
    }

    public void setTheoLocation(int i, int theoLoc) {
        fourExpAa[i].setTheoLocation(theoLoc);
        // update toString and hashCode
        String toString = fourExpAa[0].toString() + "-" + fourExpAa[1].toString() + "-" + fourExpAa[2].toString();
        hashCode = toString.hashCode();
    }

    public int compareTo(FourExpAA other) {
        return Double.compare(fourExpAa[0].getHeadLocation(), other.fourExpAa[0].getHeadLocation());
    }

    public ExpAA[] getExpAAs() {
        return fourExpAa;
    }

    public String getPtmFreeAAString() {
        return ptmFreeAAString;
    }

    public double getTotalIntensity() {
        return totalIntensity;
    }

    public double getHeadLocation() {
        return fourExpAa[0].getHeadLocation();
    }

    public double getTailLocation() {
        return fourExpAa[fourExpAa.length - 1].getTailLocation();
    }

    public FourExpAA clone() throws CloneNotSupportedException {
        super.clone();
        return new FourExpAA(fourExpAa[0].clone(), fourExpAa[1].clone(), fourExpAa[2].clone(), fourExpAa[3].clone(), 1.0);
    }

    public int size() {
        return fourExpAa.length;
    }

    public ExpAA get(int i) {
        return fourExpAa[i];
    }

    public void setRegionIdx(int regionIdx) {
        this.regionIdx = regionIdx;
    }

    public int getRegionIdx() {
        return regionIdx;
    }
}
