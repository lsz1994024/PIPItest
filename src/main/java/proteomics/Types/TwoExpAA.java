package proteomics.Types;

public class TwoExpAA implements Comparable<TwoExpAA> {

    private final ExpAA[] twoExpAa;
    private int hashCode;
    private final double totalIntensity;
    private final String ptmFreeAAString;
    private int regionIdx;

    public TwoExpAA(ExpAA aa1, ExpAA aa2, Double noise) {
        twoExpAa = new ExpAA[]{aa1, aa2};
        String toString = twoExpAa[0].toString() + "-" + twoExpAa[1].toString();
        hashCode = toString.hashCode();

        StringBuilder sb = new StringBuilder(5);
        for (ExpAA aa : twoExpAa) {
            sb.append(aa.getPtmFreeAA());
        }
        ptmFreeAAString = sb.toString();

        double intensity = twoExpAa[0].getHeadIntensity();
        for (ExpAA aa : twoExpAa) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity * intensity / noise;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        return (other instanceof TwoExpAA) && (this.hashCode() == other.hashCode());
    }

    public boolean approximateEquals(TwoExpAA other, double tolerance) {
        for (int i = 0; i < this.size(); ++i) {
            if (!this.get(i).approximateEquals(other.get(i), tolerance)) {
                return false;
            }
        }
        return true;
    }

//    public void setTheoLocation(int i, int theoLoc) {
//        twoExpAa[i].setTheoLocation(theoLoc);
//        // update toString and hashCode
//        String toString = twoExpAa[0].toString() + "-" + twoExpAa[1].toString() + "-" + twoExpAa[2].toString();
//        hashCode = toString.hashCode();
//    }

    public int compareTo(TwoExpAA other) {
        return Double.compare(twoExpAa[0].getHeadLocation(), other.twoExpAa[0].getHeadLocation());
    }

//    public ExpAA[] getExpAAs() {
//        return twoExpAa;
//    }

    public String getPtmFreeAAString() {
        return ptmFreeAAString;
    }

    public double getTotalIntensity() {
        return totalIntensity;
    }

    public double getHeadLocation() {
        return twoExpAa[0].getHeadLocation();
    }

    public double getTailLocation() {
        return twoExpAa[twoExpAa.length - 1].getTailLocation();
    }

    public TwoExpAA clone() throws CloneNotSupportedException {
        super.clone();
        return new TwoExpAA(twoExpAa[0].clone(), twoExpAa[1].clone(), 1.0);
    }

    public int size() {
        return twoExpAa.length;
    }

    public ExpAA get(int i) {
        return twoExpAa[i];
    }

    public void setRegionIdx(int regionIdx) {
        this.regionIdx = regionIdx;
    }

    public int getRegionIdx() {
        return regionIdx;
    }
}
