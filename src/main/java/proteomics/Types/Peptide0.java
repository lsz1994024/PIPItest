package proteomics.Types;


import ProteomicsLibrary.Types.SparseBooleanVector;

public class Peptide0 {

    public final SparseBooleanVector code;
    public final SparseBooleanVector code2;
    public final SparseBooleanVector code4;
    public final boolean isTarget;
    public final String[] proteins;
    public final char leftFlank;
    public final char rightFlank;

    public Peptide0(SparseBooleanVector code, SparseBooleanVector code2, SparseBooleanVector code4, boolean isTarget, String[] proteins, char leftFlank, char rightFlank) {
        this.code = code;
        this.code2 = code2;
        this.code4 = code4;
        this.isTarget = isTarget;
        this.proteins = proteins;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
    }
}
