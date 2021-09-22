package proteomics;

import ProteomicsLibrary.Score;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import ProteomicsLibrary.Binomial;
import proteomics.Search.CalSubscores;
import proteomics.Search.Search;
import proteomics.Segment.InferSegment;
import ProteomicsLibrary.PrepareSpectrum;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import proteomics.Spectrum.PreSpectra;
import proteomics.Types.*;
import sun.reflect.generics.tree.Tree;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.sql.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

public class PIPIWrap implements Callable<PIPIWrap.Entry> {

    private final BuildIndex buildIndex;
    private final MassTool massTool;
    private final double ms1Tolerance;
    private final double leftInverseMs1Tolerance;
    private final double rightInverseMs1Tolerance;
    private final int ms1ToleranceUnit;
    private final double ms2Tolerance;
    private final double minPtmMass;
    private final double maxPtmMass;
    private final int localMaxMs2Charge;
    private final Map<String, Peptide0> peptide0Map;
    private final JMzReader spectraParser;
    private final double minClear;
    private final double maxClear;
    private final ReentrantLock lock;
    private final String scanId;
    private final int precursorCharge;
    private final double precursorMass;
    private final InferPTM inferPTM;
    private final PrepareSpectrum preSpectrum;
    private final String sqlPath;
    private final Binomial binomial;


    public PIPIWrap(BuildIndex buildIndex, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double ms2Tolerance, double minPtmMass, double maxPtmMass, int localMaxMs2Charge, JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanId, int precursorCharge, double precursorMass, InferPTM inferPTM, PrepareSpectrum preSpectrum, String sqlPath, Binomial binomial) {
        this.buildIndex = buildIndex;
        this.massTool = massTool;
        this.ms1Tolerance = ms1Tolerance;
        this.leftInverseMs1Tolerance = leftInverseMs1Tolerance;
        this.rightInverseMs1Tolerance = rightInverseMs1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.localMaxMs2Charge = localMaxMs2Charge;
        this.spectraParser = spectraParser;
        this.minClear = minClear;
        this.maxClear = maxClear;
        this.lock = lock;
        this.scanId = scanId;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.inferPTM = inferPTM;
        this.preSpectrum = preSpectrum;
        this.sqlPath = sqlPath;
        this.binomial = binomial;
        peptide0Map = buildIndex.getPeptide0Map();
    }

    @Override
    public Entry call() throws Exception {
        Map<Double, Double> rawPLMap;
        try {
            lock.lock();
            // Reading peak list.
            rawPLMap = spectraParser.getSpectrumById(scanId).getPeakList();
        } finally {
            lock.unlock();
        }

        // preprocess peak list, remove noise and normalize the peaks then find the topN peaks in each 100 mz range
//        System.out.println("raw PlMap");
//        for (Map.Entry<Double, Double> pl :rawPLMap.entrySet() )
//        {
//            System.out.println(pl.getKey()+","+pl.getValue());
//        }
//        System.out.println("min"+minClear+"max"+maxClear);
        TreeMap<Double, Double> plMap = preSpectrum.preSpectrumTopNStyle(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, PreSpectra.topN);  // remove reporter ions and keep the local top peaks
        TreeMap<Double, Double> NoiseIntensity = preSpectrum.NoiseCalculate(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, PreSpectra.topN);
//        TreeMap<Double, Double> plMap = new TreeMap<>(rawPLMap);
        if (plMap.isEmpty()) {
            return null;
        }

        if (NoiseIntensity.isEmpty()) {
            return null;
        }

        // Coding the spectrum
        InferSegment inferSegment = buildIndex.getInferSegment();
        // infer the tag and add virtual peak, expAaLists is the list of the inferred tag list of spectrum
        List<ThreeExpAA> expAaLists = inferSegment.inferSegmentLocationFromSpectrum(precursorMass, plMap, NoiseIntensity);
//        for (ThreeExpAA tag3 : expAaLists)
//        {
//
//            System.out.print(tag3.printTag()+ "  ");
//        }
        List<TwoExpAA> expAaLists2 = inferSegment.inferSegmentLocationFromSpectrum2(precursorMass, plMap, NoiseIntensity);
        List<FourExpAA> expAaLists4 = inferSegment.inferSegmentLocationFromSpectrum4(precursorMass, plMap, NoiseIntensity);
        if (!expAaLists.isEmpty()) {
            // get the spectrum-coded vector, if there are two same tags in expAaList, the scanCode only uses the one with the higher intensity value
            // scanCode is the intensity and the corresponding vector bit position
            SparseBooleanVector scanCode = inferSegment.generateSegmentIntensityVector(expAaLists);  // the expAaList regards PTT and TTP as different things, scanCode regards them as the same.
            SparseBooleanVector scanCode2 = inferSegment.generateSegmentIntensityVector2(expAaLists2);
            SparseBooleanVector scanCode4 = inferSegment.generateSegmentIntensityVector4(expAaLists4);
            // Begin search.
//            if(scanId.equals("4771")){
//                System.out.println(scanId);
//            }
            Connection sqlC = DriverManager.getConnection(sqlPath);
            Statement sqlS = sqlC.createStatement();
            ResultSet sqlR = sqlS.executeQuery(String.format(Locale.US, "SELECT scanNum FROM spectraTable WHERE scanId='%s'", scanId));
            int scanN = 0;
            if (sqlR.next()) {
                scanN = sqlR.getInt("scanNum");
                sqlR.close();
                sqlS.close();
            }

            Search search = new Search(buildIndex, precursorMass, scanCode, scanCode2, scanCode4, NoiseIntensity, scanN, massTool, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, localMaxMs2Charge);

            // prepare the spectrum
            SparseBooleanVector expProcessedPL;
            if (PIPI.useXcorr) {
                expProcessedPL = preSpectrum.prepareXCorr(plMap, false);
            } else {
                expProcessedPL = preSpectrum.digitizePL(plMap);
            }

            // dynamic programming
            double localMS1ToleranceL = -1 * ms1Tolerance;
            double localMS1ToleranceR = ms1Tolerance;
            if (ms1ToleranceUnit == 1) {
                localMS1ToleranceL = (precursorMass * leftInverseMs1Tolerance) - precursorMass;
                localMS1ToleranceR = (precursorMass * rightInverseMs1Tolerance) - precursorMass;
            }

            // infer PTM using the new approach
            TreeSet<Peptide> peptideSet = new TreeSet<>(Collections.reverseOrder());
            Map<String, TreeSet<Peptide>> modSequences = new TreeMap<>();
            for (Peptide peptide : search.getPTMOnlyResult()) {
                Peptide0 peptide0 = peptide0Map.get(peptide.getPTMFreePeptide());
                PeptidePTMPattern peptidePTMPattern = inferPTM.tryPTM(expProcessedPL, plMap, precursorMass, peptide.getPTMFreePeptide(), peptide.isDecoy(), peptide.getNormalizedCrossCorr(), peptide0.leftFlank, peptide0.rightFlank, peptide.getGlobalRank(), precursorCharge, localMaxMs2Charge, localMS1ToleranceL, localMS1ToleranceR);
                if (!peptidePTMPattern.getPeptideTreeSet().isEmpty()) {
                    for (Peptide tempPeptide : peptidePTMPattern.getPeptideTreeSet()) {
                        if (tempPeptide.getScore() > 0) {
                            if (peptideSet.size() < 100) {
                                peptideSet.add(tempPeptide);
                            } else if (tempPeptide.getScore() > peptideSet.last().getScore()) {
                                peptideSet.pollLast();
                                peptideSet.add(tempPeptide);
                            }
                        }
                    }
                    // record scores with different PTM patterns for calculating PTM delta score.
                    modSequences.put(peptidePTMPattern.ptmFreePeptide, peptidePTMPattern.getPeptideTreeSet());
                }
            }

            // Calculate Score for PTM free peptide
            for (Peptide peptide : search.getPTMFreeResult()) {
                double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrix(), precursorCharge, expProcessedPL);
                if (score > 0) {
                    peptide.setScore(score);
                    peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMs2Charge, peptide.getIonMatrix(), ms2Tolerance));
                    if (peptideSet.size() < 100) {
                        peptideSet.add(peptide);
                    } else if (peptide.getScore() > peptideSet.last().getScore()) {
                        peptideSet.pollLast();
                        peptideSet.add(peptide);
                    }
                }
            }

            if (!peptideSet.isEmpty()) {
                Peptide[] peptideArray = peptideSet.toArray(new Peptide[0]);
                Peptide topPeptide = peptideArray[0];
                TreeSet<Peptide> ptmPatterns = null;
                if (topPeptide.hasVarPTM()) {
                    ptmPatterns = modSequences.get(topPeptide.getPTMFreePeptide());
                }
                new CalSubscores(topPeptide, ms2Tolerance, plMap, precursorCharge, ptmPatterns, binomial);

                Connection sqlConnection = DriverManager.getConnection(sqlPath);
                Statement sqlStatement = sqlConnection.createStatement();
                ResultSet sqlResultSet = sqlStatement.executeQuery(String.format(Locale.US, "SELECT scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, isDecoy, score FROM spectraTable WHERE scanId='%s'", scanId));
                if (sqlResultSet.next()) {
                    int scanNum = sqlResultSet.getInt("scanNum");
                    String scanId = sqlResultSet.getString("scanId");
                    int precursorCharge = sqlResultSet.getInt("precursorCharge");
                    double precursorMass = sqlResultSet.getDouble("precursorMass");
                    String mgfTitle = sqlResultSet.getString("mgfTitle");
                    int isotopeCorrectionNum = sqlResultSet.getInt("isotopeCorrectionNum");
                    double ms1PearsonCorrelationCoefficient = sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient");

                    double deltaLCn = 1;
                    if (peptideArray.length > 4) {
                        deltaLCn = (peptideArray[0].getScore() - peptideArray[4].getScore()) / peptideArray[0].getScore();
                    }
                    double deltaCn = 1;
                    if (peptideArray.length > 1) {
                        deltaCn = (peptideArray[0].getScore() - peptideArray[1].getScore()) / peptideArray[0].getScore();
                    }

                    String otherPtmPatterns = "-";
                    if (ptmPatterns != null) {
                        List<String> tempList = new LinkedList<>();
                        Iterator<Peptide> ptmPatternsIterator = ptmPatterns.iterator();
                        ptmPatternsIterator.next();
                        while (ptmPatternsIterator.hasNext()) {
                            Peptide temp = ptmPatternsIterator.next();
                            tempList.add(String.format(Locale.US, "%s-%.4f", temp.getPtmContainingSeq(buildIndex.returnFixModMap()), temp.getScore())); // Using 4 decimal here because it is write the the result file for checking. It is not used in scoring or other purpose.
                        }
                        otherPtmPatterns = String.join(";", tempList);
                    }

                    String PTMFreeCandidates = "";
                    if (search.getPTMFreeResult() != null){
                        List<String> tempList = new LinkedList<>();
                        for (Peptide temp : search.getPTMFreeResult()) {
                            tempList.add(String.format(Locale.US, "%s", temp.getPTMFreePeptide()));
                        }
                        PTMFreeCandidates = String.join(";", tempList);
                    }

                    String PTMOnlyCandidates = "";
                    if (search.getPTMOnlyResult() != null){
                        List<String> tempList = new LinkedList<>();
                        for (Peptide temp : search.getPTMOnlyResult()) {
                            tempList.add(String.format(Locale.US, "%s", temp.getPTMFreePeptide()));
                        }
                        PTMOnlyCandidates = String.join(";", tempList);
                    }

                    Entry entry = new Entry(scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, buildIndex.getLabelling(), topPeptide.getPtmContainingSeq(buildIndex.returnFixModMap()), topPeptide.getTheoMass(), topPeptide.isDecoy() ? 1 : 0, topPeptide.getGlobalRank(), topPeptide.getNormalizedCrossCorr(), topPeptide.getScore(), deltaLCn, deltaCn, topPeptide.getMatchedPeakNum(), topPeptide.getIonFrac(), topPeptide.getMatchedHighestIntensityFrac(), topPeptide.getExplainedAaFrac(), otherPtmPatterns, topPeptide.getaScore(), PTMFreeCandidates, PTMOnlyCandidates);

                    sqlResultSet.close();
                    sqlStatement.close();
                    return entry;
                } else {
                    throw new NullPointerException(String.format(Locale.US, "There is no record %s in the spectraTable.", scanId));
                }
            } else {
                return null;
            }
        } else {
            return null;
        }
    }


    public class Entry {

        final int scanNum;
        final String scanId;
        final int precursorCharge;
        final double precursorMass;
        final String mgfTitle;
        final int isotopeCorrectionNum;
        final double ms1PearsonCorrelationCoefficient;
        final String labelling;
        public final String peptide;
        final double theoMass;
        final int isDecoy;
        final int globalRank;
        final double normalizedCorrelationCoefficient;
        public final double score;
        final double deltaLCn;
        final double deltaCn;
        final int matchedPeakNum;
        final double ionFrac;
        final double matchedHighestIntensityFrac;
        final double explainedAaFrac;
        final String otherPtmPatterns; // It has 4 decimal because it is write the the result file for checking. It is not used in scoring or other purpose.
        final String aScore;
        final String PTMFreeCandidates;
        final String PTMOnlyCandidates;

        Entry(int scanNum, String scanId, int precursorCharge, double precursorMass, String mgfTitle, int isotopeCorrectionNum, double ms1PearsonCorrelationCoefficient, String labelling, String peptide, double theoMass, int isDecoy, int globalRank, double normalizedCorrelationCoefficient, double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac, double explainedAaFrac, String otherPtmPatterns, String aScore, String PTMFreeCandidates, String PTMOnlyCandidates) {
            this.scanNum = scanNum;
            this.scanId = scanId;
            this.precursorCharge = precursorCharge;
            this.precursorMass = precursorMass;
            this.mgfTitle = mgfTitle;
            this.isotopeCorrectionNum = isotopeCorrectionNum;
            this.ms1PearsonCorrelationCoefficient = ms1PearsonCorrelationCoefficient;
            this.labelling = labelling;
            this.peptide = peptide;
            this.theoMass = theoMass;
            this.isDecoy = isDecoy;
            this.globalRank = globalRank;
            this.normalizedCorrelationCoefficient = normalizedCorrelationCoefficient;
            this.score = score;
            this.deltaLCn = deltaLCn;
            this.deltaCn = deltaCn;
            this.matchedPeakNum = matchedPeakNum;
            this.ionFrac = ionFrac;
            this.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
            this.explainedAaFrac = explainedAaFrac;
            this.otherPtmPatterns = otherPtmPatterns;
            this.aScore = aScore;
            this.PTMFreeCandidates = PTMFreeCandidates;
            this.PTMOnlyCandidates = PTMOnlyCandidates;
        }
    }
}
