package proteomics.Search;

import org.omg.Messaging.SYNC_WITH_TRANSPORT;
import proteomics.Index.BuildIndex;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import proteomics.Types.*;
import proteomics.Search.KMeans.*;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;


import java.io.IOException;
import java.util.*;
import java.util.Collections;
import java.util.stream.Collector;
import java.util.EnumSet;
import java.util.Objects;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;


public class Search {

    private static final int rankNum = 20;

    private List<Peptide> ptmOnlyResult = new LinkedList<>();
    private List<Peptide> ptmFreeResult = new LinkedList<>();


    public Search(BuildIndex buildIndex, double precursorMass, SparseBooleanVector scanCode, SparseBooleanVector scanCode2, SparseBooleanVector scanCode4, TreeMap<Double, Double> NoiseIntensity, int scanN, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double minPtmMass, double maxPtmMass, int localMaxMs2Charge) throws IOException {
//        String groundtruth = null;
//        try {
//            BufferedReader reader = new BufferedReader(new FileReader("02A.pin.selected_psm.csv"));
//            reader.readLine();
//            String line = null;
//            while((line=reader.readLine())!=null){
//                String item[] = line.split(",");
//
//                String last = item[0];
//                int value = Integer.parseInt(last);
//                if (value == scanN){
//                    groundtruth = item[1];
//                }
//            }
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//        System.out.println(scanCode.getKeys().size());
        PriorityQueue<ResultEntry> ptmFreeQueue = new PriorityQueue<>(rankNum * 2);
        PriorityQueue<ResultEntry> ptmOnlyQueue = new PriorityQueue<>(rankNum * 2);
        ArrayList<Peptide0> tempPTMFreeList = new ArrayList<>();
        ArrayList<Peptide0> tempPTMOnlyList = new ArrayList<>();
        double scanNormSquare = scanCode.norm2square();
        double scanNormSquare2 = scanCode2.norm2square();
        double scanNormSquare4 = scanCode4.norm2square();
        // sort the tags by intensity and SNR
        Map<Integer, Double> resultsMap = sortScanCode(scanCode);
        Map<Integer, Double> resultsMap2 = sortScanCode(scanCode2);
        Map<Integer, Double> resultsMap4 = sortScanCode(scanCode4);

        double leftTol = ms1Tolerance;
        double rightTol = ms1Tolerance;
        if (ms1ToleranceUnit == 1) {
            leftTol = precursorMass - (precursorMass * leftInverseMs1Tolerance);
            rightTol = (precursorMass * rightInverseMs1Tolerance) - precursorMass;
        }
        double leftMass = Math.max(precursorMass + minPtmMass - leftTol, buildIndex.getMinPeptideMass());
        double rightMass = Math.min(precursorMass + maxPtmMass + rightTol, buildIndex.getMaxPeptideMass());

        if (leftMass >= rightMass) {
            return;
        }

        Map<String, Peptide0> peptide0Map = buildIndex.getPeptide0Map();
        TreeMap<Double, Set<String>> massPeptideMap = buildIndex.getMassPeptideMap();

        NavigableMap<Double, Set<String>> subMassPeptideMap = massPeptideMap.subMap(leftMass, true, rightMass, true);

//        Boolean flllaaaggg = false;
        if (!subMassPeptideMap.isEmpty()) {
            for (double mass : subMassPeptideMap.keySet()) {
                for (String sequence : massPeptideMap.get(mass)) {
                    Peptide0 peptide0 = peptide0Map.get(sequence);

                    List codeKeys = new ArrayList(peptide0.code.getKeys());
                    List codeKeys2 = new ArrayList(peptide0.code2.getKeys());
                    List codeKeys4 = new ArrayList(peptide0.code4.getKeys());
                    List resultsKeys = new ArrayList(resultsMap.keySet());
                    List resultsKeys2 = new ArrayList(resultsMap2.keySet());
                    List resultsKeys4 = new ArrayList(resultsMap4.keySet());
                    codeKeys.retainAll(resultsKeys);
                    codeKeys2.retainAll(resultsKeys2);
                    codeKeys4.retainAll(resultsKeys4);
//                    codeKeys.retainAll(firstN);
                    double overallScore = 0;
//                    if(codeKeys.size() >= 3){
//                        System.out.println(codeKeys.size());
//                    }
                    double score = 0.0;
                    double nominator;
                    double temp1 = Math.sqrt(peptide0.code.norm2square() * scanNormSquare);
                    if (temp1 > 1e-6) {
                        nominator = peptide0.code.dot(scanCode);
                        score = nominator / temp1 * codeKeys.size() / (sequence.length() - 2);
                    }

                    double score2 = 0.0;
                    double nominator2;
                    double temp2 = Math.sqrt(peptide0.code2.norm2square() * scanNormSquare2);
                    if (temp2 > 1e-6) {
                        nominator2 = peptide0.code2.dot(scanCode2);
                        score2 = nominator2 / temp2 * codeKeys2.size() / (sequence.length() - 2);
                    }

                    if(codeKeys4.size() > 0){
                        double score4 = 0.0;
                        double nominator4;
                        double temp4 = Math.sqrt(peptide0.code4.norm2square() * scanNormSquare4);
                        if (temp4 > 1e-6) {
                            nominator4 = peptide0.code4.dot(scanCode4);
                            score4 = nominator4 / temp4 * codeKeys4.size() / (sequence.length() - 2);
                        }
                        overallScore = (score * score2 * score4) / ((score * score2 * score4) + (1 - score) * (1 - score2) * (1 - score4));
                    }
                    else{
                        overallScore = (score * score2) / ((score * score2) + (1 - score) * (1 - score2));
                    }

//                    if (codeKeys4.size() == 0 && overallScore != 0){
//                        System.out.println(overallScore);
//                    }

                    // score is cross correlation coefficient
                    double deltaMass = mass - precursorMass; // caution: the order matters under ms1ToleranceUnit == 1 situation

                    if (peptide0.isTarget) {
                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                            // in the mass tolerant range, PTM-free
                            if (ptmFreeQueue.size() < rankNum || codeKeys4.size() > 0) {
                                ptmFreeQueue.add(new ResultEntry(peptide0, overallScore, sequence, false));
                            }
                            else {
//                                if ((codeKeys.size() >= 2 && codeKeys2.size() >= 3) || (codeKeys4.size() >= 1 && codeKeys2.size() >= 4)) {
//                                    ptmFreeQueue.add(new ResultEntry(peptide0, overallScore, sequence, false));
//                                }
                                if (overallScore > ptmFreeQueue.peek().score) {
//                                    if (ptmFreeQueue.peek().score < 0.05){
                                        ptmFreeQueue.poll();
                                        ptmFreeQueue.add(new ResultEntry(peptide0, overallScore, sequence, false));
//                                    }
//                                    else{
//                                        ptmFreeQueue.add(new ResultEntry(peptide0, score, sequence, false));
//                                    }
                                }
                            }
                        }

                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                            // out of the mass tolerant range, PTM-only
                            if (ptmOnlyQueue.size() < rankNum || codeKeys4.size() > 0) {
                                ptmOnlyQueue.add(new ResultEntry(peptide0, overallScore, sequence, false));
                            }
                            else {
//                                if ((codeKeys.size() >= 2 && codeKeys2.size() >= 3) || (codeKeys4.size() >= 1 && codeKeys2.size() >= 4)) {
//                                    ptmFreeQueue.add(new ResultEntry(peptide0, overallScore, sequence, false));
//                                }
                                if (overallScore > ptmOnlyQueue.peek().score) {
//                                    if (ptmOnlyQueue.peek().score < 0.15){
                                        ptmOnlyQueue.poll();
                                        ptmOnlyQueue.add(new ResultEntry(peptide0, overallScore, sequence, false));
//                                    }
//                                    else{
//                                        ptmOnlyQueue.add(new ResultEntry(peptide0, score, sequence, false));
//                                    }
                                }
                            }
                        }
                    } else {
                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                            // PTM-free
                            if (ptmFreeQueue.size() < rankNum || codeKeys4.size() > 0) {
                                ptmFreeQueue.add(new ResultEntry(peptide0, overallScore, sequence, true));
                            }
                            else {
//                                if ((codeKeys.size() >= 2 && codeKeys2.size() >= 3) || (codeKeys4.size() >= 1 && codeKeys2.size() >= 4)) {
//                                    ptmFreeQueue.add(new ResultEntry(peptide0, overallScore, sequence, true));
//                                }
                                if (overallScore > ptmFreeQueue.peek().score) {
//                                    if (ptmFreeQueue.peek().score < 0.05){
                                        ptmFreeQueue.poll();
                                        ptmFreeQueue.add(new ResultEntry(peptide0, overallScore, sequence, true));
//                                    }
//                                    else{
//                                        ptmFreeQueue.add(new ResultEntry(peptide0, score, sequence, true));
//                                    }
                                }
                            }
                        }

                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                            // PTM-only
                            if (ptmOnlyQueue.size() < rankNum || codeKeys4.size() > 0) {
                                ptmOnlyQueue.add(new ResultEntry(peptide0, overallScore, sequence, true));
                            }
                            else {
//                                if ((codeKeys.size() >= 2 && codeKeys2.size() >= 2) || (codeKeys4.size() >= 1 && codeKeys2.size() >= 4)) {
//                                    ptmFreeQueue.add(new ResultEntry(peptide0, overallScore, sequence, true));
//                                }
                                if (overallScore > ptmOnlyQueue.peek().score) {
//                                    if (ptmOnlyQueue.peek().score < 0.15){
                                        ptmOnlyQueue.poll();
                                        ptmOnlyQueue.add(new ResultEntry(peptide0, overallScore, sequence, true));
//                                    }
//                                    else{
//                                        ptmOnlyQueue.add(new ResultEntry(peptide0, score, sequence, true));
//                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

//        if (ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty()) {
//            System.out.println("empty");
//        }

//        PriorityQueue<ResultEntry> ptmFreeQueue_copy = new PriorityQueue<>(ptmFreeQueue);
//        PriorityQueue<ResultEntry> ptmOnlyQueue_copy = new PriorityQueue<>(ptmOnlyQueue);
//        // Queue to ArrayList
//        if (!(ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty())) {
//            tempPTMFreeList = convert(ptmFreeQueue);
//            tempPTMOnlyList = convert(ptmOnlyQueue);
//        }
//
//        KMeans kmeans = new KMeans();
//        if (tempPTMFreeList.size() > (subMassPeptideMap.size() * 0.02)){         // if the candidate list is too long, use kmeans to make it short
//            int init_k = (int) Math.floor(tempPTMFreeList.size() / (rankNum * 4));
////            System.out.println("Free");
//            tempPTMFreeList = kmeans.KMeans(tempPTMFreeList, scanCode, rankNum * 4, init_k);
//        }
//
//        if (tempPTMOnlyList.size() > (subMassPeptideMap.size() * 0.02)){
//            int init_k = (int) Math.floor(tempPTMOnlyList.size() / (rankNum * 4));
////            System.out.println("Only");
//            tempPTMOnlyList = kmeans.KMeans(tempPTMOnlyList, scanCode, rankNum * 4, init_k);
//        }
//
//        // ArrayList to Queue
//        PriorityQueue<ResultEntry> tempFreeQueue = new PriorityQueue<>(rankNum * 2);
//        PriorityQueue<ResultEntry> tempOnlyQueue = new PriorityQueue<>(rankNum * 2);
//        if (!(tempPTMFreeList.isEmpty() && tempPTMOnlyList.isEmpty())) {
//            tempFreeQueue = convertback(tempPTMFreeList, ptmFreeQueue_copy);
//            tempOnlyQueue = convertback(tempPTMOnlyList, ptmOnlyQueue_copy);
//        }

//        if (!(tempFreeQueue.isEmpty() && tempOnlyQueue.isEmpty())) {
////            System.out.println(1);
//            ptmFreeResult = convertResult(tempFreeQueue, massTool, localMaxMs2Charge);
//            ptmOnlyResult = convertResult(tempOnlyQueue, massTool, localMaxMs2Charge);
//        }

//        if (!flllaaaggg){
//            System.out.println("NOOOOOOO");
//        }

        if (!(ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty())) {
            ptmFreeResult = convertResult(ptmFreeQueue, massTool, localMaxMs2Charge);
            ptmOnlyResult = convertResult(ptmOnlyQueue, massTool, localMaxMs2Charge);
        }
    }

    private PriorityQueue<ResultEntry> convertback(ArrayList<Peptide0> inputList, PriorityQueue<ResultEntry> inputQueue) {
        PriorityQueue<ResultEntry> resultQueue = new PriorityQueue<>(rankNum * 2);
        while (!inputQueue.isEmpty()) {
            ResultEntry temp = inputQueue.poll();
            if (inputList.contains(temp.getPeptide0())){
                resultQueue.add(temp);
            }
        }
        return resultQueue;
    }

    private Map<Integer, Double> sortScanCode(SparseBooleanVector scanCode){
        // get the vector
        Map<Integer, Double> scanCodeList = scanCode.sparseVector;
        // sort the vector
        ValueComparator bvc =  new ValueComparator(scanCodeList);
        TreeMap<Integer, Double> sortedMap = new TreeMap<Integer, Double>(bvc);
        sortedMap.putAll(scanCodeList);
        // get the first 10
        List<Map.Entry<Integer,Double>> firstN = sortedMap.entrySet().stream().limit(20).collect(Collectors.toList());
        // convert list to map
        Map<Integer, Double> resultsMap = new HashMap<>();
        for (int i = 0; i < firstN.size(); i ++) {
            resultsMap.put(firstN.get(i).getKey(), firstN.get(i).getValue());
        }
        return resultsMap;
    }

    private ArrayList<Peptide0> convert(PriorityQueue<ResultEntry> inputQueue) {
        ArrayList<Peptide0> result = new ArrayList<>(10);
        while (!inputQueue.isEmpty()) {
            ResultEntry temp = inputQueue.poll();
            result.add(temp.getPeptide0());
        }
        return result;
    }

    private List<Peptide> convertResult(PriorityQueue<ResultEntry> inputQueue, MassTool massTool, int localMaxMs2Charge) {
        List<Peptide> peptideList = new LinkedList<>();
        int globalRank = inputQueue.size();
        while (!inputQueue.isEmpty()) {
            ResultEntry temp = inputQueue.poll();
            peptideList.add(new Peptide(temp.peptide, temp.isDecoy(), massTool, localMaxMs2Charge, temp.score, globalRank));
            -- globalRank;
        }
        return peptideList;
    }

    public List<Peptide> getPTMOnlyResult() {
        return ptmOnlyResult;
    }

    public List<Peptide> getPTMFreeResult() {
        return ptmFreeResult;
    }
}

class ValueComparator implements Comparator<Integer> {

    Map<Integer, Double> base;
    public ValueComparator(Map<Integer, Double> base) {
        this.base = base;
    }

    // Note: this comparator imposes orderings that are inconsistent with equals.
    public int compare(Integer a, Integer b) {
        if (base.get(a) >= base.get(b)) {
            return -1;
        } else {
            return 1;
        } // returning 0 would merge keys
    }
}