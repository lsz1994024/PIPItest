package proteomics.Search;

import ProteomicsLibrary.Types.SparseBooleanVector;
//import ProteomicsLibrary.Types.SparseVector;
import proteomics.Types.Peptide0;
import weka.core.stopwords.Null;

import java.lang.reflect.Array;
import java.util.*;


public class KMeans {

    static int MAX_ITER = 100;    // Maximum number of iterations
    private ArrayList<Peptide0> targetCluster;
    private ArrayList<LinkedList<Peptide0>> centroidPeptideCluster;
    private int K;

    public ArrayList<Peptide0> KMeans(ArrayList<Peptide0> peptide0List, SparseBooleanVector scanCode, int numOfResults, int init_k) {

        ArrayList<SparseBooleanVector> Centroids = find_random_Centroids(peptide0List, init_k);  // init find random centroids

        boolean flag;  // if the centroids are updated
        int iter_num = 0;
        do{
            // calculate the distance between every peptide0 and every centroids, distribute peptides to clusters
            centroidPeptideCluster = locate_peptide_to_centroids(Centroids, peptide0List);  // this is OK, except if distance between peptide and all of centroid is 0, locate to cluster 0.
            flag = true;
            // update the centroids
            ArrayList<SparseBooleanVector> new_Centroids = update_centroids(centroidPeptideCluster); // this is OK.
            // if the centroids are not changed, exit the loop, else continue
            for(int i = 0; i < centroidPeptideCluster.size(); i++){
                if(!(new_Centroids.get(i).compare(Centroids.get(i)))){
                    break;
                }
                else{
                    flag = false;
                }
            }
            Centroids = new_Centroids;
            ++ iter_num;
        }while(flag && iter_num < MAX_ITER);

        double final_dist = 0.0;
        double scancCodeNormSquare = scanCode.norm2square();
        int temp_i = 0;
        for(int i = 0; i < init_k; i++){
            double temp_dist = dist(scanCode, scancCodeNormSquare, Centroids.get(i));
            if (temp_dist > final_dist){
                final_dist = temp_dist;
                temp_i = i;
            }
            targetCluster = new ArrayList<>(centroidPeptideCluster.get(temp_i));
        }

        return targetCluster;
    }

    private ArrayList<SparseBooleanVector> update_centroids(ArrayList<LinkedList<Peptide0>> centroidPeptideCluster) {
        ArrayList<SparseBooleanVector> newCentroids = new ArrayList<SparseBooleanVector>();
        double value = 0.0;

        for (int i = 0; i < centroidPeptideCluster.size(); i++) {
            SparseBooleanVector temp_vector = new SparseBooleanVector();
            // within each cluster, find new centroid by calculating mean value
            for (Peptide0 peptide : centroidPeptideCluster.get(i)) {
                for(Integer key : peptide.code.getKeys()){
                    value = peptide.code.get(key) / centroidPeptideCluster.get(i).size();
                    if(value > 1e-6){
                        temp_vector.add(key, value);
                    }
                }
            }
            newCentroids.add(temp_vector);
        }
        return newCentroids;
    }

    private ArrayList<LinkedList<Peptide0>> locate_peptide_to_centroids(ArrayList<SparseBooleanVector> Centroids, ArrayList<Peptide0> peptide0List) {
        int init_k = Centroids.size();
//        double[][] dist = new double[init_k][peptide0List.size()];
        int temp_i = 0;
        ArrayList<LinkedList<Peptide0>> array =  new ArrayList<LinkedList<Peptide0>>();
        for(int z = 0; z < init_k; z ++){
            array.add(new LinkedList<Peptide0>());
        }
        boolean no_cluster;
        for(int j = 0; j < peptide0List.size(); j++){
            no_cluster = false;
            double max_Dist = 0.0;
            double pepNormSquare = peptide0List.get(j).code.norm2square();
            for(int i = 0; i < init_k; i++){
                double tempDist = dist(peptide0List.get(j).code,  pepNormSquare, Centroids.get(i));
                if (tempDist > max_Dist){
                    max_Dist = tempDist;
                    temp_i = i;
                    no_cluster = true;
                }
            }
            if (no_cluster){
                array.get(temp_i).add(peptide0List.get(j));
            }
        }
        return array;
    }
    
    private ArrayList<SparseBooleanVector> find_random_Centroids(ArrayList<Peptide0> peptide0List, int init_k){
        ArrayList<SparseBooleanVector> Centroids = new ArrayList<SparseBooleanVector>();
        Random r = new Random();
        ArrayList<Integer> a = new ArrayList<>();
        for(int i = 0; i < init_k ; i++){
            int temp;
            do{
                temp = r.nextInt(peptide0List.size());
            }while(a.contains(temp));

            a.add(temp);
            Centroids.add(peptide0List.get(temp).code);
        }
        return Centroids;
    }

    private double dist(SparseBooleanVector centroidCode, Double centroidNormSquare, SparseBooleanVector peptide_code) {
        double score = 0;
        double temp1 = Math.sqrt(peptide_code.norm2square() * centroidNormSquare);
        if (temp1 > 1e-6) {
            score = peptide_code.dot(centroidCode) / temp1;
        }
        return score;
    }
}
