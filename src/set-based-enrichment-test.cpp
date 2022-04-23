#include <Rcpp.h>
#include <Rmath.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>

using namespace Rcpp;
using namespace std;

class Key {
public:
  int a, b;
  Key(int a_a, int a_b) : a(a_a), b(a_b){}
  
  bool operator < (const Key& k) const {
    if(a == k.a){
      return b < k.b;
    }
    return a < k.a;
  }
};

/*
 * Sample k elements from length(S)=n to length(R)=k. (Reservoir sampling) (Unused yet.)
 * Should replace int j = 42 with some meaningful random generator.
 * https://en.wikipedia.org/wiki/Reservoir_sampling
 */
void genRandomSample(int* S, int* R, int n, int k){
  // fill the reservoir array
  for(int i = 0; i < k; i++){
    R[i] = S[i];
  }
  
  // replace elements with gradually decreasing probability
  for(int i = k; i < n; i++){
    int j = 42;
    if(j < k){
      R[j] = S[i];
    }
  }
}

int intersectSize(unordered_set<int>& s, NumericVector& v){
  int sum = 0;
  for(int k = 0; k < v.size(); k++){
    if(s.count(v[k])>0){
      ++sum;
    }  
  }
  return sum;
}

void computeFrequencies(list<int>* reversed_DB, int* categoryFrequency, NumericVector selectedGeneIds){
  for(int i = 0; i < selectedGeneIds.size(); i++){
    int geneId = selectedGeneIds[i];
    list<int> list_of_category_id = reversed_DB[geneId];
    for(list<int>::iterator it = list_of_category_id.begin(); it != list_of_category_id.end(); ++it){
      categoryFrequency[*it]++;
    }
  }
}


void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

///////////////////////////////////////////////////////////////////////////////////////////////

/**
 * list_of_all_genes: [gene1, gene2, ... geneN] (list of strings)
 * The order of the list_of_all_genes define the ids.
 * Then, the 
 *  id of gene1 is 1, 
 *  id of gene 2 is 2, 
 *  and so on.
 * 
 */
// [[Rcpp::export]]
DataFrame enrichment_test_simulation(Rcpp::List DB_of_categories, 
                                     std::vector<std::string> list_of_all_genes, 
                                     std::vector<std::string> pool, 
                                     int selectSize, int steps, 
                                     unsigned int randomSeed){
  // PARAMETERS
  //
  // Rcpp::List DB_of_categories
  // vector<string> list_of_all_genes
  // vector<string> pool
  // int selectSize
  // int steps
  // unsigned int randomSeed   // to initialise the random generator of the current thread. It is important for multithreading
  
  
  
  // For multithreading and testing.
  // It is important to initialise the random seed for each thread (if we use multithread processing)
  
  if(randomSeed!=0){set_seed(randomSeed);}
  
  // Association: gene-id -> list(category ids)
  list<int> *reversed_DB;
  // Reverse index of the list_of_all_genes array. gene-name -> gene-id
  map<string, int> geneDict;
  
  
  NumericVector poolIds(pool.size());
  //NumericMatrix pvals(DB_of_categories.size(), steps);
  
  // Compute the gene dictionary
  for(size_t i = 0; i < list_of_all_genes.size(); i++){
    geneDict[list_of_all_genes.at(i)] = i;
  }
  
  // transform the pool to pool ids
  for(size_t i = 0; i < pool.size(); i++){
    poolIds[i] = geneDict[pool[i]];
  }
  
  // geneId -> list of category ids has the gene geneId
  // preparing the gene index "map"
  // The geneIndex assigns a list of categoriy ids for each gene
  reversed_DB = new list<int>[list_of_all_genes.size()];
  for(int catId = 0; catId < DB_of_categories.size(); catId++){
    // Gene names in the current category
    List currCategory = (List) DB_of_categories[catId];
    for(int i = 0; i < currCategory.size(); i++){
      string geneName = currCategory[i];
      int geneId = geneDict[geneName];
      reversed_DB[geneId].push_back(catId);
    }
  }
  
  
  // --------------------------------
  
  
  
  // Association: category id -> the size of the intersection of the category and the pool.
  int *categoryPoolFrequency = new int[DB_of_categories.size()]();
  computeFrequencies(reversed_DB, categoryPoolFrequency, poolIds);
  
  
  // Association: (intersection of pool, intersection of select) -> count
  map<Key, int> intersectStat;
  
  // Do the experiments.
  for(int step = 0; step < steps; step++){
    NumericVector randomSelectIds = sample(poolIds, selectSize);
    
    // Association: category id -> the size of the intersection of the category and the select.
    int *categorySelectFrequency = new int[DB_of_categories.size()](); // this variable will be the return value
    computeFrequencies(reversed_DB, categorySelectFrequency, randomSelectIds);
    
    // 
    for(int catId = 0; catId < DB_of_categories.size(); catId++){
      int intersectRandomSelectSize = categorySelectFrequency[catId]; // intersect(select , DB_category)
      int intersectPoolSize = categoryPoolFrequency[catId];           // intersect(pool , DB_category)
      
      Key k(intersectPoolSize, intersectRandomSelectSize);
      
      if(intersectStat.count(k) == 0){
        intersectStat.insert(pair<Key,int>(k, 1)); 
      }else{
        intersectStat.find(k)->second++;
      }
    }
    delete[] categorySelectFrequency;
  }
  
  // console output - I do not want it 
  //  cout << "Done." << endl;
  
  // create data frame from the results
  
  IntegerVector poolIntersectSize;
  IntegerVector selectIntersectSize;
  IntegerVector count;
  
  std::map<Key,int>::iterator it = intersectStat.begin();
  for (it=intersectStat.begin(); it!=intersectStat.end(); ++it){
    poolIntersectSize.push_back(it->first.a);
    selectIntersectSize.push_back(it->first.b);
    count.push_back(it->second);
    // Rcout << "(" << it->first.a << "," << it->first.b << "," << it->second << ")\n";
  }
  
  delete[] reversed_DB;
  delete[] categoryPoolFrequency;
  
  return DataFrame::create(_["pool.intersect"] = poolIntersectSize, _["select.intersect"] = selectIntersectSize, _["count"] = count);
}

///////////////////////////////////////////////////////////////////////////////////////////////

