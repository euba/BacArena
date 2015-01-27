#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
void diffuseNaiveCpp(Rcpp::NumericMatrix y, bool donut){
  int n = y.ncol();
  int m = y.nrow();
  int l = 0;
  while(l++ < n*m ){ // one repeat for each cell
    // okay this is ugly but who cares ;P  (needed for using same seed as in R)
    int j = static_cast<int>(round(Rcpp::as<double>(Rcpp::runif(1, 0, n-1))));
    int i = static_cast<int>(round(Rcpp::as<double>(Rcpp::runif(1, 0, m-1))));
    double min = y(i,j);
    std::vector<std::pair<int,int> > list_min; //list containing minima
    for(int l=-1; l<=1; l++){
      for(int o=-1; o<=1; o++){
        int pos_i;
        int pos_j;
        if (donut==true){
          pos_i = (i + l) % m;
          pos_j = (j + o) % n;
          //if (pos_i>=0 && pos_i<n && pos_j>=0 && pos_j<m){ //boundary check
          if (pos_i == -1) pos_i = m-1;
          if (pos_j == -1) pos_j = n-1;
        }
        else{ // in case of bounds (no donut) ignore unavaible neighbours
          pos_i = i + l;
          pos_j = j + o;
          if (pos_i == -1 || pos_i == m) pos_i = i;
          if (pos_j == -1 || pos_j == n) pos_j = j;
        }
        if (y(pos_i,pos_j) < min) {
          //std::cout << y(pos_i,pos_j);
          min = y(pos_i,pos_j);
          list_min.clear(); // clear list with old minima
          std::pair <int,int> new_min(pos_i,pos_j);
          list_min.push_back(new_min);
        }
        else if (y(pos_i,pos_j) == min) {
          std::pair <int,int> new_min(pos_i,pos_j);
          list_min.push_back(new_min);
        }
       }
     }
     // choose randomly one of the found minima
     int rnd_element = static_cast<int>(round(Rcpp::as<double>(Rcpp::runif(1, 0, list_min.size()-1))));
     int min_i = list_min[rnd_element].first;
     int min_j = list_min[rnd_element].second;
     double mean = (y(i,j) + y(min_i,min_j)) / 2.0;
     y(i,j) = mean;
     y(min_i, min_j) = mean;
  }
}