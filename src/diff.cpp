#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void diffuseGrajdeanuCpp(Rcpp::NumericMatrix y, bool donut, double mu){
  // According to Grajdeanu (2007) MODELING DIFFUSION IN A DISCRETE ENVIRONMENT
  int n = y.ncol();
  int m = y.nrow();
  //int z = 0;
  // double mu = 1.0; // diff const
  Rcpp::NumericMatrix y_old = Rcpp::clone(y);
  
  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      double A = 0, d = 0, sum_A = 0, sum_y = 0;
      int neigh = 0;
      for(int l=-1; l<=1; l++){
        for(int o=-1; o<=1; o++){
          int pos_i; //neighbour with position i,j
          int pos_j;
          if (donut==true){
            pos_i = (i + l) % m;
            pos_j = (j + o) % n;
            if (pos_i == -1) pos_i = m-1;
            if (pos_j == -1) pos_j = n-1;
          }
          else{ // in case of bounds (no donut) ignore unavaible neighbours
            pos_i = i + l;
            pos_j = j + o;
            if (pos_i == -1 || pos_i == m) pos_i = i;
            if (pos_j == -1 || pos_j == n) pos_j = j;
          }
          if(i!=pos_i or j!=pos_j){ //don't get it again
            neigh++;
            d = 1.0; // cell distance
            sum_A += exp( -pow(d,2) / mu );
            sum_y += ( y_old(pos_i,pos_j) - y_old(i,j) ) * exp( -pow(d,2)/mu );
          }
        }
      }
      A = 1.0 / sum_A;
      y(i,j) = y_old(i,j) + A * sum_y;
    }
  }
}


// [[Rcpp::export]]
void diffuseNaiveCpp(Rcpp::NumericMatrix y, bool donut){
  int n = y.ncol();
  int m = y.nrow();
  int z = 0;
  while(z++ < n*m ){ // one repeat for each cell
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