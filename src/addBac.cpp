// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::depends(RcppEigen)]]

//#include <RcppEigen.h> // sparse matrix
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List addBacCpp(arma::sp_mat occmat, DataFrame orgdat, int amount, double growth, int type, int ptype) {
//Eigen::SparseMatrix<int> addBacCpp(Eigen::SparseMatrix<int> occmat, DataFrame orgdat, int amount, int type) {
//void addBacCpp(Eigen::SparseMatrix<int> occmat, DataFrame orgdat, int amount, int type) {

  //cout << occmat;
  int n = occmat.n_cols;
  int m = occmat.n_rows;
  //int n=100;
  //int m=100;
  int lastind = orgdat.nrows();
  
  //access columns of dataframe
  DoubleVector  growth_vec = orgdat["growth"];
  IntegerVector type_vec = orgdat["type"];
  IntegerVector ptype_vec = orgdat["phenotype"];
  IntegerVector x_vec = orgdat["x"];
  IntegerVector y_vec = orgdat["y"];
  for(int i=1; i<=amount; i++){
    growth_vec.push_back(growth);
    ptype_vec.push_back(ptype);
    type_vec.push_back(type);
    int x = static_cast<int>(round(as<double>(runif(1, 0, n-1))));
    int y = static_cast<int>(round(as<double>(runif(1, 0, m-1))));
    //cout << x, y;
    while(occmat(x,y) != 0){
      x = static_cast<int>(round(as<double>(runif(1, 0, n-1))));
      y = static_cast<int>(round(as<double>(runif(1, 0, m-1))));
    }
    occmat(x,y) = type;
    x_vec.push_back(x);
    y_vec.push_back(y);
  }
  DataFrame new_orgdat = DataFrame::create(_["growth"]= growth_vec, _["type"]= type_vec, _["ptype"]= ptype_vec, _["x"]= x_vec, _["y"]= y_vec);
  //cout << occmat;
  List l; l["occmat"]=occmat; l["orgdat"]=new_orgdat;
  return(l);
}


