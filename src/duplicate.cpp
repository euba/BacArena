#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
DataFrame duplicateCpp(DataFrame orgdat, int n, int m, List cellweight, IntegerMatrix occupyM){
  DataFrame orgdat_c = clone(orgdat);
  
  IntegerVector x_vec      = orgdat_c["x"];
  IntegerVector y_vec      = orgdat_c["y"];
  DoubleVector  growth_vec = orgdat_c["growth"];
  IntegerVector type_vec   = orgdat_c["type"];
  IntegerVector pheno_vec  = orgdat_c["phenotype"];
  
  int size = x_vec.size();

  //create occupy matrix
  IntegerMatrix occ(n,m);
  int z = 0;
  while(z < size){
    occ(x_vec[z]-1, y_vec[z]-1) = 1;
    z++;
  }

  z = 0;
  while(z < size){
    int sel = static_cast<int>(round(as<double>(runif(1, 0, size-1))));

    //check if growth is high enough to duplicate
    if(growth_vec[sel]>cellweight[type_vec[sel]-1]){
      int old_x = x_vec[sel];
      int old_y = y_vec[sel];      
      std::vector<std::pair<int,int> > list_free; 
      for(int i=-1; i<=1; i++){
        for(int j=-1; j<=1; j++){
          int pos_x = i + old_x;
          int pos_y = j + old_y;
          if(pos_x > 0 && pos_y > 0 && pos_x <= n && pos_y <= m){
            if(occ(pos_x-1, pos_y-1) == 0 && (pos_x != old_x || pos_y != old_y)){
              std::pair <int,int> free_pos(pos_x,pos_y);
              list_free.push_back(free_pos);
            }
          }    
        }
      }
      if(list_free.size() > 0){
        int rnd_element = static_cast<int>(round(Rcpp::as<double>(Rcpp::runif(1, 0, list_free.size()-1))));
        int new_x = list_free[rnd_element].first;
        int new_y = list_free[rnd_element].second;
        if(occupyM(new_x-1, new_y-1)==0){
          //std::cout<< "duplicating: " << growth_vec[sel] << old_x << "," << old_y << " to " << new_x << "," << new_y << "\n";
          occ(new_x-1, new_y-1) = 1;
          double new_growth = growth_vec[sel]/2.0;
          growth_vec[sel] = new_growth;
          growth_vec.push_back(new_growth);
          x_vec.push_back(new_x);
          y_vec.push_back(new_y);
          type_vec.push_back(type_vec[sel]);
          pheno_vec.push_back(pheno_vec[sel]);
        }
        else{
          //std::cout<<"do not duplicating: " << growth_vec[sel] << old_x << "," << old_y << " to " << new_x << "," << new_y << " because of " << occupyM(new_x-1, new_y-1) << "\n";
        }
      }
    }
    z++;
  }
  return(DataFrame::create(_["growth"]= growth_vec, _["type"]= type_vec, _["phenotype"]=pheno_vec, _["x"]=x_vec, _["y"]=y_vec));
}