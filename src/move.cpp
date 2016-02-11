#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
void movementCpp(DataFrame orgdat, int n, int m){
  //ordat data.frame will be changed! (pointer)
  IntegerVector x_vec = orgdat["x"];
  IntegerVector y_vec = orgdat["y"];

  int size = x_vec.size();
  
  //create occupy matrix
  IntegerMatrix occ(n,m);
  int z = 0;
  while(z < size){
    occ(x_vec[z]-1, y_vec[z]-1) = 1;
    z++;
  }
  //std::cout<<occ;
  
  z = 0;
  while(z < size){
    int sel = static_cast<int>(round(as<double>(runif(1, 0, size-1))));
    int old_x = x_vec[sel];
    int old_y = y_vec[sel];
    
    std::vector<std::pair<int,int> > list_free; 
    for(int i=-1; i<=1; i++){
      for(int j=-1; j<=1; j++){
        int pos_x = i + old_x;
        int pos_y = j + old_y;
        if(pos_x > 0 && pos_y > 0 && pos_x <= n && pos_y <= m){
          if(occ(pos_x-1, pos_y-1) == 0 && pos_x != old_x && pos_y != old_y){
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
      occ(old_x-1, old_y-1) = 0;
      occ(new_x-1, new_y-1) = 1;
      x_vec[sel] = new_x;
      y_vec[sel] = new_y;
      //std::cout<<old_x<<","<<old_y<<"-->"<<new_x<<","<<new_y<<"\n";
    }
    z++;
  }
}