# Diffusion with correct math model
src_math_diffusion <- '
Rcpp::List xlist(A);
int s = xlist.size();
int seed2 = as<int>(seed);  
srand (seed2); /* initialize random seed: */

for(int k=0; k<s; k++) {
  SEXP ll = xlist[k];
  Rcpp::NumericMatrix y(ll);
  
  int n = y.ncol();
  int m = y.nrow();
  A = 1



}
'




# Diffusion by mean with lowest neighbour
# randomly choosen cells
#
src_diffusion <- '
  Rcpp::List xlist(A);
  Rcpp::NumericMatrix y(A);
  
  int n = y.ncol();
  int m = y.nrow();
  std::cout << n;
  std::cout << m;

  int l = 0;
  while(l++ < n*m ){ // repeat as long as cells exist

    // okay this is ugly but who cares ;P  (needed for using same seed as in R)
    int i = static_cast<int>(round(Rcpp::as<double>(Rcpp::runif(1, 0, n-1))));
    int j = static_cast<int>(round(Rcpp::as<double>(Rcpp::runif(1, 0, m-1))));
    std::cout << i << "," << j << "  ";
    double min = y(i,j);
    std::vector<std::pair<int,int> > list_min; //list containing minima
    for(int l=-1; l<=1; l++){
      for(int o=-1; o<=1; o++){
        int pos_i = (i + l) % n;
        int pos_j = (j + o) % m;
        //if (pos_i>=0 && pos_i<n && pos_j>=0 && pos_j<m){ //boundary check
        if (pos_i == -1) pos_i = n-1;
        if (pos_j == -1) pos_j = m-1;

        if (y(pos_i,pos_j) < min) {
          std::cout << y(pos_i,pos_j);
          min = y(pos_i,pos_j);
          list_min.clear(); // clear list with old minima
          std::pair <int,int> new_min(pos_i,pos_j);
          list_min.push_back(new_min);
        }
        else if (y(pos_i,pos_j) == min) {
          std::pair <int,int> new_min(pos_i,pos_j);
          list_min.push_back(new_min);
        }
        //}
      }
    }
    // choose randomly one of the found minima
    int rnd_element = static_cast<int>(round(Rcpp::as<double>(Rcpp::runif(1, 0, list_min.size()-1))));
    int min_i = list_min[rnd_element].first;
    int min_j = list_min[rnd_element].second;
    
    double m = (y(i,j) + y(min_i,min_j)) / 2.0;

    y(i,j) = m;
    y(min_i, min_j) = m;
  }
'



# Diffusion by mean with lowest neighbour
# synchron updating!
#
src_diffusion3 <- '
  Rcpp::List xlist(A);
  int n = xlist.size();
  for(int k=0; k<n; k++) {
    SEXP ll = xlist[k];
    Rcpp::NumericMatrix y(ll);
    //copy matrix to diffuse successive
    //Rcpp::NumericMatrix tmp = Rcpp::clone(y);
    //std::cout << y(1,1);
    int n = y.ncol();
    int m = y.nrow();
    for(int i=0; i<n; i++){
      for(int j=0; j<m; j++){
        double min = y(i,j);
        std::vector<std::pair<int,int> > list_min; //list containing minima
        for(int l=-1; l<=1; l++){
          for(int o=-1; o<=1; o++){
            int pos_i = i + l;
            int pos_j = j + o;
            if(pos_i>=0 && pos_i<n && pos_j>=0 && pos_j<m){ //boundary check
              if(y(pos_i,pos_j) < min) {
                min = y(pos_i,pos_j);
                list_min.clear(); // clear list with old minima
                std::pair <int,int> new_min(pos_i,pos_j);
                list_min.push_back(new_min);
              }
              else if(y(pos_i,pos_j) == min) {
                std::pair <int,int> new_min(pos_i,pos_j);
                list_min.push_back(new_min);
              }
            }
          }
        }
        std::random_shuffle ( list_min.begin(), list_min.end() ); //randomize order
        // choose randomly one of the found minima
        
        int min_i = list_min.front().first;
        int min_j = list_min.front().second;
        double m = (y(i,j) + y(min_i,min_j)) / 2.0;

        y(i,j) = m;
        y(min_i, min_j) = m;
      }
    }
  }
  //return(A);
'



# Diffusion by moving average
src_diffusion2 <- '
  Rcpp::List xlist(A);
  int n = xlist.size();
  for(int k=0; k<n; k++) {
    SEXP ll = xlist[k];
    Rcpp::NumericMatrix y(ll);
    //copy matrix to diffuse sucsessive
    Rcpp::NumericMatrix tmp = Rcpp::clone(y);
    //std::cout << y(1,1);
    int n = tmp.ncol();
    int m = tmp.nrow();
    for(int i=0; i<n; i++){
      for(int j=0; j<m; j++){
        //get neighbours
        int sum = 0;
        int nb = 0;
        for(int l=-1; l<=1; l++){
          for(int o=-1; o<=1; o++){
            int pos_i = i + l;
            int pos_j = j + o;
            if(pos_i>=0 && pos_i<n && pos_j>=0 && pos_j<m){
              nb++;
              sum += tmp(pos_i, pos_j);
            }
          }
        }
        y(i,j) = (double)sum/(double)nb;
      }
    }
  }
  //return(A);
'

src_movement <- '
  Rcpp::DataFrame l(input_frame);
  Rcpp::IntegerVector x = l["x"];
  Rcpp::IntegerVector y = l["y"];
  //int seed2 = as<int>(seed);
  int n = as<int>(length_n);
  int m = as<int>(length_m);
  
  x=x-1;  // index in c++ starts at 0
  y=y-1;
  
  /* initialize random seed: */
  //srand (seed2);
  
  //int n = 10;
  //int m = 10;
  
  //std::cout << l.length();
  //std::cout << as<std::vector>(l(0))[1] << std::endl;
  
  std::set<std::pair<int,int> > s;
  for(int i=0; i < x.size(); i++){
  s.insert(std::make_pair(x(i), y(i)));
  }

  for(int i=0; i < x.size(); i++){
    int a = (x(i) + rand() % 3 - 1);
    int b = (y(i) + rand() % 3 - 1);
    if(a == -1) a = n-1;
    if(b == -1) b = m-1;
    if(a == n) a = 0;
    if(b == m) b = 0;
 
    if(s.find(std::make_pair(a,b)) == s.end()){ // if empty go for it!
      s.insert(std::make_pair(a,b));
      s.erase(std::make_pair(x(i),y(i)));
      //std::cout << "move: (" << x(i)-1 << "," << y(i)-1 << ") -> (" << a+1 << "," << b+1 << ")" << std::endl;
      x(i) = a; 
      y(i) = b;
    }
  }  
  
  x=x+1; // index in R starts at 1
  y=y+1;
  Rcpp::DataFrame new_l = Rcpp::DataFrame::create(Rcpp::Named("x")=x, Rcpp::Named("y")=y, Rcpp::Named("type")=l[2]);
  return(Rcpp::List::create(Rcpp::Named("df")=new_l));
  '