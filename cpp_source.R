# Diffusion by mean with lowest neighbour
# randomly choosen cells
#
src_diffusion <- '
Rcpp::List xlist(A);
int n = xlist.size();
srand (time(NULL)); /* initialize random seed: */
  for(int k=0; k<n; k++) {
    SEXP ll = xlist[k];
    Rcpp::NumericMatrix y(ll);
    //copy matrix to diffuse successive
    //Rcpp::NumericMatrix tmp = Rcpp::clone(y);
    //std::cout << y(1,1);
    int n = y.ncol();
    int m = y.nrow();
    int l = 0;
    while(l++ < n*m ){ // repeat as long as cells exist
      int i = rand() % n;
      int j = rand() % m;

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
  const Rcpp::NumericMatrix  source(A);
  Rcpp::NumericMatrix tmp = Rcpp::clone(source);

   /* initialize random seed: */
  srand (time(NULL));
  /* generate secret number between 1 and 10: */
  //iSecret = rand() % 3;
  
  int n = tmp.nrow();
  int m = tmp.ncol();

  for (int i = 0; i < n; i++){
    for (int j = 1; j < m; j++){
      if(source(i,j) != 0){
        int a = (i + rand() % 3 - 1) % n; // get an integer between [-1:1]
        int b = (j + rand() % 3 - 1) % m; // get an integer between [-1:1]
        if (a == -1) a = n -1; //ugly. we are not satisfied with this...
        if (b == -1) b = m -1;
        if(tmp(a,b) == 0){ // if empty go for it!
          tmp(a,b) = 1;
          tmp(i,j) = 0;
        }
      }
    }
  }
  return tmp;
'