struct MyException : public std::exception{
   std::string s;
  MyException(std::string ss) : s(ss) { cout << ss << endl;}
   ~MyException() throw () {} // Updated
   const char* what() const throw() { return s.c_str(); }
};

MatrixXi corresponding_indices(MatrixXi &x_alpha, MatrixXi &alpha, MatrixXi &beta){
  MatrixXi x_beta(beta.size(),1);
  for(int n=0; n<beta.size(); n++){
    for(int m=0; m<alpha.size(); m++){
      if(beta(n)==alpha(m)){
	x_beta(n) = x_alpha(m);
	break;
      }
    }
  }
  return x_beta;
}

//template <class T>
bool contains(MatrixXi &A, int b){
  for(int i=0; i<A.size(); i++){
    if(A(i)==b)
      return true;
  }
  return false;
}

//template <class T>
bool contains(MatrixXi &A, MatrixXi &B){
  for(int i=0; i<B.size(); i++)
    if(!contains(A,B(i))){
      return false;
    }
  return true;
}

template <class T>
bool contains(std::vector<T> A, T b){
  for(int i=0; i<A.size(); i++)
    if(A[i]==b)
      return true;
  return false;
}

template <class T>
bool contains(std::vector<T> A, std::vector<T> B){
  for(int i=0; i<B.size(); i++)
    if(!contains(A,B(i)))
      return false;
  return true;
}
