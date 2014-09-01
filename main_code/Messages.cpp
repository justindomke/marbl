class Messages{
public:
  Messages(const std::vector<MatrixXi> &cliques_in, const MatrixXi &nvals_in, const int nnodes_in, const MatrixXd & ent, const MatrixXi & cliquetype);
  Messages();
  ~Messages();
  Messages(const Messages& m0); // Copy Constructor
  Messages& operator=( const Messages& m0 );
  void initialize(const std::vector<MatrixXi> &cliques_in, const MatrixXi &nvals_in, const int nnodes_in, const MatrixXd & ent, const MatrixXi & cliquetype);
  MatrixXd& msg(int, int);
  double& msg(int, int, MatrixXi&);
  double& msg(int, int, int);
  double& msg(int, int, int, int);
  double& par_msg(int, int, MatrixXi&);
  double& par_msg(int, int, int);
  int kid_config(int, int, int, int);
  int kid_config(int, int, int);
  int find_kid(int, int);
  const std::vector<int>& get_pars(int);
  const std::vector<int>& get_kids(int);
  int vec2config(int,MatrixXi);
  MatrixXi config2vec(int,int);
  std::vector<MatrixXi> cliques;
  std::vector<MatrixXd> *lambda;
  void print_messages();
  MatrixXi nkids;
  MatrixXi npars;
  int nconfigs(int);
  MatrixXi nvals;
  int      nnodes;
  MatrixXd ent;
  MatrixXi cliquetype;
  int compute_nconfigs(int c);
  int *nconfigs_mat;
private:
  std::vector<int> *pars;
  std::vector<int> *kids;
  int compute_kid_config(int, int, int, int);
  std::vector<MatrixXi> *kid_configs;
  std::map<int,int> *kid_finder;
  int ncliques;
};

Messages& Messages::operator=( const Messages& m0 ) {
  //cout << "assignment operator called" << endl;
  this->initialize(m0.cliques, m0.nvals, m0.nnodes, m0.ent, m0.cliquetype);
  return *this;
}


Messages::Messages(const Messages& m0){ // Copy Constructor
   //counts = new unsigned int[a.size] (); // Add this to allocate memory for counts
  //Messages myobj = 
  //*this = Messages
  cout << "copy constructor called!" << endl;
  this->initialize(m0.cliques, m0.nvals, m0.nnodes, m0.ent, m0.cliquetype);
}

// dummy construactor
Messages::Messages(){
};

inline int Messages::nconfigs(int c){
  assert(0 <= c && c <= ncliques);
  return nconfigs_mat[c];
}

int Messages::compute_nconfigs(int c){
  int mult = 1;
  //cout << "cliquesize: " << cliques[c].size() << endl;
  assert(cliques[c].size() > 0 );

  for(int n=0; n<cliques[c].size(); n++){
    int i = cliques[c](n);
    mult *= nvals(i);
  }
  assert(mult>0);
  return mult;
}

int Messages::vec2config(int c, MatrixXi vec){
  int conf = 0;
  int mult = 1;
  if(vec.size() != cliques[c].size())
    throw MyException("vec size incorrect in vec2config");

  for(int n=0; n<cliques[c].size(); n++){
    int i = cliques[c](n);
    if(vec(n)<0 || vec(n) >= nvals(i)){
      cout << "c: " << c << " cliques(c): " << cliques[c].transpose() << " vec: " << vec.transpose() << " len: " << vec.size() << endl;
      throw MyException("incorrect vec argument in vec2config");
    }

    conf += vec(n)*mult;
    mult *= nvals(i);
  }
  return conf;
}

MatrixXi Messages::config2vec(int c, int config){
  MatrixXi vec = MatrixXi::Zero(cliques[c].size(),1);
  int mult = nconfigs(c);
  for(int n = cliques[c].size()-1; n>=0; n--){
    int i   = cliques[c](n);
    mult    = mult / nvals(i);
    vec(n)  = config / mult;
    config -= mult*vec(n);
    if(vec(n)<0 || vec(n) >= nvals(i))
      throw MyException("incorrect config argument in config2vec");
  }
  return vec;
}

void Messages::initialize(const std::vector<MatrixXi> &cliques_in, const MatrixXi &nvals_in, const int nnodes_in, const MatrixXd & ent_in, const MatrixXi & cliquetype_in){
  cliques = cliques_in;
  nvals   = nvals_in;
  nnodes  = nnodes_in;
  ent     = ent_in;
  cliquetype = cliquetype_in;

  ncliques = cliques.size();  

  nconfigs_mat = new int[ncliques];
  for(int c=0; c<ncliques; c++)
      nconfigs_mat[c] = compute_nconfigs(c);

  // compute all parents and children
  //
  // the naive way would be just to check all pairs of cliques, but this would be slow
  // instead, check all pairs that intersect on at least 1 node

  std::vector<std::vector<int>> i2c(nnodes);
  for(int c=0; c<cliques.size(); c++)
    for(int n=0; n<cliques[c].size(); n++){
      i2c[cliques[c](n)].push_back(c);
    }
    
  /*
  for(int i=0; i<nnodes; i++){
    cout << "node " << i << " cliques: ";
    for(int j=0; j<i2c[i].size(); j++)
        cout << i2c[i][j] << " ";
    cout << endl;
    }*/

  pars = new std::vector<int>[ncliques];
  kids = new std::vector<int>[ncliques];

  for(int i=0; i<nnodes; i++){
    for(int n=0; n<i2c[i].size(); n++){
      for(int m=n+1; m<i2c[i].size(); m++){
	int c = i2c[i][n];
	int d = i2c[i][m];
      
	//cout << cliques[c] << endl << "---" << endl << cliques[d] << endl << "rez=" << contains(cliques[c],cliques[d]) << endl;

	if(contains(cliques[c],cliques[d])){
	  if(!contains(kids[c],d))
	    kids[c].push_back(d);
	  if(!contains(pars[d],c))
	    pars[d].push_back(c);
	}
	
	if(contains(cliques[d],cliques[c])){
	  if(!contains(kids[d],c))
	    kids[d].push_back(c);
	  if(!contains(pars[c],d))
	    pars[c].push_back(d);	  
	}
      }
    }
  }

  nkids = MatrixXi(cliques.size(),1);
  npars = MatrixXi(cliques.size(),1);
  // populate nkids and npars
  for(int c=0; c<cliques.size(); c++){
    nkids(c) = kids[c].size();
    npars(c) = pars[c].size();
  }

  // put everything into kid_finder
  kid_finder = new map<int,int>[ncliques];
  for(int alpha=0; alpha<cliques.size(); alpha++){
    for(int n=0; n<get_kids(alpha).size(); n++){
      int beta = get_kids(alpha)[n];
      kid_finder[alpha][beta] = n;
    }
  }

  // precompute configuration mappers
  kid_configs = new std::vector<MatrixXi>[ncliques];
  for(int alpha=0; alpha<ncliques; alpha++){
    for(int n=0; n<nkids(alpha); n++){
      int beta = kids[alpha][n];
      kid_configs[alpha].push_back(MatrixXi(nconfigs(alpha),1));
      for(int i=0; i<nconfigs(alpha); i++)
	kid_configs[alpha][n](i) = compute_kid_config(alpha,beta,i,n);
    }
  }

  /*
  for(int c=0; c<cliques.size(); c++){
    cout << "pars[" << c << "] ";
    for(int n=0; n<pars[c].size(); n++)
      cout << " " << pars[c][n];
    cout << endl;
  }

  for(int c=0; c<cliques.size(); c++){
    cout << "kids[" << c << "] ";
    for(int n=0; n<kids[c].size(); n++)
      cout << " " << kids[c][n];
    cout << endl;
    }*/

  // check that there arent two equivalent regions that are both parents of each other
  for(int c=0; c<cliques.size(); c++){
    for(int n=0; n<get_pars(c).size(); n++){
      int d = get_pars(c)[n];
      if(contains(get_pars(d),c)){
	cout << "ERROR: regions " << c << " and " << d << " contain each other" << endl;
	throw MyException("ERROR: regions contain each other");
      }
      assert(!contains(get_pars(d),c));
    }
  }

  // create messages
  lambda = new std::vector<MatrixXd>[cliques.size()];
  //lambda = new std::vector<double*>[cliques.size()];
  for(int c=0; c<cliques.size(); c++)
    for(int n=0; n<get_kids(c).size(); n++){
      int b = kids[c][n];
      lambda[c].push_back(MatrixXd::Zero(nconfigs(b),1));
      //cout << "lambda[" << c << "] to " << b << "  has size " << lambda[c][n].size() << endl;
      //lambda[c].push_back(new double[nconfigs(b)]);
    }
}

Messages::Messages(const std::vector<MatrixXi> &cliques_in, const MatrixXi &nvals_in, const int nnodes_in, const MatrixXd & ent_in, const MatrixXi & cliquetype_in){
  initialize(cliques_in, nvals_in, nnodes_in, ent_in, cliquetype_in);
}

Messages::~Messages(){
  //nconfigs_mat = new int[ncliques];
  //pars = new std::vector<int>[ncliques];
  //kids = new std::vector<int>[ncliques];
  //kid_finder = new map<int,int>[ncliques];
  //kid_configs = new std::vector<MatrixXi>[ncliques];
  //lambda = new std::vector<MatrixXd>[cliques.size()];
  delete[] nconfigs_mat;
  delete[] pars;
  delete[] kids;
  delete[] kid_finder;
  delete[] kid_configs;
  delete[] lambda;
  //cout << "deleting... " << endl;
}  

const std::vector<int>& Messages::get_pars(int c){
  assert(c>=0 && c < ncliques);
  return pars[c];
}

const std::vector<int>& Messages::get_kids(int c){
  assert(c>=0 && c < ncliques);
  return kids[c];
}

int Messages::find_kid(int alpha, int beta){
  /*
  int n=0;
  while(n<get_kids(alpha).size()){
    if(get_kids(alpha)[n]==beta)
      break;
    n++;
  }
  */
  int n2 = kid_finder[alpha][beta];
  //assert(n < get_kids(alpha).size());
  //assert(n2 == n);
  assert(kids[alpha][n2]==beta);
  return n2;
}

MatrixXd& Messages::msg(int alpha, int beta){
  //assert(beta==-1);
  int n = find_kid(alpha,beta);
  return lambda[alpha][n];
}


double& Messages::msg(int alpha, int beta, MatrixXi& vec){
  int n = find_kid(alpha,beta);
  // need to get config in beta corresponding to vals
  MatrixXi vec2 = corresponding_indices(vec,cliques[alpha],cliques[beta]);
  int config2   = vec2config(beta,vec2);
  return lambda[alpha][n](config2);
}

int Messages::compute_kid_config(int alpha, int beta, int config0, int n){
  assert(kids[alpha][n]==beta);
  
  // decompose config2 into its values (in alpha)
  int config = config0;
  int config3 = 0;
  int mult = nconfigs(alpha);
  for(int n = cliques[alpha].size()-1; n>=0; n--){
    int i   = cliques[alpha](n);
    mult    = mult / nvals(i);
    int xi  = config / mult;
    config -= mult*xi;
    assert(xi >= 0 && xi < nvals(i));
    
    // now set m so that cliques[beta](m) == i
    // (if that is possible)
    int mult2 = 1;
    int m = 0;
    while(m < cliques[beta].size()){
      if(cliques[beta](m)==i)
	break;
      mult2 *= nvals(cliques[beta](m));
      m++;
    }
    if(m==cliques[beta].size())
      continue;
    config3 += mult2*xi;
  }  
  
  #ifdef DEBUG
  MatrixXi vec  = config2vec(alpha,config0);
  MatrixXi vec2 = corresponding_indices(vec,cliques[alpha],cliques[beta]);
  int config2   = vec2config(beta,vec2);
  assert(config2 == config3);
  #endif
  
  return config3;
}

int Messages::kid_config(int alpha, int beta, int config, int n){
  assert(kids[alpha][n]==beta);
  return kid_configs[alpha][n](config);
}

int Messages::kid_config(int alpha, int beta, int config0){
  int n = find_kid(alpha,beta);
  return kid_config(alpha, beta, config0, n);
}

double& Messages::msg(int alpha, int beta, int config, int n){
  assert(kids[alpha][n]==beta);
  int config2 = kid_config(alpha,beta,config, n);

  return lambda[alpha][n](config2);
}

double& Messages::msg(int alpha, int beta, int config){
  int n = find_kid(alpha,beta);
  /*
  MatrixXi vec  = config2vec(alpha,config);
  MatrixXi vec2 = corresponding_indices(vec,cliques[alpha],cliques[beta]);
  int config2   = vec2config(beta,vec2);
  return lambda[alpha][n](config2);
  */
  //int config2 = kid_config(alpha,beta,config);
  //return lambda[alpha][n](config2);
  return msg(alpha, beta, config, n);
}

double& Messages::par_msg(int alpha, int beta, MatrixXi& vec){
  int n      = find_kid(beta,alpha);
  int config = vec2config(alpha,vec);
  return lambda[beta][n](config);
}

double& Messages::par_msg(int alpha, int beta, int config){
  int n = find_kid(beta,alpha);
  return lambda[beta][n](config);
}

void Messages::print_messages(){
  for(int alpha=0; alpha<cliques.size(); alpha++){
    for(int n=0; n<get_kids(alpha).size(); n++){
      int beta = get_kids(alpha)[n];
      cout << "msg from " << alpha << " to " << beta << ": ";
      for(int i=0; i<nconfigs(beta); i++)
	cout << par_msg(beta,alpha,i) << " ";
      cout << endl;
    }
  }

  /*
  for(int alpha=0; alpha<cliques.size(); alpha++){
    for(int n=0; n<get_kids(alpha).size(); n++){
      int beta = get_kids(alpha)[n];
      cout << "msg from " << alpha << " to " << beta << ": ";
      cout << lambda[alpha][n].transpose();
      //for(int i=0; i<nconfigs(beta); i++)
      //cout << par_msg(beta,alpha,i) << " ";
      cout << endl;
    }
    }*/

}
