void read_command_line(int argc, char * argv[], int &where_m, int &where_d, int &ndata, int & niters, string & opt_alg, MatrixXd & opt_params, bool & have_W, int & where_w, int & where_wout, double & reg, bool doprint){
  //cout << "argc: " << argc << endl;
  where_m = 0;
  while(where_m < argc && strcmp(argv[where_m],"-m")!=0 )
    where_m++;
  if(where_m == argc)
    throw MyException("Error! no model flag found!");
  int where_m_end = where_m+1;
  while(where_m_end < argc && argv[where_m_end][0] != '-')
    where_m_end++;

  where_d = 0;
  while(where_d < argc && strcmp(argv[where_d],"-d")!=0 )
    where_d++;
  if(where_d == argc)
    throw MyException("Error! no data flag found!");
  int where_d_end = where_d+1;
  while(where_d_end < argc && argv[where_d_end][0] != '-')
    where_d_end++;

  int where_i = 0;
  niters = 10;
  while(where_i < argc && strcmp(argv[where_i],"-i")!=0)
    where_i++;
  if(where_i >= argc-1){
    if(doprint)
      cout << "using default of " << niters << " iters" << endl;
  }else{
    niters = stoi(argv[where_i+1]);
    if(doprint)
      cout << "using " << niters << " iters" << endl;
  }

  int where_r = 0;
  reg = .01;
  while(where_r < argc && strcmp(argv[where_r],"-r")!=0)
    where_r++;
  if(where_r >= argc-1){
    if(doprint)
      cout << "using default regularizer of " << reg << endl;
  }else{
    reg = stod(argv[where_r+1]);
    if(doprint)
      cout << "using regularizer of " << reg << endl;
  }


  int where_a = 0;
  opt_alg = "lbfgs";
  while(where_a < argc && strcmp(argv[where_a],"-a")!=0)
    where_a++;
  if(where_a >= argc-1){
    if(doprint)
      cout << "using default algorithm of " << opt_alg << endl;
  }else{
    opt_alg = argv[where_a+1];
    if(doprint)
      cout << "using algorithm: " << opt_alg << endl;
    // find next (or when)
    int where_a_end = where_a+1;
    while(where_a_end < argc && argv[where_a_end][0] != '-')
      where_a_end++;
    int len = where_a_end-where_a-2;
    opt_params = MatrixXd(len,1);
    for(int i=0; i<len; i++)
      opt_params(i) = strtod(argv[where_a+2+i],NULL);
    if(doprint)
      cout << "opt_params: " << opt_params.transpose() << endl;
  }

  ndata   = where_d_end-where_d-1;
  int nmodels = where_m_end-where_m-1;
  
  if(ndata != nmodels)
    throw MyException("Error! number of data not equal to number of models!");

  /*
  cout << "here are the model files: ";
  for(int i=where_m+1; i < where_m_end; i++)
    cout << argv[i] << " ";
  cout << endl;
  cout << "here are the data files: ";
  for(int i=where_d+1; i < where_d_end; i++)
    cout << argv[i] << " ";
  cout << endl;
  */

  where_w = 0;
  while(where_w < argc && strcmp(argv[where_w],"-w")!=0 )
    where_w++;
  if(where_w >= argc-1){
    have_W = false;
    if(doprint)
      cout << "using default (zero) weights" << endl;
  }
  else{
    have_W = true;
    if(doprint)
      cout << "reading W from " << argv[where_w+1] << endl;
  }

  where_wout = 0;
  while(where_wout < argc && strcmp(argv[where_wout],"-wout")!=0 )
    where_wout++;
  if(where_wout >= argc-1){
    if(doprint)
      cout << "writing W to default of W.txt" << endl;
  }
  else{
    if(doprint)
      cout << "writing W to " << argv[where_wout+1] << endl;
  }


}

std::vector<MatrixXd> read_theta(string fname){
  // read in theta
  int i;
  ifstream ftheta(fname);
  if (!ftheta) 
    throw MyException("There was a problem opening " + fname + " for reading.");
  
  std::vector<MatrixXd> theta;
  while(ftheta >> i){
    MatrixXd mytheta(i,1);
    for(int n=0; n<i; n++)
      ftheta >> mytheta(n);
    theta.push_back(mytheta);
  }
  return theta;
}

void sanity_checks(std::vector<MatrixXi> &cliques,int nnodes){
  // do some checks: min node is zero, etc.
  MatrixXi node_counts = MatrixXi::Zero(nnodes,1);
  for(int c=0; c<cliques.size(); c++){
    for(int n=0; n<cliques[c].size(); n++){
      int i = cliques[c](n);
      if(i < 0 || i >= nnodes)
	throw MyException("ERROR: node numbers should be between 0 and nnodes-1.  (Found node " + to_string(i) + " with nnodes=" + to_string(nnodes) + ")");
      node_counts(i)++;
    }
  }
  if(node_counts.minCoeff() == 0)
    throw MyException("ERROR: all nodes must appear in one clique");
}

tuple<std::vector<MatrixXi>, int, MatrixXi, MatrixXd, MatrixXi> read_model(string fname){
  int i;
  ifstream infile(fname);
  
  int nnodes = 0;
  if (!infile)
    throw MyException("There was a problem opening " + fname + "for reading.");
  
  // read in first line (nnodes and nvals)
  infile >> nnodes;
  //nvals = MatrixXi(nnodes,1);
  MatrixXi nvals(nnodes,1);
  for(int i=0; i<nnodes; i++)
    infile >> nvals(i);

  if(nvals.minCoeff()<0)
    throw MyException("Error: model nvals must all be non-negative integers");

  // read in entropy values
  int ncliques = 0;
  infile >> ncliques;
  MatrixXd ent(ncliques,1);
  for(int i=0; i<ncliques; i++)
    infile >> ent(i);

  //cout << "nnodes: " << nnodes << "  ncliques: " << ncliques << endl;

  // read in clique types
  int ncliques2 = 0;
  infile >> ncliques2;
  assert(ncliques2==ncliques);
  MatrixXd cliquetype0(ncliques,1);
  for(int i=0; i<ncliques; i++)
    infile >> cliquetype0(i);
  MatrixXi cliquetype = cliquetype0.cast<int>();
  //cout << "here is cliquetype: " << cliquetype.transpose() << endl;

  // do some sanity checks on clique types
  if(cliquetype.minCoeff()<0)
    throw MyException("Error: model factor types must all be non-negative integers");
  int ntypes = cliquetype.maxCoeff()+1;
  if(ntypes>100)
    cout << "WARNING: model is specified with 100 types.  This is OK in principle, if it is what you really intend." << endl;
  MatrixXi typecounts = MatrixXi::Zero(ntypes,1);
  for(int i=0; i<cliquetype.size(); i++)
    typecounts(cliquetype(i))++;
  for(int i=0; i<typecounts.size(); i++)
    if(typecounts(i)==0){
      cout << "WARNINING: factor type " << i << " has no elements.  This is OK in principle if it is what you want" << endl;
      cout << "You would only typically want this if some factor types only occur in a few graphs" << endl;
    }
  

  // read in other lines (clique structure)
  std::vector<MatrixXi> cliques;
  while(infile >> i){
    MatrixXi mynodes(i,1);
    for(int n=0; n<i; n++)
      infile >> mynodes(n);
    cliques.push_back(mynodes);
  }

  if(cliques.size() != ncliques)
    throw MyException("Error: " + fname + " is not self-consistent in terms of the number of cliques.");

  sanity_checks(cliques,nnodes);
  return make_tuple(cliques, nnodes, nvals, ent, cliquetype);
}

/*
std::vector<tuple<std::vector<MatrixXi>, int, MatrixXi, MatrixXd, MatrixXi>> read_models(string fname){
  int i;
  ifstream infile(fname);  

  int nmodels;
  if (!infile)
    throw MyException("There was a problem opening " + fname + "for reading.");
  
  std::vector<tuple<std::vector<MatrixXi>, int, MatrixXi, MatrixXd, MatrixXi>> models;

  infile >> nmodels;
  //cout << "nmodels: " << nmodels << endl;
  for(int n=0; n<nmodels; n++){
    // read in first line (nnodes and nvals)
    int nnodes;
    infile >> nnodes;
    //cout << "nnodes: " << nnodes << endl;
    //nvals = MatrixXi(nnodes,1);
    MatrixXi nvals(nnodes,1);
    for(int i=0; i<nnodes; i++)
      infile >> nvals(i);
    //cout << "nvals: " << nvals.transpose() << endl;
    
    // read in entropy values
    int ncliques = 0;
    infile >> ncliques;
    //cout << "ncliques: " << ncliques << endl;
    MatrixXd ent(ncliques,1);
    for(int i=0; i<ncliques; i++)
      infile >> ent(i);
    
    // read in clique types
    int ncliques2 = 0;
    infile >> ncliques2;
    //cout << "ncliques2: " << ncliques2 << endl;
    assert(ncliques2==ncliques);
    MatrixXd cliquetype0(ncliques,1);
    for(int i=0; i<ncliques; i++)
      infile >> cliquetype0(i);
    MatrixXi cliquetype = cliquetype0.cast<int>();
    //cout << "here is cliquetype: " << cliquetype.transpose() << endl;
    
    //while(infile >> i){
    // read in other lines (clique structure)
    std::vector<MatrixXi> cliques;
    for(int v=0; v<ncliques; v++){
      infile >> i;
      MatrixXi mynodes(i,1);
      for(int n=0; n<i; n++)
	infile >> mynodes(n);
      cliques.push_back(mynodes);
    }
    sanity_checks(cliques,nnodes);
    models.push_back(make_tuple(cliques, nnodes, nvals, ent, cliquetype));
  }
  return models;
}
*/

void write_marginals(string fname, std::vector<MatrixXd> &mu){
  ofstream mufile(fname);
  for(int alpha=0; alpha<mu.size(); alpha++){
    mufile << mu[alpha].size() << " ";
    for(int n=0; n<mu[alpha].size(); n++)
      mufile << mu[alpha](n) << " ";
    mufile << endl;
  }
}

void write_gradient(double L, std::vector<MatrixXd> & dtheta){
  ofstream ofile("results.txt");
  ofile.precision(25);
  ofile << L << endl;
  for(int alpha=0; alpha < dtheta.size(); alpha++){
    ofile << dtheta[alpha].size();
    for(int i=0; i<dtheta[alpha].size(); i++)
      ofile << " " << dtheta[alpha](i);
    ofile << endl;
  }
  return;
}


void read_datum(string fname, std::vector<MatrixXd> & x, MatrixXi & y){
    int ndata;
    int ncliques;
    int nfeats;
    int nnodes;
    ifstream ftheta(fname);
    if (!ftheta) 
      throw MyException("There was a problem opening " + fname + " for reading.");
    
    ftheta >> ncliques;
    //cout << "ncliques[" << i << "]: " << ncliques << endl;
    for(int j=0; j<ncliques; j++){
      ftheta >> nfeats;
      x.push_back(MatrixXd(nfeats,1));
      for(int k=0; k<nfeats; k++)
	ftheta >> x[j](k);
    }
    ftheta >> nnodes;
    //cout << "nnodes[" << i << "]: " << nnodes << endl;
    y = MatrixXi(nnodes,1);
    for(int j=0; j<nnodes; j++)
      ftheta >> y(j);
}

// this one just reads a single datum
void read_data(string fname, std::vector<MatrixXd> & x){
  //x = std::vector<MatrixXd>();

    int ncliques;
    int nfeats;
    int nnodes;
    ifstream ftheta(fname);
    if (!ftheta) 
      throw MyException("There was a problem opening " + fname + " for reading.");
    
    ftheta >> ncliques;
    for(int j=0; j<ncliques; j++){
      ftheta >> nfeats;
      x.push_back(MatrixXd(nfeats,1));
      for(int k=0; k<nfeats; k++)
	ftheta >> x[j](k);
    }
}


/*
// read CRF data from a file
// the format of the file is
// #data
// #cliques in first datum
// #features in first clique, feat1, feat2, ..., featN
// #features in 2nd clique, feat1, feat2, ..., featN
// ...
// #features in Nth clique, feat1, feat2, ..., featN
// #output vars, y1, y2, ..., yN
// #cliques in 2nd datum
// ...
void read_data(string fname, std::vector<std::vector<MatrixXd>> & x, std::vector<MatrixXi> & y){
    x = std::vector<std::vector<MatrixXd>>();
    y = std::vector<MatrixXi>();

    int ndata;
    int ncliques;
    int nfeats;
    int nnodes;
    ifstream ftheta(fname);
    if (!ftheta) 
      throw MyException("There was a problem opening " + fname + " for reading.");
    
    ftheta >> ndata;
    cout << "ndata: " << ndata << endl;
    for(int i=0; i<ndata; i++){
      ftheta >> ncliques;
      //cout << "ncliques[" << i << "]: " << ncliques << endl;
      x.push_back(std::vector<MatrixXd>());
      for(int j=0; j<ncliques; j++){
	ftheta >> nfeats;
	//if(j % 1000 == 1)
	//  cout << "nfeats[" << j << "]: " << nfeats << endl;
	x[i].push_back(MatrixXd(nfeats,1));
	for(int k=0; k<nfeats; k++)
	  ftheta >> x[i][j](k);
	//if(j % 100 == 1)
	//  cout << "x[" << i << "][" << j << "]: " << x[i][j].transpose() << endl;
      }
      ftheta >> nnodes;
      //cout << "nnodes[" << i << "]: " << nnodes << endl;
      y.push_back(MatrixXi(nnodes,1));
      for(int j=0; j<nnodes; j++)
	ftheta >> y[i](j);
      //if(y[i](0)<0 || y[i](0)>1)	

      //cout << "y[" << i << "]: " << y[i].transpose() << endl;      
      //cout << "nnodes: " << nnodes << "  y min/max: " << y[i].minCoeff() << "/" << y[i].maxCoeff() << endl;

    }
    //while(ftheta >> i){  
}
*/

// read CRF data from a file
// the format of the file is
// #data
// #cliques in first datum
// #features in first clique, feat1, feat2, ..., featN
// #features in 2nd clique, feat1, feat2, ..., featN
// ...
// #features in Nth clique, feat1, feat2, ..., featN
// #output vars, y1, y2, ..., yN
// #cliques in 2nd datum
// ...
/*
void read_data(string fname, std::vector<std::vector<MatrixXd>> & x){
    x = std::vector<std::vector<MatrixXd>>();

    int ndata;
    int ncliques;
    int nfeats;
    int nnodes;
    ifstream ftheta(fname);
    if (!ftheta) 
      throw MyException("There was a problem opening " + fname + " for reading.");
    
    ftheta >> ndata;
    for(int i=0; i<ndata; i++){
      ftheta >> ncliques;
      x.push_back(std::vector<MatrixXd>());
      for(int j=0; j<ncliques; j++){
	ftheta >> nfeats;
	x[i].push_back(MatrixXd(nfeats,1));
	for(int k=0; k<nfeats; k++)
	  ftheta >> x[i][j](k);
	//cout << "x[" << i << "][" << j << "]: " << x[i][j].transpose() << endl;
      }
    }
    }*/




void write_params(string fname, const std::vector<MatrixXd> & W){
  ofstream of(fname);
  if(!of)
    throw MyException("There was a problem opening " + fname + " for writing.");
  
  int N = W.size();  

  of << N << endl;
  for(int n=0; n<N; n++){
    of << W[n].rows() << " " << W[n].cols() << endl;
    for(int j=0; j<W[n].cols(); j++)
      for(int i=0; i<W[n].rows(); i++){
	of << W[n](i,j) << " ";
     of << endl;
    }
  }
}

std::vector<MatrixXd> read_params(string fname){
  ifstream ifi(fname);
  if(!ifi)
    throw MyException("There was a problem opening " + fname + " for reading.");
  
  int N, rows, cols;
  ifi >> N;

  std::vector<MatrixXd> W;

  for(int n=0; n<N; n++){
    ifi >> rows;
    ifi >> cols;
    W.push_back(MatrixXd(rows,cols));
    for(int j=0; j<W[n].cols(); j++)
      for(int i=0; i<W[n].rows(); i++){
	ifi >> W[n](i,j);
    }
  }
  return W;
}
