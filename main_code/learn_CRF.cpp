#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <map>
#include <unordered_map>
#include "math.h"
#include <tuple>
#include <stack>

//#define DEBUG
#define NDEBUG

#include "eigenstuff.cpp"
#include "util.cpp"
#include "io.cpp"
#include "Messages.cpp"
#include "testing.cpp"
#include "parchild.cpp"
#include "opt.cpp"
#include "CRF.cpp"

// TODO:
// write different datums and models to separate files, so writing/reading can be multithreaded
// alternative binary format
// stochastic gradient training

using namespace std;

int main(int argc, char * argv[]){
  Eigen::initParallel();

  // parse the command line options.  plan is:
  // 1) find the first and last location of the models
  // 2) find the first and last location of the data

  cout << "argc: " << argc << endl;
  int where_m = 0;
  while(where_m < argc && strcmp(argv[where_m],"-m")!=0 )
    where_m++;
  if(where_m == argc)
    throw MyException("Error! no model flag found!");
  int where_m_end = where_m+1;
  while(where_m_end < argc && argv[where_m_end][0] != '-')
    where_m_end++;

  int where_d = 0;
  while(where_d < argc && strcmp(argv[where_d],"-d")!=0 )
    where_d++;
  if(where_d == argc)
    throw MyException("Error! no data flag found!");
  int where_d_end = where_d+1;
  while(where_d_end < argc && argv[where_d_end][0] != '-')
    where_d_end++;

  int where_i = 0;
  int niters = 10;
  while(where_i < argc && strcmp(argv[where_i],"-i")!=0)
    where_i++;
  if(where_i >= argc-1)
    cout << "using default of " << niters << " iters" << endl;
  else{
    niters = stoi(argv[where_i+1]);
    cout << "using " << niters << " iters" << endl;
  }

  int where_a = 0;
  string opt_alg = "lbfgs";
  MatrixXd opt_params;
  while(where_a < argc && strcmp(argv[where_a],"-a")!=0)
    where_a++;
  if(where_a >= argc-1)
    cout << "using default algorithm of " << opt_alg << endl;
  else{
    opt_alg = argv[where_a+1];
    cout << "using algorithm: " << opt_alg << endl;
    // find next (or when)
    int where_a_end = where_a+1;
    while(where_a_end < argc && argv[where_a_end][0] != '-')
      where_a_end++;
    int len = where_a_end-where_a-2;
    opt_params = MatrixXd(len,1);
    for(int i=0; i<len; i++)
      opt_params(i) = strtod(argv[where_a+2+i],NULL);
    cout << "opt_params: " << opt_params.transpose() << endl;
  }
  
  int where_w = 0;
  bool have_W;
  while(where_w < argc && strcmp(argv[where_w],"-w")!=0 )
    where_w++;
  if(where_w >= argc-1){
    have_W = false;
    cout << "using default (zero) weights" << endl;
  }
  else{
    have_W = true;
    cout << "reading W from " << argv[where_w+1] << endl;
  }

  //cout << "where_m: " << where_m << " where_m_end: " << where_m_end << endl;
  //cout << "where_d: " << where_d << " where_d_end: " << where_d_end << endl;
  cout << "here are the model files: ";
  for(int i=where_m+1; i < where_m_end; i++)
    cout << argv[i] << " ";
  cout << endl;
  cout << "here are the data files: ";
  for(int i=where_d+1; i < where_d_end; i++)
    cout << argv[i] << " ";
  cout << endl;

  int ndata   = where_d_end-where_d-1;
  int nmodels = where_m_end-where_m-1;
  
  if(ndata != nmodels)
    throw MyException("Error! number of data not equal to number of models!");

  std::vector<std::vector<MatrixXd>> x(ndata);
  std::vector<MatrixXi> y(ndata);
  cout << "reading data...";
  cout.flush();
  #pragma omp parallel for
  for(int i=0; i<ndata; i++)
    read_1data(argv[where_d+1+i],x[i],y[i]);
  cout << "done." << endl;
      
  //int len = y[1].rows();

  cout << "reading models and creating message structures... "; 
  cout.flush();
  std::vector<Messages> m(nmodels);
  #pragma omp parallel for
    for(int n=0; n<nmodels; n++){
      std::vector<MatrixXi> cliques; int nnodes; MatrixXi nvals; MatrixXd ent; MatrixXi cliquetype;
      tie(cliques, nnodes, nvals, ent, cliquetype) = read_model(argv[where_m+1+n]);
      m[n] = Messages(cliques, nvals, nnodes, ent, cliquetype);      
      //m.push_back(Messages(cliques, nvals, nnodes, ent, cliquetype));
    }
    cout << "done." << endl;

    // do checks on y
    for(int n=0; n<m.size(); n++)
      for(int i=0; i<m[n].nnodes; i++)
	assert(y[n](i) >= -1 && y[n](i) < m[n].nvals(i));

    // do checks on m
    for(int n=0; n<m.size(); n++){
      for(int c=0; c<m[n].cliques.size(); c++){
	assert(m[n].compute_nconfigs(c) > 0);
	assert(m[n].nconfigs(c) > 0);
      }
    }
      

    cout << "messages built" << endl;

    std::vector<MatrixXd> W;
    if(have_W){
      W = read_params(argv[where_w+1]);
      W = fit_CRF(x, y, m, opt_alg, opt_params, niters, W);
    } else
      W = fit_CRF(x, y, m, opt_alg, opt_params, niters);
    
    write_params("W.txt",W);

  //cout << ent << endl;

  //auto theta = read_theta("theta.txt");

  //test_config_tforms(m);
  //test_indexing(m);
  //test_msg_indexing(m);

  //auto mu = infer_parchild(m,theta,ent);
  //write_marginals(mu);

  // re-initialize messages
  //m = Messages(cliques, nvals, nnodes);

  //double L;
  //std::vector<MatrixXd> dtheta;
  //infer_parchild_bprop(m,theta,L,dtheta);  
  //write_gradient(L, dtheta);

  /*
  std::vector<MatrixXd> x;
  MatrixXi y;
  read_datum(
  */

  return 0;
}
