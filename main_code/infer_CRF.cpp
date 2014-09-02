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

using namespace std;

// g++-4.8 -O3 -std=c++11 infer.cpp -Ieigen313/ -o infer

// plan:
// new read_data function for x only
// then read in W
// then do W2theta
// then call inference
// then return the marginals
//
// all this should be vectorized to work over many datums

int main(int argc, char * argv[]){

  int where_m = 0;
  while(where_m < argc && strcmp(argv[where_m],"-m")!=0 )
    where_m++;
  if(where_m == argc)
    throw MyException("Error! no model flag found! (-m)");

  int where_d = 0;
  while(where_d < argc && strcmp(argv[where_d],"-d")!=0 )
    where_d++;
  if(where_d == argc)
    throw MyException("Error! no data flag found! (-d)");

  int where_w = 0;
  while(where_w < argc && strcmp(argv[where_w],"-w")!=0 )
    where_w++;
  if(where_w == argc)
    throw MyException("Error! no weights flag found! (-w)");

  int where_mu = 0;
  while(where_mu < argc && strcmp(argv[where_mu],"-mu")!=0 )
    where_mu++;
  if(where_mu == argc)
    throw MyException("Error! no marginals flag found! (-mu)");

  int where_i = 0;
  int niters = 10;
  while(where_i < argc && strcmp(argv[where_i],"-i")!=0)
    where_i++;
  if(where_i >= argc-1){
    cout << "using default of " << niters << " iters" << endl;
  }else{
    niters = stoi(argv[where_i+1]);
    cout << "using " << niters << " iters" << endl;
  }

  //string where = "/Volumes/ramdisk/";
  //string where = "./";

  std::vector<MatrixXd> x;
  read_data(argv[where_d+1], x);

  /*
  for(int i=0; i<x.size(); i++)
    cout << "x["<<i<<"]: " << x[i].transpose() << endl;
  */

  std::vector<MatrixXi> cliques; int nnodes; MatrixXi nvals; MatrixXd ent; MatrixXi cliquetype;
  tie(cliques, nnodes, nvals, ent, cliquetype) = read_model(argv[where_m+1]);
  Messages m = Messages(cliques, nvals, nnodes, ent, cliquetype);
    
  //cout << ent << endl;

  auto W = read_params(argv[where_w+1]);

  // check that W maps to correct sizes
  for(int c=0; c<cliques.size(); c++){
    int ctype = cliquetype(c);
    if(ctype < 0 || ctype >= W.size())
      throw MyException("Clique type outside of valid range.  (Must be between 0 and # weight matrices)");
  }

  // check that data aligns with W and factors
  if(x.size()!=cliques.size())
    throw MyException("Error: Input data has different number of cliques than model.");
  for(int i=0; i<x.size(); i++){
    int ctype = cliquetype(i);
    if(x[i].size() != W[ctype].cols())
      throw MyException("Error: size of data for factor " + to_string(i) + " does not match size of corresponding weight vector (W[" + to_string(ctype) + "])");
  }


  /*
  cout << "here is W:" << endl;
  cout << W.size() << endl;
  cout << W[0] << endl << endl;
  cout << W[1] << endl;
  */

  if(m.cliques.size() != x.size())
    cout << "ERROR!: size of input does not match model!" << endl;

  auto theta = W2theta(W, m, x);

  /*
  cout << "here is theta:" << endl;
  cout << theta[0] << endl;
  */

  //auto theta = read_theta("theta.txt");

  //test_config_tforms(m);
  //test_indexing(m);
  //test_msg_indexing(m);

  auto mu = infer_parchild(m,theta,ent,niters);
  write_marginals(argv[where_mu+1],mu);

  return 0;
}
