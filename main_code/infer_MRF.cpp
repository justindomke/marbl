#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <map>
#include <unordered_map>
#include "math.h"
#include <tuple>
#include <stack>

#define DEBUG
//#define NDEBUG

#include "eigenstuff.cpp"
#include "util.cpp"
#include "io.cpp"
#include "Messages.cpp"
#include "testing.cpp"
#include "parchild.cpp"
#include "opt.cpp"
#include "CRF.cpp"

using namespace std;

int main(int argc, char * argv[]){

  int where_m = 0;
  while(where_m < argc && strcmp(argv[where_m],"-m")!=0 )
    where_m++;
  if(where_m == argc)
    throw MyException("Error! no model flag found! (-m)");

  int where_f = 0;
  while(where_f < argc && strcmp(argv[where_f],"-f")!=0 )
    where_f++;
  if(where_f == argc)
    throw MyException("Error! no factors flag found! (-f)");

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

  std::vector<MatrixXi> cliques; int nnodes; MatrixXi nvals; MatrixXd ent; MatrixXi cliquetype;
  tie(cliques, nnodes, nvals, ent, cliquetype) = read_model(argv[where_m+1]);
  Messages m = Messages(cliques, nvals, nnodes, ent, cliquetype);

  auto theta = read_theta(argv[where_f+1]);

  //test_config_tforms(m);
  //test_indexing(m);
  //test_msg_indexing(m);

  auto mu = infer_parchild(m,theta,ent,niters);
  write_marginals(argv[where_mu+1],mu);

  return 0;
}
