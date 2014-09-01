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

  //string where = "/Volumes/ramdisk/";
  string where = "./";

  std::vector<MatrixXd> x;
  read_data(where + "data.txt", x);
  //return 0;

  //cout << "here is x:" << endl;
  //cout << x[0][0] << endl;

  //int nnodes;
  //MatrixXi nvals;
  //auto cliques = read_model(nnodes, nvals);
  std::vector<MatrixXi> cliques; int nnodes; MatrixXi nvals; MatrixXd ent; MatrixXi cliquetype;
  tie(cliques, nnodes, nvals, ent, cliquetype) = read_model(where + "model.txt");
  Messages m = Messages(cliques, nvals, nnodes, ent, cliquetype);

  //cout << ent << endl;

  auto W = read_params("Win.txt");

  // cout << "here is W:" << endl;
  // cout << W.size() << endl;
  // cout << W[0] << endl << endl;
  // cout << W[1] << endl;

  if(m.cliques.size() != x.size())
    cout << "ERROR!: size of input does not match model!" << endl;

  auto theta = W2theta(W, m, x);

  //cout << "here is theta:" << endl;
  //cout << theta[0] << endl;

  //auto theta = read_theta("theta.txt");

  //test_config_tforms(m);
  //test_indexing(m);
  //test_msg_indexing(m);

  auto mu = infer_parchild(m,theta,ent);
  write_marginals(where + "mu.txt",mu);

  return 0;
}
