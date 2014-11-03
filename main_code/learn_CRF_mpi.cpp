#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
//#include <map>
//#include <unordered_map>
#include "math.h"
#include <tuple>
#include <stack>
#include <mpi.h>
#include <sys/time.h>
#include <unistd.h>

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

// TO install liblbfgs
// download it, untar it, go to directory
// ./configure --prefix=/home/jcdicsa/.liblbfgs
// make
// make install

// IMPORTANT: export OMP_STACKSIZE=10m
// g++-4.8 -O3 -fopenmp -std=c++11 learn.cpp -llbfgs -Ieigen313/ -o learn

// LINUX COMPILE
// g++ -fopenmp -O2 -std=c++11 learn.cpp -Wl,-Bstatic -I/home/jcdicsa/.liblbfgs/include -L/home/jcdicsa/.liblbfgs/lib -llbfgs -Wl,-Bdynamic -Ieigen313/ -o learn


typedef std::vector<MatrixXd> features;

MatrixXd serialize(const std::vector<MatrixXd> A){
  // first, compute a length
  int numel = 0;
  for(int i=0; i<A.size(); i++)
    numel += A[i].rows()*A[i].cols();
  MatrixXd B(numel,1);
  int where=0;
  for(int i=0; i<A.size(); i++){
      for(int n=0; n<A[i].rows()*A[i].cols(); n++){
	B(where) = A[i](n);
	where++;
      }
    }	
  return B;
}

void deserialize(MatrixXd B, std::vector<MatrixXd> & A){
  int where = 0;
  for(int i=0; i<A.size(); i++){
      for(int n=0; n<A[i].rows()*A[i].cols(); n++){
	A[i](n) = B(where);
	where++;
      }
    }	
}

std::vector<MatrixXd> init_W(const std::vector<features> & x, const std::vector<MatrixXi> y, std::vector<Messages> & m){

  // TODO: be more careful here with computing nctypes and checking all the models/data
  int nctypes = 1 + m[0].cliquetype.maxCoeff();
  VectorXi output_sizes = VectorXi::Ones(nctypes)*-1;
  for(int c=0; c<m[0].cliques.size(); c++){
    int ctype    = m[0].cliquetype(c);
    if(output_sizes(ctype)==-1){
      output_sizes(ctype) = m[0].nconfigs(c);
    }
    assert(output_sizes(ctype)==m[0].nconfigs(c));
  }

  // do somethiing similar for input sizes  
  VectorXi input_sizes = VectorXi::Ones(nctypes)*-1;
  for(int c=0; c<m[0].cliques.size(); c++){
    int ctype    = m[0].cliquetype(c);
    if(input_sizes(ctype)==-1){
      input_sizes(ctype) = x[0][c].size();
    }
  }
  for(int n=0; n<x.size(); n++){
    for(int c=0; c<m[n].cliques.size(); c++){
      int ctype    = m[n].cliquetype(c);
      assert(input_sizes(ctype)==x[n][c].size());
    }
  }

  // now, create a bunch of arrays that map from inputs to outputs
  std::vector<MatrixXd> W;
  for(int ctype=0; ctype < nctypes; ctype++){
    W.push_back(MatrixXd::Random(output_sizes(ctype),input_sizes(ctype)));
  }
  return W;
}

//std::vector<MatrixXd> W = create_W(argc, argv, where_m, where_d, have_W, where_w);
std::vector<MatrixXd> create_W(int argc, char * argv[], int where_m, int where_d, bool have_W, int where_w, int &nctypes, VectorXi & output_sizes, VectorXi & input_sizes){
  std::vector<std::vector<MatrixXd>> x(1);
  std::vector<MatrixXi> y(1);
  std::vector<Messages> m(1);
  
  // actually create W (requires reading first datum)
  int i = 0;
  read_datum(argv[where_d+1+i],x[i],y[i]);
  int n = 0;
  std::vector<MatrixXi> cliques; int nnodes; MatrixXi nvals; MatrixXd ent; MatrixXi cliquetype;
  tie(cliques, nnodes, nvals, ent, cliquetype) = read_model(argv[where_m+1+n]);
  m[n] = Messages(cliques, nvals, nnodes, ent, cliquetype);      
  
  std::vector<MatrixXd> W;
  if(have_W)
    W = read_params(argv[where_w+1]);
  else
    W = init_W(x, y, m);


  nctypes = 1 + m[0].cliquetype.maxCoeff();
  output_sizes = VectorXi::Ones(nctypes)*-1;
  for(int c=0; c<m[0].cliques.size(); c++){
    int ctype    = m[0].cliquetype(c);
    if(output_sizes(ctype)==-1){
      output_sizes(ctype) = m[0].nconfigs(c);
    }
    assert(output_sizes(ctype)==m[0].nconfigs(c));
  }

  // do somethiing similar for input sizes  
  input_sizes = VectorXi::Ones(nctypes)*-1;
  for(int c=0; c<m[0].cliques.size(); c++){
    int ctype    = m[0].cliquetype(c);
    if(input_sizes(ctype)==-1){
      input_sizes(ctype) = x[0][c].size();
    }
  }
  for(int n=0; n<x.size(); n++){
    for(int c=0; c<m[n].cliques.size(); c++){
      int ctype    = m[n].cliquetype(c);
      assert(input_sizes(ctype)==x[n][c].size());
    }
  }

  return W;
}

void check_datum(Messages &m, MatrixXi &y){
  for(int i=0; i<m.nnodes; i++)
    assert(y(i) >= -1 && y(i) < m.nvals(i));

    // do checks on m
  for(int c=0; c<m.cliques.size(); c++){
    assert(m.compute_nconfigs(c) > 0);
    assert(m.nconfigs(c) > 0);
  }
}


 void compute_loss_cold(int argc, char * argv[], MatrixXd &W_serialized, int ncalc, int *who_to_calc, double &my_L, MatrixXd &my_dW_serialized,
			std::vector<features> &x, std::vector<MatrixXi> &y, std::vector<Messages> &m, std::vector<bool> &t, std::vector<MatrixXd> & W,
			int nctypes, VectorXi &output_sizes, VectorXi &input_sizes, double &gradnorm, int iter_mult){
   // parse the command line
   string opt_alg; MatrixXd opt_params; bool have_W; double reg;
   int where_m, where_m_end, where_d, where_d_end, niters, ndata, where_w, where_wout;
   read_command_line(argc, argv, where_m, where_d, ndata, niters, opt_alg, opt_params, have_W, where_w, where_wout, reg, false);

    // put serialized weights in
   deserialize(W_serialized,W);
   
   // create structure for dW
   auto my_dW = std::vector<MatrixXd>();
   for(int ctype=0; ctype < nctypes; ctype++)
     my_dW.push_back(MatrixXd::Zero(output_sizes(ctype),input_sizes(ctype)));

   for(int i=0; i<ncalc; i++){
     int n = who_to_calc[i];

     if(!t[n]){
      read_datum(argv[where_d+1+n],x[n],y[n]);
      
      std::vector<MatrixXi> cliques; int nnodes; MatrixXi nvals; MatrixXd ent; MatrixXi cliquetype;
      tie(cliques, nnodes, nvals, ent, cliquetype) = read_model(argv[where_m+1+n]);
      m[n] = Messages(cliques, nvals, nnodes, ent, cliquetype);      
      check_datum(m[n],y[n]);
      t[n] = true;
    }

     auto theta = W2theta(W,m[n],x[n]);
	  
     // compute configurations (to make loss simpler)
     MatrixXi y_configs = -1*MatrixXi::Ones(m[n].cliques.size(),1);
     for(int alpha=0; alpha<m[n].cliques.size(); alpha++){
       if(m[n].cliques[alpha].size()>1) continue;
       
       MatrixXi y_alpha = MatrixXi(m[n].cliques[alpha].size(),1);
       for(int jj=0; jj<y_alpha.size(); jj++)
	 y_alpha(jj) = y[n](m[n].cliques[alpha](jj));
       
       if(y_alpha.minCoeff()<0) continue; // if any variables hidden, forget it
       
       //int loc = m[n].cliques[alpha](0);
       y_configs(alpha) = m[n].vec2config(alpha,y_alpha);
     }
     
     // create loss function for this datum
     auto loss = [y_configs](const std::vector<MatrixXd> & mu, double & L, std::vector<MatrixXd> & dlogmu){
       L = 0;
       for(int alpha=0; alpha<mu.size(); alpha++){
	 if(y_configs(alpha)==-1) continue; // skip when y_configs == -1
	 
	 L -= log(mu[alpha](y_configs(alpha)));
	 dlogmu[alpha].setZero();
	 dlogmu[alpha](y_configs(alpha)) = -1;
	 
	 if(L!=L) cout << "badmu: " << mu[alpha].transpose() << endl;
       }
     };
     
     std::vector<MatrixXd> dtheta;
     
     double mini_L = 0;
     infer_parchild_bprop(m[n],theta,mini_L,dtheta,loss,niters*iter_mult);
     my_L += mini_L;

     gradnorm += y[n].rows()*y[n].cols();
     
     // now, propagate dtheta back to W
     for(int c=0; c<m[n].cliques.size(); c++){
       int ctype = m[n].cliquetype(c);
       my_dW[ctype] += dtheta[c]*x[n][c].transpose();
     }
     
     my_dW_serialized = serialize(my_dW);
   }
 }


int main(int argc, char * argv[]){

  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if(world_rank >= 1){

    // parse the command line
    string opt_alg; MatrixXd opt_params; bool have_W; double reg;
    int where_m, where_m_end, where_d, where_d_end, niters, ndata, where_w, where_wout;
    read_command_line(argc, argv, where_m, where_d, ndata, niters, opt_alg, opt_params, have_W, where_w, where_wout, reg, false);
    
    // create structure for weights
    int nctypes; VectorXi output_sizes, input_sizes;
    std::vector<MatrixXd> W = create_W(argc, argv, where_m, where_d, have_W, where_w, nctypes, output_sizes, input_sizes);

    // create a place to store any data and message structures this worker might need
    std::vector<features> x(ndata);
    std::vector<MatrixXi> y(ndata);
    std::vector<Messages> m(ndata);
    std::vector<bool>     t(ndata); // constructed or not?
    for(int i=0; i<ndata; i++)
      t[i] = false;

    int master_rank = 0;
    // repeat as follows:
    // recieve number of data to process
    // if it is -1, quit
    // otherwise, receive a list of data and a parameter gradient, compute the loss, return it
    int n, W_numel, ncalc;
    while(true){
      MPI_Recv(&ncalc, 1, MPI_INT, master_rank, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      if(ncalc==-1){
	MPI_Finalize();
	//return 0;
	exit(0);
      }

      int who_to_calc[ncalc];

      // list of data to process
      MPI_Recv(&who_to_calc, ncalc, MPI_INT, master_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // get the size of W
      MPI_Recv(&W_numel, 1, MPI_INT, master_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      double W_serialized_data[W_numel];

      // recieve it
      MPI_Recv(&W_serialized_data, W_numel, MPI_DOUBLE, master_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // put it into a matrix
      MatrixXd W_serialized(W_numel,1);
      MatrixXd my_dW_serialized = MatrixXd::Zero(W_numel,1); 
      for(int i=0; i<W_numel; i++)
	W_serialized(i) = W_serialized_data[i];

      // get iter mult
      int iter_mult = 0;
      MPI_Recv(&iter_mult, 1, MPI_INT, master_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      double my_L = 0;
      double gradnorm = 0;
      compute_loss_cold(argc,argv,W_serialized,ncalc,who_to_calc,my_L,my_dW_serialized,x,y,m,t,W,nctypes,output_sizes,input_sizes,gradnorm,iter_mult);

      // now send back the loss
      MPI_Send(&my_L, 1, MPI_DOUBLE, master_rank, 0, MPI_COMM_WORLD);

      // and send back the gradient
      MPI_Send(my_dW_serialized.data(), W_numel, MPI_DOUBLE, master_rank, 0, MPI_COMM_WORLD);

      // and send back gradnorm
      MPI_Send(&gradnorm,1,MPI_DOUBLE,master_rank,0,MPI_COMM_WORLD);
    }
  }

    struct timeval start, end;
    long mtime, seconds, useconds;    
    gettimeofday(&start, NULL);

  // parse the command line options.  plan is:
  // 1) find the first and last location of the models
  // 2) find the first and last location of the data

  // get information from the command line
  string opt_alg; MatrixXd opt_params; bool have_W; double reg;
  int where_m, where_m_end, where_d, where_d_end, niters, ndata, where_w, where_wout;
  read_command_line(argc, argv, where_m, where_d, ndata, niters, opt_alg, opt_params, have_W, where_w, where_wout, reg, true);

  // either read W from a file, or create it
  int nctypes; VectorXi output_sizes, input_sizes;
  std::vector<MatrixXd> W = create_W(argc, argv, where_m, where_d, have_W, where_w, nctypes, output_sizes, input_sizes);

  // divide up the data amongst the workers
  int nworkers = min(world_size-1,ndata);
  cout << "nworkers: " << nworkers << endl;
  vector<vector<int>> who_to_calc(nworkers);
  int who = 0;
  for(int n=0; n<ndata; n++){
    who_to_calc[who].push_back(n);
    who = (who + 1) % nworkers;
  }  

  int iter_mult;
  auto erisk = [&](const std::vector<MatrixXd> & W, double & L, std::vector<MatrixXd> & dW){
    MatrixXd myL(ndata,1);
    L  = 0;
    double gradnorm = 0;
    dW = std::vector<MatrixXd>();
    for(int ctype=0; ctype<nctypes; ctype++)
      dW.push_back(MatrixXd::Zero(output_sizes(ctype),input_sizes(ctype)));
    MatrixXd dW_serialized = serialize(dW);
    int numel_dW = dW_serialized.rows();
    
    // assign work to everyone
    for(int worker_rank=0; worker_rank<nworkers; worker_rank++){
      MatrixXd W_serialized     = serialize(W);

      int ncalc = who_to_calc[worker_rank].size();
      // send # data to calculate loss on
      MPI_Send(&ncalc, 1, MPI_INT, worker_rank+1, 0, MPI_COMM_WORLD);
      // send data
      MPI_Send(&who_to_calc[worker_rank][0], ncalc, MPI_INT, worker_rank+1, 0, MPI_COMM_WORLD);
      // send size of W
      MPI_Send(&numel_dW, 1, MPI_INT, worker_rank+1, 0, MPI_COMM_WORLD);
      // send W itself
      MPI_Send(W_serialized.data(), numel_dW, MPI_DOUBLE, worker_rank+1, 0, MPI_COMM_WORLD);
      // send iter_mult
      MPI_Send(&iter_mult, 1, MPI_INT, worker_rank+1, 0, MPI_COMM_WORLD);
    }

    // gather results
    for(int worker_rank=0; worker_rank<nworkers; worker_rank++){
      double my_L = 0;
      MatrixXd my_dW_serialized = MatrixXd(numel_dW,1);
      double my_gradnorm;

      // recieve loss
      MPI_Recv(&my_L, 1, MPI_DOUBLE, worker_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // receive gradient
      MPI_Recv(my_dW_serialized.data(), numel_dW, MPI_DOUBLE, worker_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // receieve gradborm
      MPI_Recv(&my_gradnorm, 1, MPI_DOUBLE, worker_rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      L             += my_L;
      dW_serialized += my_dW_serialized;
      gradnorm      += my_gradnorm;
    }

    L             /= gradnorm;
    dW_serialized /= gradnorm;

    deserialize(dW_serialized,dW);

    //regularization yo yo yo
    double reg_penalty = 0;
    for(int ctype=0; ctype<W.size(); ctype++){
      reg_penalty +=   reg*W[ctype].array().pow(2).sum();
      dW[ctype]   += 2*reg*W[ctype];
    }
    L += reg_penalty;
    cout << "reg penalty: " << reg_penalty << endl;

    gettimeofday(&end, NULL);
    seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    cout << "elapsed time: " << mtime/1000.0 << " seconds" << endl;
  };
  
  /*
  cout << "initial weights:" << endl;
  for(int ctype=0; ctype < nctypes; ctype++)
    cout << " W["<<ctype<<"]:"<<endl<<W[ctype]<<endl;
  */

  if(opt_alg != "lbfgs")
    throw new MyException("optimization algoritm must be lbfgs");

  /*
  iter_mult = 0;
  if(!have_W){
    cout << "since initial weights are zero, initializing with a search with iter=0" << endl;
    W = lbfgs(erisk, W, false);
    }*/
  iter_mult = 1;

  for(int restarts=0; restarts < 4; restarts++){
    cout << "STARTING SEARCH # " << restarts << endl;
    W = lbfgs(erisk, W, false);
  }

  // tell slaves to exit
  int exit_flag = -1;
  for(int worker_rank=0; worker_rank<world_size-1; worker_rank++)
    MPI_Send(&exit_flag, 1, MPI_INT, worker_rank+1, 0, MPI_COMM_WORLD);

  /*
  cout << "weights after optimization:" << endl;
  for(int ctype=0; ctype < nctypes; ctype++)
    cout << " W["<<ctype<<"]:"<<endl<<W[ctype]<<endl;
  */    

  if(where_wout >= argc-1)
    write_params("W.txt",W);
  else
    write_params(argv[where_wout+1],W);
  
    MPI_Finalize();
    exit(0);
    //return 0;
}
