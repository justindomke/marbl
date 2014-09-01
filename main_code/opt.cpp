//#include <iostream>
//#include "eigenstuff.cpp"
#include "lbfgs.h"
//#include <numeric_limits>
#include <random>

template <class T>
MatrixXd gradient_descent(T fun, MatrixXd x){
  double f;
  MatrixXd g;
  for(int i=0; i<100; i++){
    fun(x,f,g);
    if(i % 10 == 0)
      cout << "iter: " << i << " f: " << f << " x: " << x.transpose() << endl;
    x -= .1*g;
  }

  return x;
}

template <class T>
void check_grad(const T & fun, MatrixXd x){
  int N = x.size();
  
  MatrixXd g0 = MatrixXd(N,1);
  double L;
  fun(x,L,g0);

  cout << "g0: " << g0.transpose() << endl;

  // space for evaluating functions
  MatrixXd tmp = MatrixXd(N,1);
  double L1, L2;

  //double e = 1e-5;
  // fancy automatic size selection
  double e = pow(std::numeric_limits<double>::epsilon(),1.0/3)*(1+x.array().abs().maxCoeff())/100;

  cout << "using e = " << e << endl;

  if(N <= 10){
    cout << "because vector is low dimensional, testing all directions... ";
    MatrixXd g2 = MatrixXd(N,1);
    for(int i=0; i<N; i++){
      MatrixXd x1 = x;
      x1(i) -= e;
      fun(x1,L1,tmp);

      MatrixXd x2 = x;
      x2(i) += e;
      fun(x2,L2,tmp);

      g2(i) = (L2-L1)/(2*e);
    }

    double err = (g0.array()-g2.array()).abs().maxCoeff();
    cout << "max error in gradient: " << err << endl;
  } else{
    cout << "because vector is high dimensional, using 5 random directions... ";
    int M = 5;
    MatrixXd g1 = MatrixXd(M,1);
    MatrixXd g2 = MatrixXd(M,1);
    for(int i=0; i<M; i++){
      MatrixXd r = MatrixXd::Random(N,1);

      g1(i) = (g0.array()*r.array()).sum();

      MatrixXd x1 = x;
      x1 -= e*r;
      fun(x1,L1,tmp);

      MatrixXd x2 = x;
      x2 += e*r;
      fun(x2,L2,tmp);

      g2(i) = (L2-L1)/(2*e);
    }

    double err = (g1.array()-g2.array()).abs().maxCoeff();
    cout << "max error in gradient: " << err << endl;
    cout << "g1: " << g1.transpose() << endl;
    cout << "g2: " << g2.transpose() << endl;
    assert(err < 1e-3);
  }
  return;
}

template <class T>
MatrixXd lbfgs(T & fun, MatrixXd x0, bool do_grad_check){
  if(do_grad_check)
    check_grad(fun,x0);

  int N = x0.size();

  int i, ret = 0;
  lbfgsfloatval_t fx;
  lbfgsfloatval_t *x = lbfgs_malloc(N);
  lbfgs_parameter_t param;

  for(int i=0; i<N; i++){
    x[i] = x0(i);
  }

  auto evaluate = [](void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step){
    int i;
    lbfgsfloatval_t fx = 0.0;

    // extract input vector
    MatrixXd my_x(n,1);
    for(int i=0; i<n; i++)
      my_x(i) = x[i];

    //cout << "my_x: " << my_x.transpose() << endl;

    // call function
    // need to use this void pointer horror, since for some reason having the lambda just capture fun causes tons of trouble
    T fun = *static_cast<T*>(instance);
    double my_f = 0;
    MatrixXd my_g = MatrixXd::Zero(n,1);
    fun(my_x, my_f, my_g);

    // guard function to deal with NaN badness
    bool hasNaN = (my_f != my_f);
    for(int j=0; j<my_g.size(); j++)
      if(my_g(j)!=my_g(j))
	hasNaN = true;
    if(hasNaN){
      my_f = 1e9;
      my_g = MatrixXd::Zero(n,1);	
    }

    // store in lbfgs format
    fx = my_f;
    for(int i=0; i<n; i++){
      g[i] = my_g(i);
    }

    return fx;
  };

  auto progress = [](void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls){
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    /*
    cout << "x: ";
    for(int i=0; i<n; i++)
      cout << x[i] << " ";
    cout << endl;
    printf("\n");
    */
    return 0;
  };

  lbfgs_parameter_init(&param);

  //param.epsilon = 1e-20;
  //param.delta   = 1e-20;
  //param.min_step = .05;
  param.m        = 60;
  //param.ftol    = 1e-6;
  //param.gtol    = .01;
  param.gtol = .99;
  //param.xtol    = 1e-10;
  param.max_iterations = 250;
  //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO; // harms perf.

  ret = lbfgs(N, x, &fx, evaluate, progress, &fun, &param);

  /* Report the result. */
  printf("L-BFGS optimization terminated with status code = %d\n", ret);
  printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
  
  MatrixXd x_out(N,1);
  for(int i=0; i<N; i++)
    x_out(i) = x[i];

  lbfgs_free(x);
  return x_out;
}

template <class T>
MatrixXd lbfgs(T & fun, MatrixXd x0){
  return lbfgs(fun,x0,false);
}


template <class T>
std::vector<MatrixXd> lbfgs(const T &fun, std::vector<MatrixXd> x, bool do_grad_check){
  int N = 0;
  for(int i=0; i<x.size(); i++){
    N += x[i].size();
  }

  std::vector<MatrixXd> g;
  for(int i=0; i<x.size(); i++){
    g.push_back(MatrixXd(x[i].size(),1));
  }

  MatrixXd params = MatrixXd(N,1);
  int where = 0;
  for(int i=0; i<x.size(); i++){
    params.block(where,0,x[i].size(),1) = Map<MatrixXd>(x[i].data(), x[i].size(), 1);
    where += x[i].size(); // this was missing, and things seemed to work fine, disturbingly...
  }

  auto fun2 = [&](const MatrixXd & params, double & f, MatrixXd & grad){
    // first, unpack params into x0
    int where2 = 0;
    for(int i=0; i<x.size(); i++){
        for(int j=0; j<x[i].size(); j++)
            x[i](j) = params(where2+j);
        where2 += x[i].size();
    }

    // now, call original function
    fun(x,f,g);
    
    // now, put g into grad
    where2 = 0;
    for(int i=0; i<x.size(); i++){
        for(int j=0; j<x[i].size(); j++)
            grad(where2+j) = g[i](j);
      where2 += x[i].size();
    }    
  };

  params = lbfgs(fun2, params, do_grad_check);

  where = 0;
  for(int i=0; i<x.size(); i++){
    Map<MatrixXd>(x[i].data(), x[i].size(), 1) = params.block(where,0,x[i].size(),1);
    where += x[i].size();
  }
  
  return x;
}

template <class T>
std::vector<MatrixXd> lbfgs(const T &fun, std::vector<MatrixXd> x){
  return lbfgs(fun, x, false);
}


bool good_iter_to_print(int i){
  if(i==0) return true;
  int a = pow(10,floor(log10(i*1.0)));
  return (i%a)==0;
}



template <class T>
MatrixXd sgd(T & fun, MatrixXd x0, int ndata, int niter, double step){
  // could do this using a lambda
  //if(do_grad_check)
  //  check_grad(fun,x0);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, ndata-1);

  int N      = x0.size();
  MatrixXd x = x0;
  double f   = 0;
  MatrixXd g = MatrixXd::Zero(N,1);

  cout << "niter: " << niter << " step: " << step << endl;

  //int niter   = 10000;
  //double step = .025;
  //step        = .001;

  double f_smooth = 0;

  MatrixXd r = MatrixXd::Zero(N,1).array()+.0001;
  MatrixXd p = MatrixXd::Zero(N,1);

  for(int iter=0; iter < niter; iter++){
    int i = dis(gen); // random between 0 and ndata

    cout << "r: " << r.transpose() << endl;

    fun(x, f, g, i);
    //x -= step*g;

    r += g.array().pow(2).matrix();
    //MatrixXd g_scaled = (g.array()/r.array().sqrt()).matrix();
    MatrixXd g_scaled = g;
    double lambda = .1;
    p = (1-lambda)*p + lambda*g_scaled;
    x -= step*p;

    double alpha = max(.1,1.0/(iter+1));
    f_smooth = (1-alpha)*f_smooth + alpha*f;

    if(1 || good_iter_to_print(iter))
      cout << "iter: " << iter << "  f_smooth: " << f_smooth*ndata << endl;
  }

  return x;
}


template <class T>
std::vector<MatrixXd> sgd(const T &fun, std::vector<MatrixXd> x, int ndata, int niter, double step){
  int N = 0;
  for(int i=0; i<x.size(); i++){
    N += x[i].size();
  }

  std::vector<MatrixXd> g;
  for(int i=0; i<x.size(); i++){
    g.push_back(MatrixXd(x[i].size(),1));
  }

  MatrixXd params = MatrixXd(N,1);
  int where = 0;
  for(int i=0; i<x.size(); i++){
    params.block(where,0,x[i].size(),1) = Map<MatrixXd>(x[i].data(), x[i].size(), 1);
    where += x[i].size(); // this was missing, and things seemed to work fine, disturbingly...
  }

  auto fun2 = [&](const MatrixXd & params, double & f, MatrixXd & grad, int i){
    // first, unpack params into x0
    int where2 = 0;
    for(int i=0; i<x.size(); i++){
      cout << "x[" << i << "]: " << endl << x[i] << endl;
      for(int j=0; j<x[i].size(); j++)
	x[i](j) = params(where2+j);
      where2 += x[i].size();
    }

    // now, call original function
    fun(x,f,g,i);
    
    // now, put g into grad
    where2 = 0;
    for(int i=0; i<x.size(); i++){
        for(int j=0; j<x[i].size(); j++)
            grad(where2+j) = g[i](j);
      where2 += x[i].size();
    }    
  };

  params = sgd(fun2, params, ndata, niter, step);

  where = 0;
  for(int i=0; i<x.size(); i++){
    Map<MatrixXd>(x[i].data(), x[i].size(), 1) = params.block(where,0,x[i].size(),1);
    where += x[i].size();
  }
  
  return x;
}


void test_sgd(){
  int D      = 4;
  int ndata  = 150;
  MatrixXd r = -1 + 2*MatrixXd::Random(D,ndata).array();
  MatrixXd v = -1 + 2*MatrixXd::Random(D,1).array();
  auto fun = [&](const MatrixXd & x, double & f, MatrixXd & g, int i){
    f  = exp((x.array()*r.col(i).array()).sum());
    g  = exp((x.array()*r.col(i).array()).sum())*r.col(i);

    //f += (x.array()-v.array()).pow(2).sum();
    //g += 2*(x-v);
    return;
  };

  auto fun_all = [&](const MatrixXd & x, double & f, MatrixXd & g){
    double my_f   = 0;
    MatrixXd my_g = MatrixXd::Zero(g.rows()*g.cols(),1);
      cout << "g: " << g.transpose() << " my_g: " << my_g.transpose() << endl;
    for(int i=0; i<ndata; i++){ 
      fun(x,my_f,my_g,i);
      f += my_f;
      g += my_g;
    }
    return;
  };

  MatrixXd x = MatrixXd::Random(D,1);

  check_grad(fun_all, x);

  //auto xp = lbfgs(fun, x);
  //gradient_descent(fun, x);
  auto xp  = sgd(fun,x,ndata, 100, .1);
  auto xpp = lbfgs(fun_all,x);

  double f, fp, fpp;
  MatrixXd g(D,1);

  cout << "x  : " << x.transpose()  << endl;
  cout << "xp : " << xp.transpose() << endl;
  cout << "xpp: " << xpp.transpose() << endl;
  
  fun_all(x,  f,  g);
  fun_all(xp ,fp ,g);
  fun_all(xpp,fpp,g);

  cout << "f(x): " << f << " f(xp): " << fp << " f(xpp): " << fpp << endl;

}


void test_optimization(){
  int D = 20;
  MatrixXd r = MatrixXd::Random(D,1);
  MatrixXd v = MatrixXd::Random(D,1);
  auto fun = [&](const MatrixXd & x, double & f, MatrixXd & g){
    //f = x.array().pow(2).sum();
    //g = 2*x;
    f  = exp((x.array()*r.array()).sum());
    g  = exp((x.array()*r.array()).sum())*r;
    //f += x.array().abs().sum();
    //g += x.array()/x.array().abs();
    f += (x.array()-v.array()).pow(2).sum();
    g += 2*(x-v);
    return;
  };
  MatrixXd x = MatrixXd::Random(D,1);

  check_grad(fun, x);

  auto xp = lbfgs(fun, x);
  //gradient_descent(fun, x);

  double f;
  double fp;
  MatrixXd g;
  fun(x, f, g);
  fun(xp,fp,g);

  cout << "f(x): " << f << " f(xp): " << fp << endl;
}


void test_optimization2(){
  int D = 10;
 
  std::vector<MatrixXd> r, v, x;
  for(int i=0; i<100; i++){
    r.push_back(MatrixXd::Random(D,1));
    v.push_back(MatrixXd::Random(D,1));
    x.push_back(MatrixXd::Random(D,1));
  }

  auto fun = [&](const std::vector<MatrixXd> & x, double & f, std::vector<MatrixXd> & g){
    f = 0;
    for(int i=0; i<x.size(); i++){
      f    += exp((x[i].array()*r[i].array()).sum());
      g[i]  = exp((x[i].array()*r[i].array()).sum())*r[i];
      f    += (x[i].array()-v[i].array()).pow(2).sum();
      g[i] += 2*(x[i]-v[i]);
    }
    return;
  };

  auto xp = lbfgs(fun, x, true);

  /*
  double f;
  double fp;
  MatrixXd g;
  fun(x, f, g);
  fun(xp,fp,g);*/

  //cout << "f(x): " << f << " f(xp): " << fp << endl;
}


/*
int main(){
  //test_optimization2();
  test_sgd();

  return 0;
  }*/
