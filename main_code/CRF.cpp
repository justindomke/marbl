
std::vector<MatrixXd> W2theta(const std::vector<MatrixXd> & W, const Messages & m, const std::vector<MatrixXd> & x){
  std::vector<MatrixXd> theta;
  assert(m.cliques.size() == x.size());
  for(int c=0; c<m.cliques.size(); c++){
      int ctype = m.cliquetype(c);
      //cout << "c=" << c << endl;
      //cout << "W[ctype].size: " << W[ctype].rows() << " x " << W[ctype].cols() << endl;
      //cout << "x[c].size: " << x[c].rows() << " x " << x[c].cols() << endl;
      theta.push_back(W[ctype]*x[c]);
      /*
      if(theta[c].array().abs().maxCoeff() > 1e-5){
	cout << "theta["<<c<<"]: " << theta[c].transpose() << endl;
	cout << "W: " << W[ctype] << endl;
	cout << "x[c]: " << x[c] << endl;
	}*/
  }
  return theta;
}

std::vector<MatrixXd> fit_CRF(const std::vector<std::vector<MatrixXd>> & x, const std::vector<MatrixXi> y, std::vector<Messages> & m, string opt_alg, MatrixXd opt_params, int niters, std::vector<MatrixXd> W){

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

  Eigen::initParallel();

  // now, create an evaluation function
  // n says who to evaluate
  auto erisk = [&](const std::vector<MatrixXd> & W, double & L, std::vector<MatrixXd> & dW, MatrixXi who2eval){
      L  = 0;
      dW = std::vector<MatrixXd>();

      for(int ctype=0; ctype < nctypes; ctype++)
          dW.push_back(MatrixXd::Zero(output_sizes(ctype),input_sizes(ctype)));
      
      double gradnorm = 0;
      for(int n=0; n<x.size(); n++)
	gradnorm += y[n].size();

      #pragma omp parallel for
      //for(int n=0; n<x.size(); n++){
      for(int where=0; where < who2eval.size(); where++){
	int n = who2eval(where);

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
	auto loss = [&n, y_configs, gradnorm](const std::vector<MatrixXd> & mu, double & L, std::vector<MatrixXd> & dlogmu){
	    L = 0;
	    for(int alpha=0; alpha<mu.size(); alpha++){
	      if(y_configs(alpha)==-1) continue; // skip when y_configs == -1

	      L -= log(mu[alpha](y_configs(alpha)))/gradnorm;
	      dlogmu[alpha].setZero();
	      dlogmu[alpha](y_configs(alpha)) = -1/gradnorm;

	      if(L!=L) cout << "badmu: " << mu[alpha].transpose() << endl;
	    }
          };
          
          double myL;
          std::vector<MatrixXd> dtheta;
	  
          infer_parchild_bprop(m[n],theta,myL,dtheta,loss,niters);

	  // now, propagate dtheta back to W
	  #pragma omp critical
	  {	    
	    L += myL;
	    for(int c=0; c<m[n].cliques.size(); c++){
              int ctype = m[n].cliquetype(c);
              dW[ctype] += dtheta[c]*x[n][c].transpose();
	    }
	  }
      }

      //regularization yo yo yo
      double reg = .01;
      for(int ctype=0; ctype<W.size(); ctype++){
	L         +=   reg*W[ctype].array().pow(2).sum();
	dW[ctype] += 2*reg*W[ctype];
      }
  };
  
  auto erisk_all = [&](const std::vector<MatrixXd> & W, double & L, std::vector<MatrixXd> & dW){
    MatrixXi who2eval(x.size(),1);
    for(int n=0; n<x.size(); n++)
      who2eval(n) = n;
    erisk(W, L, dW, who2eval);
  };

  auto erisk_1 = [&](const std::vector<MatrixXd> & W, double & L, std::vector<MatrixXd> & dW, int i){
    MatrixXi who2eval(1,1);
    who2eval(0) = i;
    erisk(W, L, dW, who2eval);
  };

  cout << "initial weights:" << endl;
  for(int ctype=0; ctype < nctypes; ctype++)
    cout << " W["<<ctype<<"]:"<<endl<<W[ctype]<<endl;
  
  if(opt_alg == "lbfgs"){
    cout << "calling lbfgs:" << endl;
    W = lbfgs(erisk_all, W, false);
  } else if(opt_alg == "sgd"){
    int opt_iters    = opt_params(0);
    double opt_stepsize = opt_params(1);    
    W = sgd(erisk_1, W, x.size(), opt_iters, opt_stepsize);
  } else{
    throw new MyException("unknown optimization algorithm: " + opt_alg);
  }

  cout << "weights after optimization:" << endl;
  for(int ctype=0; ctype < nctypes; ctype++)
    cout << " W["<<ctype<<"]:"<<endl<<W[ctype]<<endl;

  return W;
}

std::vector<MatrixXd> fit_CRF(const std::vector<std::vector<MatrixXd>> & x, const std::vector<MatrixXi> y, std::vector<Messages> & m, string opt_alg, MatrixXd opt_params, int niters){

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
    W.push_back(MatrixXd::Zero(output_sizes(ctype),input_sizes(ctype)));
  }

  return fit_CRF(x, y, m, opt_alg, opt_params, niters, W);

}
