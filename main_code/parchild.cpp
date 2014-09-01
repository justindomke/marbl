// need to implement:
// current_marginal
// marginalize_down

//#define MAXITER 10

void marginalize_down(Messages & m, MatrixXd & mu_alpha, int alpha, int beta, MatrixXd & mu_beta){
  assert(contains(m.cliques[alpha],m.cliques[beta]));
  assert(contains(m.get_pars(beta),alpha));
  assert(contains(m.get_kids(alpha),beta));
  assert(mu_beta.size() == m.nconfigs(beta));
  mu_beta.setZero();
  for(int i=0; i<m.nconfigs(alpha); i++){
    int j         = m.kid_config(alpha,beta,i);
    mu_beta(j)   += mu_alpha(i);
  }
  assert( abs(1 - mu_beta.sum()) < 1e-10);
}

MatrixXd marginalize_down(Messages & m, MatrixXd & mu_alpha, int alpha, int beta){
  assert(contains(m.cliques[alpha],m.cliques[beta]));
  assert(contains(m.get_pars(beta),alpha));
  assert(contains(m.get_kids(alpha),beta));

  MatrixXd mu_beta = MatrixXd(m.nconfigs(beta),1);
  marginalize_down(m, mu_alpha, alpha, beta, mu_beta);
  return mu_beta;
}

void current_marginal(Messages & m, std::vector<MatrixXd> & theta, int alpha, MatrixXd & mu_alpha){
  //MatrixXd mu_alpha = MatrixXd::Zero(m.nconfigs(alpha),1);
  // don't bother if entropy factor is zero
  assert(mu_alpha.size() == m.nconfigs(alpha));
  if(m.ent(alpha)==0){
    mu_alpha.array() *= 0.0;
    mu_alpha.array() += 1.0/(mu_alpha.rows()*mu_alpha.cols());

    MatrixXd tmp(m.nconfigs(alpha),1);
    for(int i=0; i<m.nconfigs(alpha); i++){
      tmp(i) = theta[alpha](i);
      
      for(int n=0; n<m.nkids(alpha); n++){
	int beta = m.get_kids(alpha)[n];
	tmp(i) += m.msg(alpha,beta,i,n);
      }
      for(int n=0; n<m.npars(alpha); n++){
	int gamma = m.get_pars(alpha)[n];
	tmp(i) -= m.par_msg(alpha,gamma,i);
      }
    }
    //cout << "tmp: " << tmp.transpose() << endl;

    return;
  }
  
  for(int i=0; i<m.nconfigs(alpha); i++){
    mu_alpha(i) = theta[alpha](i);
    
    for(int n=0; n<m.nkids(alpha); n++){
      int beta = m.get_kids(alpha)[n];
      mu_alpha(i) += m.msg(alpha,beta,i,n);
      }
    for(int n=0; n<m.npars(alpha); n++){
      int gamma = m.get_pars(alpha)[n];
      mu_alpha(i) -= m.par_msg(alpha,gamma,i);
    }
    mu_alpha(i) /= m.ent(alpha);
  }
  
  mu_alpha.array() -= log_sum_exp(mu_alpha);
  mu_alpha = mu_alpha.array().exp();

  //cout << " mu[" << alpha << "]: " << mu_alpha.transpose() << endl;

  //cout << mu_alpha.transpose() << endl;
  //if(!abs(mu_alpha.sum() - 1) < 1e-4){
  //  cout << " theta: " << theta[alpha].transpose() << endl;
  //  cout << " mu[" << alpha << "]: " << mu_alpha.transpose() << endl;
  //}
  assert(abs(mu_alpha.sum() - 1) < 1e-4);
}

MatrixXd current_marginal(Messages & m, std::vector<MatrixXd> & theta, int alpha){
  MatrixXd mu_alpha = MatrixXd::Zero(m.nconfigs(alpha),1);
  current_marginal(m, theta, alpha, mu_alpha);
  return mu_alpha;
}

void check_marginal_match(Messages & m, std::vector<MatrixXd> & theta, int nu){
  // debugging routine to check that a message update successfully enforced marginalization
  // recompute marginal for nu and parents
  std::vector<MatrixXd> mu_pars;

  MatrixXd mu_nu = current_marginal(m, theta, nu);
  for(int n=0; n<m.get_pars(nu).size(); n++){
    int alpha = m.get_pars(nu)[n];
    mu_pars.push_back(current_marginal(m, theta, alpha));
  }
  
  // marginalize out parents again, check for match
  std::vector<MatrixXd> mu_marg;
  //cout << "mu[" << nu << "]: " << mu_nu.transpose() << " ent: " << m.ent(nu) << endl;
  for(int n=0; n<m.get_pars(nu).size(); n++){
    int alpha = m.get_pars(nu)[n];
    mu_marg.push_back(marginalize_down(m, mu_pars[n], alpha, nu));
    //cout << "rr[" << alpha << "]: " << mu_marg[n].transpose() << endl;
    //assert( (mu_nu.array()-mu_marg[n].array()).abs().maxCoeff()  < 1e-10 );
  }
}

void current_marginal_bprop(Messages & m , std::vector<MatrixXd> & theta , int alpha, MatrixXd & mu_alpha,
			    Messages & dm, std::vector<MatrixXd> & dtheta, MatrixXd & dlogmu_alpha){
  // backprops from dmu_alpha onto dm and dtheta

  // don't bother if entropy factor is zero
  if(m.ent(alpha)==0){
    return;
  }
  
  //mu_alpha = current_marginal(m, theta, alpha, ent);

  //cout << "mu_alpha: " << mu_alpha.transpose() << endl;
  //cout << "mu_alpha: " << current_marginal(m, theta, alpha, ent).transpose() << endl;
  //assert( (mu_alpha-current_marginal(m, theta, alpha, ent)).array().maxCoeff() < 1e-5);  

  double mysum = 0;
  for(int i=0; i<m.nconfigs(alpha); i++)
    mysum += dlogmu_alpha(i);
  
  for(int i=0; i<m.nconfigs(alpha); i++){
    //mu_alpha(i) = theta[alpha](i);
    dtheta[alpha](i) += (1/m.ent(alpha))*(dlogmu_alpha(i) - mysum*mu_alpha(i));
    //dtheta[alpha](i) += (1/ent(alpha))*(dlogmu_alpha(i) );
    
    for(int n=0; n<m.nkids(alpha); n++){
      int beta = m.get_kids(alpha)[n];
      //mu_alpha(i) += m.msg(alpha,beta,i,n);
      dm.msg(alpha,beta,i,n) += (1/m.ent(alpha))*(dlogmu_alpha(i) - mysum*mu_alpha(i));
      //dm.msg(alpha,beta,i,n) += (1/ent(alpha))*(dlogmu_alpha(i) );
      }
    for(int n=0; n<m.npars(alpha); n++){
      int gamma = m.get_pars(alpha)[n];
      //mu_alpha(i) -= m.par_msg(alpha,gamma,i);
      dm.par_msg(alpha,gamma,i) -= (1/m.ent(alpha))*(dlogmu_alpha(i) - mysum*mu_alpha(i));
      //dm.par_msg(alpha,gamma,i) -= (1/ent(alpha))*(dlogmu_alpha(i) );
    }
  }
}

template <class T>
void infer_parchild_bprop(Messages & m, std::vector<MatrixXd> & theta, double & L, std::vector<MatrixXd> & dtheta, T & loss, int niters ){
  // make a set of uniform marginals

  std::vector<MatrixXd> mu;
  std::vector<MatrixXd> dlogmu;
  for(int alpha=0; alpha<m.cliques.size(); alpha++){
    mu.push_back(MatrixXd::Zero(m.nconfigs(alpha),1));
    dlogmu.push_back(MatrixXd::Zero(m.nconfigs(alpha),1));
  }

  Messages dm = Messages(m.cliques, m.nvals, m.nnodes, m.ent, m.cliquetype);

  //std::vector<MatrixXd> dtheta;
  assert(dtheta.size() == 0);
  for(int alpha=0; alpha<m.cliques.size(); alpha++){
    //cout << "mysize: " << theta[alpha].size() << endl;
    dtheta.push_back(MatrixXd::Zero(theta[alpha].size(),1));
  }

  // initialize messages using thetas
  for(int alpha=0; alpha<m.cliques.size(); alpha++){
    for(int n=0; n<m.get_pars(alpha).size(); n++){    
      int gamma = m.get_pars(alpha)[n];
      for(int x_alpha=0; x_alpha <  m.nconfigs(alpha); x_alpha++){
	//cout << "old: " << m.par_msg(alpha,gamma,x_alpha);
	m.par_msg(alpha,gamma,x_alpha) = m.ent(gamma)*theta[alpha](x_alpha);
	//cout << "  new: " << m.par_msg(alpha,gamma,x_alpha) << endl;
      }
    }
  }

  std::stack<double> mystack;

  /*
  // put messages for zero entropy nodes
  for(int alpha=0; alpha<m.cliques.size(); alpha++)
    if(ent(alpha)==0)
      for(int n=0; n<m.npars(alpha); n++){
	int gamma = m.get_pars(alpha)[n];
	for(int i=0; i<m.nconfigs(alpha); i++)
	  m.par_msg(alpha,gamma,i) = theta[alpha](i) / m.npars(alpha);
      }
  */

  //vector<MatrixXd> mu_margAR[m.cliques.size()];
  auto mu_margAR = new vector<MatrixXd>[m.cliques.size()];
  //vector<MatrixXd> dlogmu_margAR[m.cliques.size()];
  auto dlogmu_margAR = new vector<MatrixXd>[m.cliques.size()];
  for(int alpha=0; alpha<m.cliques.size(); alpha++){
    for(int n=0; n<m.get_pars(alpha).size(); n++){
      mu_margAR[alpha].push_back(MatrixXd(m.nconfigs(alpha),1));
      dlogmu_margAR[alpha].push_back(MatrixXd::Zero(m.nconfigs(alpha),1));
    }
  }

  for(int iter=0; iter<niters; iter++){
    //for(int nu=0; nu<m.cliques.size(); nu++){
    
    for(int nu0=0; nu0<m.cliques.size()*2; nu0++){
      int nu = nu0;
      if(nu0 >= m.cliques.size())
	nu = 2*m.cliques.size()-1-nu0;
  
      // compute marginal for nu and parents
      current_marginal(m, theta, nu, mu[nu]);
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	current_marginal(m, theta, alpha, mu[alpha]);
      }
      
      // marginalize out parents
      //std::vector<MatrixXd> mu_marg;
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	//mu_marg.push_back(marginalize_down(m, mu[alpha], alpha, nu));
	//mu_margAR[nu][n] = marginalize_down(m, mu[alpha], alpha, nu);
	marginalize_down(m, mu[alpha], alpha, nu, mu_margAR[nu][n]);
      }

      double n_nu = 0;
      for(int n=0; n<m.get_pars(nu).size(); n++)
	n_nu += m.ent(m.get_pars(nu)[n]);

      // update the stupid message      
      for(int n=0; n<m.get_pars(nu).size(); n++){

	int alpha = m.get_pars(nu)[n];
	
	for(int x_nu=0; x_nu < m.nconfigs(nu); x_nu++){
	  mystack.push(m.par_msg(nu,alpha,x_nu));
	  m.par_msg(nu,alpha,x_nu) += (m.ent(alpha)/(m.ent(nu)+n_nu))*(m.ent(nu) * log( mu[nu](x_nu) ));	  

	  for(int np=0; np < m.get_pars(nu).size(); np++){
	    int alphap = m.get_pars(nu)[np];
	    m.par_msg(nu,alpha,x_nu) += (m.ent(alpha)/(m.ent(nu)+n_nu))*(m.ent(alphap)* log( mu_margAR[nu][np](x_nu)));
	  }
	  m.par_msg(nu,alpha,x_nu) -= m.ent(alpha)*log(mu_margAR[nu][n](x_nu));
	}

	m.msg(alpha,nu).array() -= m.msg(alpha,nu).array().mean();
      }
    }
  }

  // compute mu for all non-zero entropy nodes
  for(int alpha=0; alpha<m.cliques.size(); alpha++)
    if(m.ent(alpha) != 0)
      mu[alpha] = current_marginal(m, theta, alpha);
  
  // now, do it for zero-entropy nodes
  for(int nu=0; nu<m.cliques.size(); nu++)
    if(m.ent(nu) == 0){
      bool updated = false;
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	if(m.ent(alpha)!=0){
	  mu[nu] = marginalize_down(m, mu[alpha], alpha, nu);
	  updated = true;
	  break;
	}
      }
      assert(updated);
    }

  // compute loss
  // call loss function that was passed in
  loss(mu,L,dlogmu);

  // push the losses for zero-entropy nodes onto the other nodes they came from
  for(int nu=0; nu<m.cliques.size(); nu++){
    if(m.ent(nu) == 0){
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	if(m.ent(alpha)!=0){
	  for(int x_alpha=0; x_alpha<m.nconfigs(alpha); x_alpha++){
	    int x_nu = m.kid_config(alpha, nu, x_alpha);
	    dlogmu[alpha](x_alpha) += dlogmu[nu](x_nu) * mu[alpha](x_alpha) / mu[nu](x_nu);
	  }
	  for(int x_nu=0; x_nu<m.nconfigs(nu); x_nu++)
	    dlogmu[nu](x_nu) = 0;
	  break;
	}
      }
    }
  }

  // backprop marginals onto messages
  for(int alpha=0; alpha<m.cliques.size(); alpha++)
    if(m.ent(alpha)!=0)
      current_marginal_bprop(m, theta, alpha, mu[alpha], dm, dtheta, dlogmu[alpha]);  

  // how backprop works:
  //
  // we are always maintaining two things:
  // dtheta
  // dm
  //
  // go through nodes in reverse order
  // for each node:
  // initialize dlogmu to zero
  // backprop dm onto dlogmu, dlogmu_MARG
  // backprop dlogmu_MARG into dlogmu
  // backprop logmu onto dm and dtheta
  // roll back mu to previous values

  for(int iter=0; iter<niters; iter++){
    //for(int nu0=0; nu0<m.cliques.size()*2; nu0++){
    for(int nu0=m.cliques.size()*2-1; nu0>=0; nu0--){
      int nu = nu0;
      if(nu0 >= m.cliques.size())
	nu = 2*m.cliques.size()-1-nu0;
  
      double n_nu = 0;
      for(int n=0; n<m.get_pars(nu).size(); n++)
	n_nu += m.ent(m.get_pars(nu)[n]);

      dlogmu[nu].setZero();
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	dlogmu[alpha].setZero();
      }
      for(int n=0; n<m.get_pars(nu).size(); n++)
	dlogmu_margAR[nu][n].setZero();

      // update the stupid message      
      for(int n=m.get_pars(nu).size()-1; n>=0; n--){
	int alpha = m.get_pars(nu)[n];

	for(int x_nu = m.nconfigs(nu)-1; x_nu>=0; x_nu--){
          //m.par_msg(nu,alpha,x_nu) += (ent(alpha)/(ent(nu)+n_nu))*(ent(nu) * log( mu[nu](x_nu) ));	  
	  dlogmu[nu](x_nu) += dm.par_msg(nu,alpha,x_nu) * (m.ent(alpha)/(m.ent(nu)+n_nu))*m.ent(nu);

	  for(int np=0; np < m.get_pars(nu).size(); np++){
	    int alphap = m.get_pars(nu)[np];
	    //m.par_msg(nu,alpha,x_nu) += (ent(alpha)/(ent(nu)+n_nu))*(ent(alphap) * log( mu_margAR[nu][np](x_nu)));
	    dlogmu_margAR[nu][np](x_nu) += dm.par_msg(nu,alpha,x_nu) * (m.ent(alpha)/(m.ent(nu)+n_nu))*m.ent(alphap);
	  }
	  //m.par_msg(nu,alpha,x_nu) -= ent(alpha)*log(mu_margAR[nu][n](x_nu));
	  dlogmu_margAR[nu][n](x_nu) -= m.ent(alpha)*dm.par_msg(nu,alpha,x_nu);

	  m.par_msg (nu, alpha, x_nu) = mystack.top(); mystack.pop();
	  // this line cost me a couple days!
	  //dm.par_msg(nu, alpha, x_nu) = 0;
	}
      }      

      // compute marginal for nu and parents
      current_marginal(m, theta, nu, mu[nu]);      
      //dlogmu[nu].setZero();
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	current_marginal(m, theta, alpha, mu[alpha]);
	//dlogmu[alpha].setZero();
      }
      
      // marginalize out parents
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	marginalize_down(m, mu[alpha], alpha, nu, mu_margAR[nu][n]);
      }

      // backprop from marginalized parent marginals to parent marginals themselves
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	for(int x_alpha=0; x_alpha<m.nconfigs(alpha); x_alpha++){
	  int x_nu = m.kid_config(alpha, nu, x_alpha);
	  dlogmu[alpha](x_alpha) += dlogmu_margAR[nu][n](x_nu)*mu[alpha](x_alpha)/mu_margAR[nu][n](x_nu);
	}
      }

      current_marginal_bprop(m, theta, nu, mu[nu], dm, dtheta, dlogmu[nu]);
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	current_marginal_bprop(m, theta, alpha, mu[alpha], dm, dtheta, dlogmu[alpha]);
      }
    }
  }

  // re-propagate for initial messages
  for(int alpha=0; alpha<m.cliques.size(); alpha++){
    for(int n=0; n<m.get_pars(alpha).size(); n++){    
      int gamma = m.get_pars(alpha)[n];
      for(int x_alpha=0; x_alpha <  m.nconfigs(alpha); x_alpha++){
	dtheta[alpha](x_alpha) += dm.par_msg(alpha,gamma,x_alpha)*m.ent(gamma);
      }
    }
  }

  delete [] mu_margAR;
  delete [] dlogmu_margAR;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ideas to speed this up:
// 1: create a flag that tells us when marginals need to be updated
// 2: (better) greedily update the marginals when messages change

std::vector<MatrixXd> infer_parchild(Messages & m, std::vector<MatrixXd> & theta, MatrixXd & ent){
  // make a set of uniform marginals
  std::vector<MatrixXd> mu;
  for(int alpha=0; alpha<m.cliques.size(); alpha++){
    mu.push_back(MatrixXd::Zero(m.nconfigs(alpha),1));
  }

  // initialize messages using thetas
  for(int alpha=0; alpha<m.cliques.size(); alpha++){
    for(int n=0; n<m.get_pars(alpha).size(); n++){    
      int gamma = m.get_pars(alpha)[n];
      for(int x_alpha=0; x_alpha <  m.nconfigs(alpha); x_alpha++){
	//cout << "old: " << m.par_msg(alpha,gamma,x_alpha);
	m.par_msg(alpha,gamma,x_alpha) = m.ent(gamma)*theta[alpha](x_alpha);
	//cout << "  new: " << m.par_msg(alpha,gamma,x_alpha) << endl;
      }
    }
  }

  /*
  // put messages for zero entropy nodes
  for(int alpha=0; alpha<m.cliques.size(); alpha++)
    if(ent(alpha)==0)
      for(int n=0; n<m.npars(alpha); n++){
	int gamma = m.get_pars(alpha)[n];
	for(int i=0; i<m.nconfigs(alpha); i++)
	  m.par_msg(alpha,gamma,i) = theta[alpha](i) / m.npars(alpha);
      }
  */

  vector<vector<MatrixXd>> mu_margAR(m.cliques.size());
  for(int alpha=0; alpha<m.cliques.size(); alpha++){
    for(int n=0; n<m.get_pars(alpha).size(); n++)
      mu_margAR[alpha].push_back(MatrixXd(m.nconfigs(alpha),1));
  }

  int MAXITER = 5;
  for(int iter=0; iter<MAXITER; iter++){
    //for(int nu=0; nu<m.cliques.size(); nu++){
    
    for(int nu0=0; nu0<m.cliques.size()*2; nu0++){
      int nu = nu0;
      if(nu0 >= m.cliques.size())
	nu = 2*m.cliques.size()-1-nu0;
  
      // compute marginal for nu and parents
      current_marginal(m, theta, nu, mu[nu]);      
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	current_marginal(m, theta, alpha, mu[alpha]);
      }
      
      // marginalize out parents
      //std::vector<MatrixXd> mu_marg;
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	//mu_marg.push_back(marginalize_down(m, mu[alpha], alpha, nu));
	mu_margAR[nu][n] = marginalize_down(m, mu[alpha], alpha, nu);
	//marginalize_down(m, mu[alpha], alpha, nu, mu_margAR[nu][n]);
      }

      double n_nu = 0;
      for(int n=0; n<m.get_pars(nu).size(); n++)
	n_nu += ent(m.get_pars(nu)[n]);

      // update the stupid message      
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	
	for(int x_nu=0; x_nu < m.nconfigs(nu); x_nu++){
	  m.par_msg(nu,alpha,x_nu) += (ent(alpha)/(ent(nu)+n_nu))*(ent(nu) * log( mu[nu](x_nu) ));	  
	  for(int np=0; np < m.get_pars(nu).size(); np++){
	    int alphap = m.get_pars(nu)[np];
	    //m.par_msg(nu,alpha,x_nu) += (ent(alpha)/(ent(nu)+n_nu))*(ent(alphap)* log( mu_marg[np](x_nu)));
	    m.par_msg(nu,alpha,x_nu) += (ent(alpha)/(ent(nu)+n_nu))*(ent(alphap)* log( mu_margAR[nu][np](x_nu)));
	  }
	  //m.par_msg(nu,alpha,x_nu) -= ent(alpha)*log(mu_marg[n](x_nu));
	  m.par_msg(nu,alpha,x_nu) -= ent(alpha)*log(mu_margAR[nu][n](x_nu));
	}
	m.msg(alpha,nu).array() -= m.msg(alpha,nu).array().mean();
      }

      #ifdef DEBUG
      check_marginal_match(m, theta, nu);
      #endif
    }
  }
  
  // compute mu for all non-zero entropy nodes
  for(int alpha=0; alpha<m.cliques.size(); alpha++)
    if(ent(alpha) != 0)
      mu[alpha] = current_marginal(m, theta, alpha);
  
  // now, do it for zero-entropy nodes
  for(int nu=0; nu<m.cliques.size(); nu++)
    if(ent(nu) == 0){
      bool updated = false;
      for(int n=0; n<m.get_pars(nu).size(); n++){
	int alpha = m.get_pars(nu)[n];
	if(ent(alpha)!=0){
	  mu[nu] = marginalize_down(m, mu[alpha], alpha, nu);
	  updated = true;
	  break;
	}
      }
      assert(updated);
    }

  return mu;
}

