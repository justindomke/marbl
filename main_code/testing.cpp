void test_config_tforms(Messages & m){
  for(int alpha=0; alpha<m.cliques.size(); alpha++){
    for(int n=0; n<m.get_kids(alpha).size(); n++){
      int beta  = m.get_kids(alpha)[n];
      //cout << "alpha: " << alpha << "  beta: " << beta << endl;

      for(int i=0; i<m.nconfigs(alpha); i++){
	//cout << "i: " << i << " vec: " << m.config2vec(alpha,i).transpose() << " and back: " << m.vec2config(alpha,m.config2vec(alpha,i)) << endl;
	MatrixXi vec = m.config2vec(alpha,i);
	int config = m.vec2config(alpha,vec);
	assert(i == config);
      }
    }
  }
}

void test_indexing(Messages & m){
  // find a node that has some kids
  for(int alpha=0; alpha < m.cliques.size(); alpha++){
    for(int n=0; n<m.get_kids(alpha).size(); n++){
      int beta = m.get_kids(alpha)[n];
      
      cout << "alpha cliques: " << m.cliques[alpha].transpose() << endl;
      cout << "beta  cliques: " << m.cliques[beta].transpose() << endl;   
      
      for(int i=0; i<m.nconfigs(alpha); i++){
	int j = m.kid_config(alpha,beta,i);
	MatrixXi x = m.config2vec(alpha,i);
	MatrixXi y = m.config2vec(beta ,j);
	cout << "i: " << i << " j: " << j << " x: " << x.transpose() << " y: " << y.transpose() << endl;
      }
    }
  }
}

/*
void test_msg_indexing(Messages & m){
  for(int alpha=0; alpha<m.cliques.size(); alpha++){
    for(int n=0; n<m.get_kids(alpha).size(); n++){
      int beta = m.get_kids(alpha)[n];
      for(int i=0; i<m.nconfigs(beta); i++){
	m.lambda[alpha][n](i) = -1 + 2*((double)rand()/(double)RAND_MAX);
      }
    }
  }

  for(int alpha=0; alpha<m.cliques.size(); alpha++){
    for(int n=0; n<m.get_kids(alpha).size(); n++){
      int beta = m.get_kids(alpha)[n];
      for(int i=0; i<m.nconfigs(beta); i++)
	assert(m.msg(alpha,beta,i) == m.msg(alpha,beta)(i));
    }
  }
}
*/
