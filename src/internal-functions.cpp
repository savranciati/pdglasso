
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

//' @title Half-vectorization
//' @description Internal function. Half-vectorization operator, without the diagonal.
//' @param M symmetric matrix.
//' @return a vector.
//' @noRd
// [[Rcpp::export]]
NumericVector half_vec(const NumericMatrix& M) {
  int p = M.nrow();
  int n = p * (p - 1) / 2;
  NumericVector out(n);
  
  int k = 0;
  for (int j = 1; j < p; j++) {
    for (int i = 0; i < j; i++) {
      out[k++] = M(i, j);
    }
  }
  
  return out;
}




//' @title From symmetric matrix to vector
//' @description Internal function. Transforms a symmetric matrix into a vector.
//' @param M symmetric matrix.
//' @return a vector.
//' @noRd
// [[Rcpp::export]]
NumericVector mat2vec(const NumericMatrix& M) {
  int p = M.nrow();
  int q = p / 2;
  int dim_hs = q * (q - 1) / 2;
  int out_size = 2 * q + 4 * dim_hs + q;
  NumericVector out(out_size);
  
  int k = 0; 
  
  // diag of top-left block
  for (int i = 0; i < q; i++) out[k++] = M(i, i);
  
  // diag of bottom-right block
  for (int i = q; i < p; i++) out[k++] = M(i, i);
  
  // upper-triangular of top-left block
  for (int j = 1; j < q; j++)
    for (int i = 0; i < j; i++)
      out[k++] = M(i, j);
  
  // upper-triangular of bottom-right block
  for (int j = q + 1; j < p; j++)
    for (int i = q; i < j; i++)
      out[k++] = M(i, j);
  
  // top-right block (rows 0..q-1, cols q..p-1)
  for (int j = q; j < p; j++)
    for (int i = 0; i < q && i < j - q; i++)
      out[k++] = M(i, j);
  
  // bottom-left block 
  for (int j = 0; j < q; j++)
    for (int i = q; i < p && i - q < j; i++)
      out[k++] = M(i, j);
  
  // 7. diag of top-right block
  for (int i = 0; i < q; i++)
    out[k++] = M(i, i + q);
  
  return out;
}




// Inverse operation with respect to mat2vec()
NumericMatrix vec2mat(const NumericVector& m) {
  int d = m.size();
  int p = (-1 + std::sqrt(1 + 8.0 * d)) / 2.0;
  int q = p / 2;
  int dim_hs = q * (q - 1) / 2;
  
  int i1 = 0;
  int i2 = i1 + q;
  int i3 = i2 + q;
  int i4 = i3 + dim_hs;
  int i5 = i4 + dim_hs;
  int i6 = i5 + dim_hs;
  int i7 = i6 + dim_hs;
  
  NumericMatrix M(p, p);
  
  // top-left block 
  for (int j = 1; j < q; j++) {
    int offset = j*(j-1)/2;
    for (int i = 0; i < j; i++) {
      M(i, j) = m[i3 + offset + i];
      M(j, i) = M(i, j);  
    }
  }
  for (int i = 0; i < q; i++) M(i, i) = m[i1 + i];
  
  // bottom-right block 
  for (int j = q+1; j < p; j++) {
    int offset = (j-q)*(j-q-1)/2;
    for (int i = q; i < j; i++) {
      M(i, j) = m[i4 + offset + (i-q)];
      M(j, i) = M(i, j);  // symmetric
    }
  }
  for (int i = q; i < p; i++) M(i, i) = m[i2 + i - q];
  
  // top-right & bottom-left block 
  
  for (int j = 1; j < q; j++) {
    int offset = j*(j-1)/2;
    for (int i = 0; i < j; i++) {
      M(i, j+q) = m[i5 + offset + i];
      M(j+q, i) = M(i, j+q);  // bottom-left mirror
    }
  }
  
  for (int i = 1; i < q; i++) {
    int offset = i*(i-1)/2;
    for (int j = 0; j < i; j++) {
      M(i, j+q) = m[i6 + offset + j];
      M(j+q, i) = M(i, j+q);  // top-right mirror
    }
  }
  
  // top-right block
  for (int i = 0; i < q; i++) {
    M(i, i+q) = m[i7 + i];
    M(i+q, i) = M(i, i+q);  
  }
  
  return M;
}




// Computes F%*%v
NumericVector F_by_vec(const NumericVector& v, int p, const std::string& acr_type) {
  int q = p / 2;
  int dim_hs = q * (q - 1) / 2;
  
  int i1 = 0, i2 = i1 + q;
  int i3 = i2 + q, i4 = i3 + dim_hs;
  int i5 = i4 + dim_hs, i6 = i5 + dim_hs;
  int i7 = i6 + dim_hs;
  
  int size_v1 = i2 - i1;
  int size_v2 = i4 - i3;
  int size_v3 = i6 - i5;
  
  // Lambda to fill blocks
  auto diff_block = [&](int start1, int start2, int len, NumericVector& out, int& k) {
    for (int i = 0; i < len; i++) out[k++] = v[start1 + i] - v[start2 + i];
  };
  
  // output sizes
  int out_size = 0;
  if (acr_type.find("V") != std::string::npos) out_size += size_v1;
  if (acr_type.find("I") != std::string::npos) out_size += size_v2;
  if (acr_type.find("A") != std::string::npos) out_size += size_v3;
  
  NumericVector out(out_size);
  int k = 0;
  
  if (acr_type.find("V") != std::string::npos) diff_block(i1, i2, size_v1, out, k);
  if (acr_type.find("I") != std::string::npos) diff_block(i3, i4, size_v2, out, k);
  if (acr_type.find("A") != std::string::npos) diff_block(i5, i6, size_v3, out, k);
  
  return out;
}







// Computes t(F)%*%v
NumericVector tF_by_vec(const NumericVector& v, int p, const std::string& acr_type) {
  int q = p / 2;
  int dim_hs = q * (q - 1) / 2;
  int dim_h  = p * (p + 1) / 2;
  
  NumericVector out(dim_h, 0.0);
  int k = 0;  
  
  // Lambda to handle a block
  auto assign_block = [&](int start, int block_size) {
    for (int i = 0; i < block_size; i++, k++) {
      out[start + i] = v[k];
      out[start + block_size + i] = -v[k];
    }
  };
  
  if (acr_type.find("V") != std::string::npos) assign_block(0, q);
  if (acr_type.find("I") != std::string::npos) assign_block(2 * q, dim_hs);
  if (acr_type.find("A") != std::string::npos) assign_block(2 * q + 2 * dim_hs, dim_hs);
  
  return out;
}




// Inner ADMM loop
NumericMatrix admm_inner(const NumericMatrix& X,
                         const NumericMatrix& U,
                         double rho1,
                         NumericVector lambda1,
                         const NumericVector& lambda2,
                         double rho2_init,
                         bool varying_rho2,
                         int max_iter_int,
                         double eps_abs,
                         double eps_rel,
                         int n_row_F,
                         const std::string& acr_type) {
  
  int p = X.nrow();
  
  // precompute X + U 
  NumericMatrix XU(p,p);
  for(int i=0;i<p;i++)
    for(int j=0;j<p;j++)
      XU(i,j) = X(i,j) + U(i,j);
  
  NumericVector b = mat2vec(XU);
  int d = b.size();
  
  NumericVector x(d, 1.0);
  NumericVector v(n_row_F, 0.0), t(n_row_F, 0.0);
  NumericVector last_v(n_row_F);
  
  NumericVector lambda(n_row_F);
  if(lambda2.size() == 1) std::fill(lambda.begin(), lambda.end(), lambda2[0]/rho1);
  else for(int i=0;i<n_row_F;i++) lambda[i] = lambda2[i]/rho1;
  
  // ADMM parameters
  double mu = 10.0, tau_inc = 2.0, tau_dec = 2.0;
  double rho2 = rho2_init;
  double alpha = rho2 / (1.0 + 2.0*rho2);
  
  // preallocate temp vectors 
  NumericVector Fb(d), rhs(n_row_F), tFrhs(d), Fx(n_row_F), Fx_tmp(n_row_F), tv(d), tF_rho2t(d);
  
  for(int kk=0; kk<max_iter_int; kk++){
    
    // x update
    Fb = F_by_vec(b, p, acr_type);
    for(int i=0;i<n_row_F;i++) rhs[i] = v[i] - t[i] - Fb[i];
    tFrhs = tF_by_vec(rhs, p, acr_type);
    for(int i=0;i<d;i++) x[i] = alpha * tFrhs[i] + b[i];
    
    // Fx update
    Fx_tmp = F_by_vec(x, p, acr_type);
    for(int i=0;i<n_row_F;i++) Fx[i] = Fx_tmp[i] + t[i];
    
    // v update (soft-threshold)
    std::copy(v.begin(), v.end(), last_v.begin());
    for(int i=0;i<n_row_F;i++){
      double z = Fx[i], th = lambda[i]/rho2;
      v[i] = (z>th)? z-th : (z<-th)? z+th : 0.0;
    }
    
    for(int i=0;i<n_row_F;i++) t[i] = Fx[i] - v[i];
    
    double r_norm = 0.0, s_norm = 0.0;
    for(int i=0;i<n_row_F;i++){
      double r = Fx_tmp[i] - v[i];
      r_norm += r*r;
      rhs[i] = v[i] - last_v[i];  // reuse rhs
    }
    r_norm = sqrt(r_norm);
    
    tv = tF_by_vec(rhs, p, acr_type);
    for(int i=0;i<d;i++) s_norm += std::pow(rho2*tv[i], 2.0);
    s_norm = sqrt(s_norm);
    
    // tolerances 
    double Fx_norm = 0.0, v_norm = 0.0;
    for(int i=0;i<n_row_F;i++){
      Fx_norm += Fx_tmp[i]*Fx_tmp[i];
      v_norm += v[i]*v[i];
    }
    double eps_pri = sqrt((double)n_row_F)*eps_abs + eps_rel * std::max(sqrt(Fx_norm), sqrt(v_norm));
    
    for(int i=0;i<n_row_F;i++) tF_rho2t[i] = t[i]*rho2;
    tF_rho2t = tF_by_vec(tF_rho2t, p, acr_type);
    double tF_norm = 0.0;
    for(int i=0;i<d;i++) tF_norm += tF_rho2t[i]*tF_rho2t[i];
    double eps_dual = sqrt((double)p*(p+1)/2.0)*eps_abs + eps_rel*sqrt(tF_norm);
    
    if(r_norm < eps_pri && s_norm < eps_dual) break;
    
    // adaptive rho2
    if(varying_rho2){
      if(r_norm > mu*s_norm){
        rho2 *= tau_inc;
        alpha = rho2 / (1.0 + 2.0*rho2);
        for(int i=0;i<n_row_F;i++) t[i] /= tau_inc;
      }
      if(s_norm > mu*r_norm){
        rho2 /= tau_dec;
        alpha = rho2 / (1.0 + 2.0*rho2);
        for(int i=0;i<n_row_F;i++) t[i] *= tau_dec;
      }
    }
  }
  
  // final threshold
  NumericVector z(d);
  
  NumericVector thr(d);
  if(lambda1.size() == 1) std::fill(thr.begin(), thr.end(), lambda1[0]/rho1);
  else for(int i=0;i<d;i++) thr[i] = lambda1[i]/rho1;
  
  for(int i=0;i<d;i++){
    double s = std::fabs(x[i]) - thr[i];
    z[i] = (s>0)? (x[i]>0 ? s : -s) : 0.0;
  }
  
  return vec2mat(z);
}




//' @title ADMM graphical lasso algorithm for coloured GGMs for paired data
//' @description Internal function. Implements a double ADMM loop.
//' @param S numeric matrix
//' @param X numeric matrix
//' @param lambda1 numeric vector
//' @param lambda2 numeric vector
//' @param rho1 scalar
//' @param rho2 scalar
//' @param varying_rho1 logical
//' @param varying_rho2 logical
//' @param max_iter integer
//' @param eps_abs scalar
//' @param eps_rel scalar
//' @param acr_type character
//' @param n_row_F integer
//' @return a list
//' @noRd
// [[Rcpp::export]]
List admm_pdglasso_internal(const NumericMatrix& S,
                            NumericMatrix X,
                            NumericVector lambda1,
                            const NumericVector& lambda2,
                            double rho1,
                            double rho2,
                            bool varying_rho1,
                            bool varying_rho2,
                            int max_iter,
                            double eps_abs,
                            double eps_rel,
                            const std::string& acr_type,
                            int n_row_F) {
  
  int p = S.nrow();
  double mu = 10.0, tau_inc = 2.0, tau_dec = 2.0;
  int n_iter = 0;
  int n_iter_rho1_update_last = 0;
  
  // preallocate 
  NumericMatrix U(p,p);
  NumericMatrix Z(p,p);
  NumericMatrix last_Z(p,p);
  NumericMatrix A(p,p);
  
  double r_norm  = 0.0;
  double s_norm  = 0.0;
  double eps_pri = 0.0;
  double eps_dual= 0.0;
  
  for(int k=0; k<max_iter; k++){
    for(int i=0;i<p;i++)
      for(int j=0;j<p;j++)
        A(i,j) = rho1*(Z(i,j)-U(i,j)) - S(i,j);
    
    // Eigen decomposition
    Eigen::Map<Eigen::MatrixXd> A_eig(as<Map<MatrixXd>>(A));
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(A_eig);
    Eigen::VectorXd dvals = solver.eigenvalues();
    Eigen::MatrixXd Q = solver.eigenvectors();
    
    // Transform eigenvalues
    dvals = (dvals.array() + (dvals.array().square() + 4.0*rho1).sqrt()) / (2.0*rho1);
    
    // Reconstruct X 
    Eigen::MatrixXd X_eig = Q * dvals.asDiagonal() * Q.transpose();
    for(int i=0;i<p;i++)
      for(int j=0;j<p;j++)
        X(i,j) = X_eig(i,j);
    
    // store Z 
    std::copy(Z.begin(), Z.end(), last_Z.begin());
    
    // inner ADMM 
    Z = admm_inner(X, U, rho1, lambda1, lambda2, rho2,
                   varying_rho2, max_iter, eps_abs, eps_rel,
                   n_row_F, acr_type);
    
    // update U
    for(int i=0;i<p;i++)
      for(int j=0;j<p;j++)
        U(i,j) += X(i,j) - Z(i,j);
    
    // residuals
    r_norm = 0.0, s_norm = 0.0;
    for(int i=0;i<p;i++)
      for(int j=i;j<p;j++){
        double dx = X(i,j)-Z(i,j);
        double dz = Z(i,j)-last_Z(i,j);
        r_norm += dx*dx;
        s_norm += dz*dz;
      }
      r_norm = sqrt(r_norm);
    s_norm = rho1 * sqrt(s_norm);
    
    // norms for tolerances 
    double X_norm=0.0,Z_norm=0.0,U_norm=0.0;
    for(int i=0;i<p;i++)
      for(int j=0;j<p;j++){
        X_norm += X(i,j)*X(i,j);
        Z_norm += Z(i,j)*Z(i,j);
        U_norm += U(i,j)*U(i,j);
      }
      X_norm = sqrt(X_norm);
    Z_norm = sqrt(Z_norm);
    U_norm = sqrt(U_norm);
    
    eps_pri  = sqrt((double)p*(p+1)/2.0)*eps_abs + eps_rel*std::max(X_norm,Z_norm);
    eps_dual = sqrt((double)p*(p+1)/2.0)*eps_abs + eps_rel*rho1*U_norm;
    
    n_iter = k+1;
    if(r_norm<eps_pri && s_norm<eps_dual) break;
    
    // --- varying rho1 ---
    if(varying_rho1){
      double scale = rho1;
      if(r_norm > mu*s_norm/scale){
        n_iter_rho1_update_last = n_iter;
        rho1 *= tau_inc;
        for(int i=0;i<p;i++)
          for(int j=0;j<p;j++)
            U(i,j) /= tau_inc;
      }
      if(s_norm/scale > mu*r_norm){
        rho1 /= tau_dec;
        for(int i=0;i<p;i++)
          for(int j=0;j<p;j++)
            U(i,j) *= tau_dec;
      }
    }
  }
  
  return List::create(
    Named("X") = X,
    Named("res_primal") = r_norm,
    Named("res_dual") = s_norm,
    Named("n_iter") = n_iter,
    Named("n_iter_rho1_update_last") = n_iter_rho1_update_last,
    Named("last_rho1") = rho1,
    Named("eps_primal") = eps_pri,
    Named("eps_dual") = eps_dual,
    Named("converged") = (n_iter<max_iter)
  );
}

