#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>

double uni_differential (double x1, double y1, double t1, double x2, double y2, double t2,
                         double sigma_square, double scale, double nu, double vx, double vy,
                         int nonstat) {

  double l1x_new, l1y_new, l2x_new, l2y_new, xlag, ylag;
  double con, expr = 0.0, cov_val;

  l1x_new = x1 - vx * t1;
  l1y_new = y1 - vy * t1;
  l2x_new = x2 - vx * t2;
  l2y_new = y2 - vy * t2;
  
  if(nonstat == 0){
    xlag = l1x_new - l2x_new;
    ylag = l1y_new - l2y_new;
    
    expr = scale * sqrt(pow(xlag, 2) + pow(ylag, 2));
  }else if(nonstat == 1){
    double bx = 0.15, by = 0.15;
    
    double dist_from_source_l1 = sqrt(pow(l1x_new - bx, 2) + pow(l1y_new - by, 2));
    double dist_from_source_l2 = sqrt(pow(l2x_new - bx, 2) + pow(l2y_new - by, 2));
    
    l1x_new = bx + (l1x_new - bx) * (1 + 2 * exp(-0.5 * pow(dist_from_source_l1, 2)));
    l1y_new = by + (l1y_new - by) * (1 + 2 * exp(-0.5 * pow(dist_from_source_l1, 2)));
    l2x_new = bx + (l2x_new - bx) * (1 + 2 * exp(-0.5 * pow(dist_from_source_l2, 2)));
    l2y_new = by + (l2y_new - by) * (1 + 2 * exp(-0.5 * pow(dist_from_source_l2, 2)));
    
    xlag = l1x_new - l2x_new;
    ylag = l1y_new - l2y_new;
    
    expr = scale * sqrt(pow(xlag, 2) + pow(ylag, 2));
    
  }else if(nonstat == 2){
    
    double omega1 = 0 + 1 * (l1x_new - .5) + 1 * (l1y_new - .5) + 3 * pow(l1x_new - .5, 2) + -1 * pow(l1y_new - .5, 2);
    double lam1_1 = exp(- pow(l1x_new, 2)) - exp(- pow(l1y_new, 2));
    double lam1_2 = sin(l1x_new) * exp(- pow(l1y_new, 2)) - sin(l1y_new) * exp(- pow(l1x_new, 2));
    
    double omega2 = 0 + 1 * (l2x_new - .5) + 1 * (l2y_new - .5) + 3 * pow(l2x_new - .5, 2) + -1 * pow(l2y_new - .5, 2);
    double lam2_1 = exp(- pow(l2x_new, 2)) - exp(- pow(l2y_new, 2));
    double lam2_2 = sin(l2x_new) * exp(- pow(l2y_new, 2)) - sin(l2y_new) * exp(- pow(l2x_new, 2));
    
    double D1_11 = exp(lam1_1) * cos(omega1) * cos(omega1) + exp(lam1_2) * sin(omega1) * sin(omega1);
    double D1_12 = exp(lam1_1) * cos(omega1) * sin(omega1) - exp(lam1_2) * sin(omega1) * cos(omega1);
    double D1_22 = exp(lam1_1) * sin(omega1) * sin(omega1) + exp(lam1_2) * cos(omega1) * cos(omega1);
    
    double D2_11 = exp(lam2_1) * cos(omega2) * cos(omega2) + exp(lam2_2) * sin(omega2) * sin(omega2);
    double D2_12 = exp(lam2_1) * cos(omega2) * sin(omega2) - exp(lam2_2) * sin(omega2) * cos(omega2);
    double D2_22 = exp(lam2_1) * sin(omega2) * sin(omega2) + exp(lam2_2) * cos(omega2) * cos(omega2);
    
    double D1_det = D1_11 * D1_22 - D1_12 * D1_12;
    double D2_det = D2_11 * D2_22 - D2_12 * D2_12;
    
    double D_11 = 0.5 * (D1_11 + D2_11);
    double D_12 = 0.5 * (D1_12 + D2_12);
    double D_22 = 0.5 * (D1_22 + D2_22);
    double D_det = D_11 * D_22 - D_12 * D_12;
    
    double DInv_11 = D_22; 
    double DInv_22 = D_11;
    double DInv_12 = - D_12; 
    
    double sigma = sqrt(sqrt(D1_det * D2_det)/D_det);
    
    sigma_square = sigma_square * pow(sigma, 2);
    
    xlag = l1x_new - l2x_new;
    ylag = l1y_new - l2y_new;
    
    expr = scale * sqrt(xlag * xlag * DInv_11 + xlag * ylag * DInv_12 + xlag * ylag * DInv_12 + ylag * ylag * DInv_22) / D_det;
    
  }
  
  con = pow(2, nu - 1) * tgamma(nu);
  con = 1.0 / con;
  con = sigma_square * con;
  
  if(expr == 0){
    cov_val = sigma_square;
  }else{
    cov_val = con * pow(expr, nu) * gsl_sf_bessel_Knu(nu, expr);
  }

  return(cov_val);
}

void CovUniLagrangian (double * x, double * y, double * t,
                       double *sigma2, double *scale, double *nu,
                       double *vx, double *vy, 
                       double * dist, int *n, int *nonstat) {
  int nn = *n, i, j, nonstatt = *nonstat;

  double NU = *nu, SIG2 = *sigma2, SCALE = *scale, VX = *vx, VY = *vy;

  for (i = 0; i < nn; i++) {
    for (j = 0; j < nn ; j++) {
      dist[nn * i + j] = uni_differential(x[i], y[i], t[i], x[j], y[j], t[j],
                                              SIG2, SCALE, NU, VX, VY, nonstatt);
    }
  }
}



