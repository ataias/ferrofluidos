#ifndef NAVIERSTOKES_HH_
#define NAVIERSTOKES_HH_

void navier_stokes_step1(int n, double dt, double mu, double rho, 
						 MatrixXd *u,     MatrixXd *v, 
						 MatrixXd *u_old, MatrixXd *v_old, 
						 MatrixXd *fx,    MatrixXd *fy
						 );

void navier_stokes_step2(int n, double dt, double mu, double rho, 
						 MatrixXd *p, //pressão 
						 MatrixXd *u_star, MatrixXd *v_start //v*, obtido ignorando-se o termo envolvendo p
						 );

#endif