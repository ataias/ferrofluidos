#include <FFH.hh>
#include <iostream>

void navier_stokes_step1(int n, double dt, double mu, double rho, MatrixXd *u, MatrixXd *v, MatrixXd *u_old, MatrixXd *v_old, MatrixXd *fx, MatrixXd *fy){
	//u and v will be modified!!!
	int i, j; double u_s, v_s, u_t, v_t;
	double dx = 1.0/(n-1.0); double dx2 = dx*dx;
	double f_s; double dtx = dt/(2.*dx);

	//As far as I am concerned, we have dirichlet boundary conditions on u and v
	for(i=1; i<n-1; i++){
		for(j=1; j<n-1; j++){
			//Those two lines are just summing the stencil around i,j
			u_s  = (*u_old)(i+1,j)+(*u_old)(i-1,j)+(*u_old)(i,j+1)+(*u_old)(i,j-1);
			v_s  = (*v_old)(i+1,j)+(*v_old)(i-1,j)+(*v_old)(i,j+1)+(*v_old)(i,j-1);

			//Those lines take an average of u around i,j considering the 
			// staggered grid
			u_t  = (*u_old)(i,j)+(*u_old)(i+1,j)+(*u_old)(i+1,j-1)+(*u_old)(i,j-1);
			u_t *= 0.25;
			v_t  = (*v_old)(i,j)+(*v_old)(i-1,j)+(*v_old)(i-1,j+1)+(*v_old)(i,j+1);
			v_t *= 0.25;

			//Now, let's compute u* in i,j
			//Order is important here!
			f_s = (*fx)(i,j)+(*fx)(i-1,j);
			(*u)(i,j)  =  mu*(u_s-4*(*u_old)(i,j))/dx2+f_s;
			(*u)(i,j) *= dt/rho; 
			(*u)(i,j) += -dtx*(*u_old)(i,j)*((*u_old)(i+1,j)-(*u_old)(i-1,j));
			(*u)(i,j) +=  dtx*v_t*((*u_old)(i,j+1)-(*u_old)(i,j-1));
			(*u)(i,j) +=  (*u_old)(i,j);
			
			//Now, let's compute v* in i,j
			f_s = (*fy)(i,j)+(*fy)(i-1,j);
			(*v)(i,j)  =  mu*(v_s-4*(*v_old)(i,j))/dx2+f_s;
			(*v)(i,j) *= dt/rho; 
			(*v)(i,j) += -dtx*u_t*((*v_old)(i+1,j)-(*v_old)(i-1,j));
			(*v)(i,j) +=  dtx*(*v_old)(i,j)*((*v_old)(i,j+1)-(*v_old)(i,j-1));
		}
	}
}