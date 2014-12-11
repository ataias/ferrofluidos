#include <FFH.hh>
#include <iostream>


void navier_stokes_step1(int n, double dt, double mu, double rho, MatrixXd *u, MatrixXd *v, MatrixXd *u_old, MatrixXd *v_old, MatrixXd *fx, MatrixXd *fy, VectorXd uB){
	//u and v will be modified!!!
	int i, j; double u_s, v_s, u_t, v_t;
	double dx = 1.0/(n-1.0); double dx2 = dx*dx;
	double dtx = dt/(2.*dx);
	MatrixXd *aux;

	//Boundary conditions, velocidades normais na malha escalonada
	for(j=1; j<n-1; j++) (*u)(1,j)=0; //esquerda
	for(j=1; j<n-1; j++) (*u)(n-1,j)=0; //direita
	for(i=1; i<n-1; i++) (*v)(i,1)=0; //embaixo
	for(i=1; i<n-1; i++) (*v)(i,n-1)=0; //em cima

		//Process internal points
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
				//value is fixed in the boundaries
				if(i!=1 && i!=(n-1)) { //in the case i==1 or n-1, it is a boundary point
					//Order is important here!
				(*u)(i,j)  =  mu*(u_s-4*(*u_old)(i,j))/dx2+(*fx)(i,j);
				(*u)(i,j) *= dt/rho; 
				(*u)(i,j) += -dtx*(*u_old)(i,j)*((*u_old)(i+1,j)-(*u_old)(i-1,j));
				(*u)(i,j) += -dtx*v_t*((*u_old)(i,j+1)-(*u_old)(i,j-1));
				(*u)(i,j) +=  (*u_old)(i,j);
				} 
				//Now, let's compute v* in i,j
				if(j!=1 && j!=(n-1)) { //in the case j==1 or n-1, it is a boundary point
				(*v)(i,j)  =  mu*(v_s-4*(*v_old)(i,j))/dx2+(*fy)(i,j);
				(*v)(i,j) *= dt/rho; 
				(*v)(i,j) += -dtx*u_t*((*v_old)(i+1,j)-(*v_old)(i-1,j));
				(*v)(i,j) += -dtx*(*v_old)(i,j)*((*v_old)(i,j+1)-(*v_old)(i,j-1));
				(*v)(i,j) += (*u_old)(i,j);
				} 
			}
		}

		//Process boundary points
		//Deveria só manter o valor anterior, (*v)(0,j)=(*v_old)(0,j) também para os outros
		//só mudará a condição de fronteira para o próxima passo de tempo, não v*!
		// for(j=1; j<n-1; j++) (*v)(0,j) = -(*v_old)(0,j); //esquerda
		// for(i=1; i<n-1; i++) (*u)(i,0) = -(*u)(i,1); //embaixo
		// for(j=1; j<n-1; j++) (*v)(n-1,j) = -(*v)(n-2,j); //direita
		// for(i=1; i<n-1; i++) (*u)(i,n-1) = 2*uB(i)-(*u)(i,n-2); // em cima

}
	
void navier_stokes_step2(int n, double dt, double mu, double rho, 
						 MatrixXd *p, //pressão, assumo que já seja alocada
						 MatrixXd *u_start, MatrixXd *v_start, //v*, obtido ignorando-se o termo envolvendo p
						 MatrixXd *fx, MatrixXd *fy
						 ){
	int i, j; double DIVij;
	double dx = 1.0/(n-1.0); double dx2 = dx*dx;
	double dtx = dt/(2.*dx);
	double sum_aux; double f_aux;

	//Processar pontos internos
	for(i=1; i<n-1; i++){
		for(j=1; j<n-1; j++){
			DIVij = (rho/dt)*((*u)(i+1,j)-(*u)(i,j)+(*v)(i,j+1)-(*v)(i,j))/dx;
			sum_aux = (*p)(i+1,j)+(*p)(i-1,j)+(*p)(i,j+1)+(*p)(i,j-1);
			(*p)(i,j) = 0.25*(sum_aux)-0.25*dx2*DIVij;
		}
	}

	//Processar fronteira i=1
	i = 1;
	for(j=1; j<n-1; j++){
		sum_aux = 2*(*u)(i,j)-5*(*u)(i+1,j)+4*(*u)(i+2,j)-(*u)(i+3,j);
		f_aux = (*fx)(i,j)+(*fx)(i+1,j);
		(*p)(i-1,j) = (*p)(i,j)-(mu/dx)*sum_aux-(rho*dx/2)*f_aux;
	}

	//Processar fronteira i=n
	i = n-2; //n é o tamanho da malha escalonada, o tamanho da malha de fato é n-1
	for(j=1; j<n-1; j++){
		sum_aux = 2*(*u)(i,j)-5*(*u)(i-1,j)+4*(*u)(i-2,j)-(*u)(i-3,j);
		f_aux = (*fx)(i,j)+(*fx)(i-1,j);
		(*p)(i+1,j) = (*p)(i,j)+(mu/dx)*sum_aux+(rho*dx/2)*f_aux;
	}


	//Processar fronteira j=1
	j = 1;
	for(i=1; i<n-1; i++){
		sum_aux = 2*(*v)(i,j)-5*(*v)(i,j+1)+4*(*v)(i,j+2)-(*v)(i,j+3);
		f_aux = (*fy)(i,j)+(*fy)(i,j+1);
		(*p)(i,j-1) = (*p)(i,j)-(mu/dx)*sum_aux-(rho*dx/2)*f_aux;
	}

	//Processar fronteira j=n
	j = n-2; //n é o tamanho da malha escalonada, o tamanho da malha de fato é n-1
	for(i=1; i<n-1; i++){
		sum_aux = 2*(*v)(i,j)-5*(*v)(i,j-1)+4*(*v)(i,j-2)-(*v)(i,j-3);
		f_aux = (*fy)(i,j)+(*fy)(i,j-1);
		(*p)(i,j+1) = (*p)(i,j)+(mu/dx)*sum_aux+(rho*dx/2)*f_aux;
	}

}