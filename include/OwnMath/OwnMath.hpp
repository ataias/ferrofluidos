/**
 * @file OwnMath.hpp
 * @author  Ataias Pereira Reis <ataiasreis@gmail.com>
 * @version 1.0.1
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * A classe Time representa um momento de tempo.
 */

#include<iostream>
#include<vector>
#include<Dense>
#include<Eigen>
#include<Sparse>
#include<SparseCore>
#include<SparseCholesky>
#include<Core>
#include<Eigenvalues>
#include<Geometry>

#include<fstream>
#include<string>
#include<cstdlib>
#include<ctime>

using namespace std;
using namespace Eigen;
/*Global Variables*/
    int N = 21;
#define AT 21
typedef Matrix<double, 50, 50> Matrix50d;
typedef Matrix<double, AT, AT> ArrayATd;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

MatrixXd dc_U(MatrixXd U,
		      double u1,
		      double u2,
		      double u3,
		      double u4
		      )
{
	/**
	 * The function dc_U is made to create a dirichlet boundary condition,
	 * one with each side constant.
	 *
	 * This considers a square matrix, with dx = dy.
	 *
	 * @param  U is the matrix to which the boundary will be changed
	 * @param  u1 is the value of column 0
	 * @param  u2 is the value of row 0
	 * @param  u3 is the value of row(N-1)
	 * @param  u4 is the value of column(N-1)
	 * @return The matrix U is returned, with changed boundary conditions
	 */
	U = MatrixXd::Zero(N,N);
	U.col(0) = VectorXd::Constant(N,u1);
	U.row(0) = RowVectorXd::Constant(N,u2);
	U.row(N-1) = RowVectorXd::Constant(N,u3);
	U.col(N-1) = VectorXd::Constant(N,u4);

	return(U);
}

MatrixXd dnc_U(MatrixXd U,
		       char file[]
		       )
{
	/**
	 * The function dcn_U is made to create a dirichlet boundary condition,
	 * one with each side. Actually it could be any condition, it just set
	 * the boundary of the matrix. Beside that, it takes the whole matrix,
	 * so includes any initial condition in the file to be read.
	 *
	 * This considers a square matrix, with dx = dy. The values are arbitrary
	 * and taken from a text file with the matrix. It does not require the number
	 * of columns of lines.
	 *
	 * @param  U is the matrix to which the boundary will be changed
	 * @param  file[] is the name of text file
	 * @return The matrix U is returned, with changed boundary conditions
	 */
	ifstream m_file;
	m_file.open(file);
	int rn = 0, cn=0; /*!< var rn, cn :
	 	 	 	 	 	 \brief rows and columns number 青い */

	/*Loop to discover number of rows and columns while reading file*/
	string line;
	while(getline(m_file, line)){
		rn++;
		}
	m_file.close();

	m_file.open(file);
	double num;
	while(!m_file.eof()){
		m_file >> num;
		cn++;
	}
	cn = cn/rn;
	m_file.close();
	U.resize(rn,cn);
	m_file.open(file);
	for (int i = 0; i < rn; i++) {
		for (int j = 0; j < cn; j++) {
			m_file >> U(i,j);
	    }
	  }
	m_file.close();
	return(U);
}

//g is the matrix of non-homogeneity
MatrixXd PoissonSparseSI(MatrixXd g,
						 MatrixXd U,
						 double dt,
						 double dx,
						 double K
						 )
{
	/**
	 * The function PoissonSparseSI is made to solve the Poisson Equation
	 * using an implicit method of solution
	 *
	 * Equation is: \f$K\nabla^2 U=g(x,y)\f$
	 * @param  g is the non-homogeneity
	 * @param  U is the matrix where the solution will be save
	 * @param  dt is the time step
	 * @param  dx is the spatial step. In y, the step is dy = dx
	 * @param  K is the coefficient in the formula.
	 * @return The matrix U is returned, as solution of the given equation
	 */
	std::vector<T> tripletList; //triplet, needed to fill sparse matrix
	long unsigned n = (N-2)*(N-2);
	tripletList.reserve(5*n);
	long unsigned i,j, k=1;
	MatrixXd ind = MatrixXd::Zero(n,2);
	VectorXd b(n);
	U = g; //Initialize U
	/*Important coefficients of matrix and vector b*/
	double alpha = dt*K/(dx*dx);
	double beta = 1+4*dt*K/(dx*dx);
	/*Starting to build A*/
	for(i = 0; i<n; i++)
		tripletList.push_back(T(i,i,-1*beta));
	for(i = 0; i<= n-2; i++)
		if( (i+1) % (N-2) != 0) tripletList.push_back(T(i,i+1,alpha));
	for(i = 0; i<= n-2; i++)
	    if( (i+1) % (N-2) != 0) tripletList.push_back(T(i+1,i,alpha));
	for(i = 0; i < (n-(N-2)); i++)
		tripletList.push_back(T(i,i+N-2,alpha));
	for(i = 0; i<(n-(N-2)); i++)
		tripletList.push_back(T(i+N-2,i,alpha));
	SpMat A(n,n);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	/*Ending A*/
	for(i = 1; i<= (N-2); i++)
	    for(j=1; j<=(N-2); j++){
	    	/*ind has the indexes of U in a matrix with two columns*/
	        ind(j+(i-1)*(N-2)-1,0)= i+1;
	        ind(j+(i-1)*(N-2)-1,1)= j+1;
	    }


	for(i = 0; i<n; i++){
	    double vb = 0; /*Value of i_th element in b*/
	    if(ind(i,0)-1 == 0) vb=vb+U(ind(i,1)-1,0);
	    if(ind(i,0)+1 == N-1) vb=vb+U(ind(i,1)-1,N-1);
	    if(ind(i,1)-1 == 0) vb=vb+U(N-1,ind(i,0)-1);
	    if(ind(i,1)+1 == N-1) vb=vb+U(0,ind(i,0)-1);
	    b(i) = -1.0*g(ind(i,0)-1,ind(i,1)-1)-1*alpha*vb;
	    if(vb !=0.0) cout << vb << "heeeeeeeeeeeeey, vb..." << endl;
	}

	cout << U(10,10) << "\n";

	Eigen::SimplicialCholesky<SpMat> chol(A);
	VectorXd u = chol.solve(b);
	for(i=1; i<=(N-2); i++)
	    for(j=1; j<=(N-2); j++){
	    	U(i,j)=u(3*i+j-4);
	    }
	return(U);
}

typedef Matrix<double, AT, AT> MatrixATd;

MatrixATd* SolveHeatEquationE(
					MatrixATd conditions,
					double alpha, // u_t = alpha*laplacian(u)
					int N_t, //Number of points in time variable
					double L, //length of interval in x and y;
					double t_max //max time
				  )
{
	double dt = t_max / (N_t);
	double dx = (1.0*L)/(AT-1);
	double r = 1.0*alpha*dt/(dx*dx);
	double h = 1.0*1-4.0*r;

	MatrixATd *U;
	U = (MatrixATd*) malloc (N_t*sizeof(*U));
	U[0] = conditions;

	for(int k=1; k<N_t; k++){
		U[k] = U[k-1];
    for(int i=1; i<AT-1; i++){
      for(int j=1; j<AT-1; j++)
        U[k](i,j) = r*U[k-1](i-1,j) + h*U[k-1](i,j) + r*U[k-1](i+1,j)+r*U[k-1](i,j-1)+r*U[k-1](i,j+1);
    }

	}
	return(U);
}

MatrixATd* SolveHeatEquationI(
					MatrixATd conditions,
					double K, // u_t = alpha*laplacian(u)
					int N_t, //Number of points in time variable
					double L, //length of interval in x and y;
					double t_max //max time
				  ){

	double dt = 1.0*t_max / (double) (N_t-1);
	double dx =  1.0*L / (double) (AT-1); //de 0 a L divido para para ter 10 de dimensão, a divisão se faz por AT-1 ao invés de AT

	MatrixATd *U;
	U = (MatrixATd*) malloc (N_t*sizeof(*U));
	U[0] = conditions;

	for(int k=1; k<N_t; k++){
        U[k] = PoissonSparseSI(U[k-1], U[k], dt, dx, K);
    }

	return(U);
}


using namespace Eigen;
using namespace std;

MatrixXd Poisson(MatrixXd g){ //seria bom checar se g é NxN

	long unsigned n = (N-2)*(N-2);
	long unsigned i,j, k=1;
	double dx = 1/N; // malha é de 0 a 1 em x e y
//	SpMat A(n,n);
	MatrixXd A = MatrixXd::Zero(n,n);
	MatrixXd U = MatrixXd::Zero(N,N);
	MatrixXd aux = MatrixXd::Zero(n,2);
//	SpMat b(n,1);
	MatrixXd b = MatrixXd::Zero(n,1);
	MatrixXd u = MatrixXd::Zero(n,1);

	for(i = 0; i<n; i++)
	    A(i,i) = -4;

	for(i = 0; i<= n-2; i++){
	    A(i,i+1) = 1;
	    if( (i+1) % (N-2) == 0) A(i,i+1) = 0;
	}
	for(i = 0; i<= n-2; i++){
	    A(i+1,i) = 1;
	    if((i+1) % (N-2) == 0) A(i+1,i) = 0;
	}

	for(i = 0; i < (n-(N-2)); i++)
	    A(i,i+N-2) = 1;

	for(i = 0; i<(n-(N-2)); i++)
	    A(i+N-2,i) = 1;

	//Condições de Contorno
	for(i = 0; i<N; i++){
	    U(i,0) = 0;
	    U(0,i) = 1;
	    U(N-1,i) = 0;
	    U(i,N-1) = 0;
	}
	// Esta matriz tem as coordenadas dos U interiores ordenadas
	// em um vetor coluna
	for(i = 1; i<= (N-2); i++)
	    for(j=1; j<=(N-2); j++){
	        aux(j+(i-1)*(N-2)-1,0)= j+1;
	        aux(j+(i-1)*(N-2)-1,1)= i+1;
	    }

	for(i = 0; i<n; i++){
	    double valor = 0;
	    double x = aux(i,0);
	    double y = aux(i,1);
	    if(x-1 == 1)
	        valor=valor+U(y-1,0);

	    if(x+1 == N)
	        valor=valor+U(y-1,N-1);
	    if(y-1 == 1)
	        valor=valor+U(N-1,x-1);
	    if(y+1 == N)
	        valor=valor+U(0,x-1);

	    b(i,0) = dx*dx*g(x,y)-1*valor;
	}

	u = A.ldlt().solve(b); //Solve system
	for(i=2; i<=(N-1); i++)
	    for(j=2; j<=(N-1); j++)
	        U(N+1-i-1,j-1) = u((i-2)*(N-2)+j-1-1);

	return(U);
}

MatrixXd PoissonSparse(MatrixXd g, MatrixXd U){

	std::vector<T> tripletList; //triplet, needed to fill sparse matrix
	long unsigned n = (N-2)*(N-2);
	tripletList.reserve(5*n);
	long unsigned i,j, k=1;
	double dx = 1.0/N; /*Mesh has domain 0 < x < 1 and 0 < y < 1*/
	MatrixXd ind = MatrixXd::Zero(n,2);
	VectorXd b(n);


	/*Starting to build A*/
	for(i = 0; i<n; i++)
		tripletList.push_back(T(i,i,-4));
	for(i = 0; i<= n-2; i++)
		if( (i+1) % (N-2) != 0) tripletList.push_back(T(i,i+1,1));
	for(i = 0; i<= n-2; i++)
	    if( (i+1) % (N-2) != 0) tripletList.push_back(T(i+1,i,1));
	for(i = 0; i < (n-(N-2)); i++)
		tripletList.push_back(T(i,i+N-2,1));
	for(i = 0; i<(n-(N-2)); i++)
		tripletList.push_back(T(i+N-2,i,1));
	SpMat A(n,n);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	/*Ending A*/

	for(i = 1; i<= (N-2); i++)
	    for(j=1; j<=(N-2); j++){
	    	/*ind has the indexes of U in a matrix with two columns*/
	        ind(j+(i-1)*(N-2)-1,0)= j+1;
	        ind(j+(i-1)*(N-2)-1,1)= i+1;
	    }

	for(i = 0; i<n; i++){
	    double vb = 0; /*Value of i_th element in b*/
	    if(ind(i,0)-1 == 1) vb=vb+U(ind(i,1)-1,0);
	    if(ind(i,0)+1 == N) vb=vb+U(ind(i,1)-1,N-1);
	    if(ind(i,1)-1 == 1) vb=vb+U(N-1,ind(i,0)-1);
	    if(ind(i,1)+1 == N) vb=vb+U(0,ind(i,0)-1);
	    b(i) = dx*dx*g(ind(i,0),ind(i,1))-1*vb;
	}

	Eigen::SimplicialCholesky<SpMat> chol(A);
	VectorXd u = chol.solve(b);
	for(i=2; i<=(N-1); i++)
	    for(j=2; j<=(N-1); j++)
	        U(N+1-i-1,j-1) = u((i-2)*(N-2)+j-1-1);
	return(U);
}
//---------------------------------

/*
 *
 *
 *
 */
void fileName(char* file)
{
	/**
	 * Method for creating a file name according to matrix size
	 *
	 * Matrix size is not a parameter of function, it is set as a global variable
	 * @param file is a char pointer where the string will be saved
	 * @return No return value.
	 */
    file[0] = 'm';
    file[1] = 'e';
    file[2] = 's';
    file[3] = 'h';
    int m = N/100;
    file[4] = '0' + m;
    m = N - 100*m;
    file[5] = '0' + m/10;
    m = m%10;
    file[6] = '0' + m;
    file[7] = '.';
    file[8] = 'm';
    file[9] = '\0';
}
void save(char file[],
	      clock_t tStart,
	      time_t *inicio,
	      MatrixXd U,
	      string equation
	      )
{
	/**
	 * Method for saving a matrix in a file .m
	 * After save, it is possible to open easily with octave
	 * and print the graph for file
	 *
	 * @param  file[] name of file where matrix will be save
	 * @param  tStart has the value of cpu time in the start, in this method
	 * the "tFinal" is calculated and subtracted from tStart
	 * @param  *inicio is a pointer that indicates the time of start, and date
	 * @param U is the matrix that will be saved
	 * @param equation is a string that receives the equation being solved,
	 * just to show in the .m file. Ex: Laplace Equation - u_xx + u_yy = f(x,y)
	 * @return No return value.
	 */
  ofstream m_file;
  m_file.open(file);
  m_file.unsetf(ios::floatfield);
  m_file.precision(15);
  m_file << "%_ Start time: " << ctime(inicio) << endl;
  m_file << "%_ Solving equation " << equation << endl;
  m_file << "%_ Square Mesh" << endl;
  m_file << "%_ Domain in x and y: 0<x<1 e 0<y<1 " << endl;
  m_file << "%_ N = " << N << ". Spent time: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s." << endl;
  m_file << "N = " << N << ";" << endl;
  m_file << endl;
  m_file << "U = [" << endl;
  m_file << U << endl;
  m_file << "];" << endl;
  m_file << "X = 0:1/(N-1):1;" << endl;
  m_file << "Y = 0:1/(N-1):1;" << endl;
  m_file << "figure;" << endl;
  m_file << "mesh(X,Y,U);" << endl;
  m_file << "hold on" << endl;
  m_file << "title('Solution of L(U) = f(x,y) with ";
  m_file << N << "x" << N << " points')" << endl;
  m_file << "xlabel('x')" << endl;
  m_file << "ylabel('y')" << endl;
  m_file << "zlabel('u(x,y)');" << endl;
  m_file << "mkdir png;" << endl;
  m_file << "hold off;" << endl;
  int a = N/100;
  int b = N - 100*a;
  int c = b%10;
  m_file << "print -dpng 'png/mesh" << a << b/10 << c << ".png';" << endl;
  time_t fim = time(NULL);
  m_file << "%_ End of operation: " << ctime(&fim) << endl;
  m_file.close();
  cout << "Mesh "<< N << "x" << N <<". Succeed to write in file. Operation time: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s." << endl;
}


void saveTmEqn(char file[],
			   clock_t tStart,
			   time_t *inicio,
			   string equation,
			   MatrixXd *U,
			   int N
			   )
{
	/**
	 * Method for saving various matrices in a file .m
	 * After save, it is possible to open easily with octave
	 * and make a video of the matrices
	 *
	 * @param  file[] name of file where matrix will be save
	 * @param  tStart has the value of cpu time in the start, in this method
	 * the "tFinal" is calculated and subtracted from tStart
	 * @param  *inicio is a pointer that indicates the time of start, and date
	 * @param equation is a string that receives the equation being solved,
	 * just to show in the .m file. Ex: Laplace Equation - u_xx + u_yy = f(x,y)
	 * @param *U is a pointer to the 3d matrix U, that will be saved
	 * @param N indicates how many two-dimensional matrices there are in matrix U
	 * @return No return value.
	 */
	/*
  ofstream m_file;
  m_file.open(file);
  m_file.unsetf(ios::floatfield);
  m_file.precision(15);
  m_file << "%_ Start time: " << ctime(inicio) << endl;
  m_file << "%_ Solving equation u_xx + u_yy = u_t" << endl;
  m_file << "%_ Square Mesh" << endl;
  m_file << "%_ Domain in x and y: 0<x<1 e 0<y<1 " << endl;
  m_file << "%_ N = " << N << ". Spent time: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s." << endl;
  m_file << "N = " << N << ";" << endl;
  m_file << endl;

  int frames = time / dt;
  for(int i=1; i<frames; i++){
	  m_file << "U(:,:," << i << ") = [" << endl;
	  m_file << U.block<> << endl;
  }
  m_file << "];" << endl;
  m_file << endl;

  m_file << "X = 0:(1/N):(1-(1/N)/8);" << endl;
  m_file << "Y = 0:(1/N):(1-(1/N)/8);" << endl;
  m_file << "figure;" << endl;
  m_file << "mesh(X,Y,U);" << endl;
  m_file << "hold on" << endl;
  m_file << "title('Solution of L(U) = f(x,y) with ";
  m_file << N << "x" << N << " points')" << endl;
  m_file << "xlabel('x')" << endl;
  m_file << "ylabel('y')" << endl;
  m_file << "zlabel('u(x,y)');" << endl;
  m_file << "mkdir png;" << endl;
  int a = N/100;
  int b = N - 100*a;
  int c = b%10;
  m_file << "print -dpng 'png/mesh" << a << b/10 << c << ".png';" << endl;
  time_t fim = time(NULL);
  m_file << "%_ End of operation: " << ctime(&fim) << endl;
  m_file.close();
  cout << "Succeed to write in file. Operation time: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s." << endl;
  */
}

