/*
 * Universidade de Brasília - UnB
 * Departamento de Matemática - MAT
 * Departamento de Engenharia Mecânica - ENM
 * Projeto de Iniciação Científica - 2012/2013
 * Aluno: Ataias Pereira Reis
 */

/*
 * Description of the problem in these lines
 */
#include<OwnMath.hpp>

using namespace std;
using namespace Eigen;

int main()
{
	/*This programs */
	time_t start;
	clock_t tStart = clock();

	/*------Time Operations-----*/
	time(&start);
	tStart = clock();
	/*--------------------------*/
	MatrixATd conditions;
	char file[] = "file1.dat";
	conditions = dnc_U(conditions, file);
	/*------------------------------------------------*/
	cout << "Teste" << endl;
	MatrixATd *U, *Uold, ER;
	int N_t = 1000;
	double t_max = 0.1;
	double L = 1;
	double dt = 1.0*t_max / (double) (N_t-1);
	double dx =  1.0*L / (double) (AT-1); //de 0 a L divido para para ter 10 de dimensão, a divisão se faz por AT-1 ao invés de AT
	U = (MatrixATd*) malloc (N_t*sizeof(*U));
	Uold = (MatrixATd*) malloc (N_t*sizeof(*Uold));
	U[0] = conditions;
	Uold[0] = conditions;
	double K = 1;
	double alpha = 1/(1/dt + 4*K/(dx*dx));
	bool STOP = false;
	double r = 0;
	VectorXd element(N_t);
	element(0)=U[0](5,5);
	for(unsigned k=0; k<N_t-1; k++){
		do{
			for(unsigned i = 1; i<AT-1; i++)
				for(unsigned j = 1; j<AT-1; j++)
					U[k+1](i,j) = r*Uold[k+1](i,j)  + (1-r)*alpha*( K*(Uold[k+1](i+1,j) + Uold[k+1](i-1,j) + Uold[k+1](i,j+1) + Uold[k+1](i,j-1))/(dx*dx) + U[k](i,j)/dt);
				ER = U[k+1]-Uold[k+1];
				if(ER.norm() < 0.0000001) STOP = true;
				Uold[k+1] = U[k+1];
		}while(!STOP);
		element(k+1)= U[k+1](5,5);
		STOP=false;
	}
	/*------------------*/

	/*-----File operations------*/
    fileName(file); /*Gives a string to file according to N*/
	save(file, tStart, &start, element, "u_xx+u_yy=f(x,y)"); /*Save file in .m file to plot graph*/
	/*--------------------------*/
	return(0);
}
