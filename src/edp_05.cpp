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
#include<stdheader.hpp>

int main()
{
	/*This programs */
	time_t start;
	clock_t tStart = clock();

	/*------Time Operations-----*/
	time(&start);
	tStart = clock();
	/*--------------------------*/
	dMatrix conditions;
	char file[] = "file1.dat";
	conditions = dReadConditions(file);
	/*------------------------------------------------*/
	cout << "Teste" << endl;
	dMatrix *U, *Uold, ER;
	int N_t = 1000;
	double t_max = 0.1;
	double L = 1;
	double dt = 1.0*t_max / (double) (N_t-1);
	double dx =  1.0*L / (double) (MATRIX_ORDER-1); //de 0 a L divido para para ter 10 de dimensão, a divisão se faz por AT-1 ao invés de AT
	U = (dMatrix*) malloc (N_t*sizeof(*U));
	Uold = (dMatrix*) malloc (N_t*sizeof(*Uold));
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
			for(unsigned i = 1; i<MATRIX_ORDER-1; i++)
				for(unsigned j = 1; j<MATRIX_ORDER-1; j++)
					U[k+1](i,j) = r*Uold[k+1](i,j)  + (1-r)*alpha*( K*(Uold[k+1](i+1,j) + Uold[k+1](i-1,j) + Uold[k+1](i,j+1) + Uold[k+1](i,j-1))/(dx*dx) + U[k](i,j)/dt);
				ER = U[k+1]-Uold[k+1];
				if(ER.norm() < 0.0000001) STOP = true;
				Uold[k+1] = U[k+1];
		}while(!STOP);
		element(k+1)= U[k+1](5,5);
		STOP=false;
	}
	/*------------------*/
	dMatrix *Uj;
	Uj = SolveHeatEquationI(conditions, K, N_t, L, t_max);
	cout << Uj[1]-U[1] << endl << endl;
	cout << Uj[2]-U[2] << endl;
	/*-----File operations------*/
    fileName(file); /*Gives a string to file according to N*/
	save(file, tStart, &start, element, "u_xx+u_yy=f(x,y)"); /*Save file in .m file to plot graph*/
	/*--------------------------*/
	return(0);
}
