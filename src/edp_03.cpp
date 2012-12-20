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

int main()
{
	/*This programs */
	time_t start;
	clock_t tStart = clock();
	N = 50;
	char file[10];
	do{
	/*------Time Operations-----*/
	time(&start);
	tStart = clock();
	/*--------------------------*/

	/*Solving Operations*/
	MatrixXd g =  MatrixXd::Constant(N,N,0.0); 	/*Build g function*/
	MatrixXd U = dc_U(U,1,2,3,4); 	/*Build Dirichlet contour condition on U*/
	U = PoissonSparse(g, U); /*Solve L(U) = g using sparse matrix*/
	/*------------------*/

	/*-----File operations------*/
    fileName(file); /*Gives a string to file according to N*/
	save(file, tStart, &start, U, "u_xx+u_yy=f(x,y)"); /*Save file in .m file to plot graph*/
	/*--------------------------*/
	N++;
	}while(N<=200);
	return(0);
}
