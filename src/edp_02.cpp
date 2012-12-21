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

	/*------Time Operations-----*/
	time_t start;
	time(&start);
	clock_t tStart = clock();
	/*--------------------------*/

	/*Solving Operations*/
	MatrixXd g =  MatrixXd::Constant(N,N,0.0); 	/*Build g function*/
	MatrixXd U = dc_U(U,1,2,3,4); 	/*Build Dirichlet contour condition on U*/
	U = PoissonSparse(g, U); /*Solve L(U) = g using sparse matrix*/
	/*------------------*/

	/*-----File operations------*/
    char file[10];
    fileName(file); /*Gives a string to file according to N*/
	save(file, tStart, &start, U, "u_xx+u_yy = f(x,y)"); /*Save file in .m file to plot graph*/
	/*--------------------------*/

	return(0);
}
