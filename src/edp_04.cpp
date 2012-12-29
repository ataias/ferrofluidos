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
	dMatrix conditions;
	char file[] = "file4.dat";
	conditions = dReadConditions(file);
	/*------------------------------------------------*/

	dMatrix *U;
	int N_t = 1000;
	U = SolveHeatEquationI(conditions, 1, N_t, 1, 0.1);
	/*------------------*/

	/*-----File operations------*/
    fileName(file); /*Gives a string to file according to N*/
/*    save("h0mesh050.m", tStart, &start, U[0], "u_xx+u_yy=f(x,y)");
    save("h01mesh050.m", tStart, &start, U[1], "u_xx+u_yy=f(x,y)");
    save("h02mesh050.m", tStart, &start, U[2], "u_xx+u_yy=f(x,y)");
    save("h03mesh050.m", tStart, &start, U[3], "u_xx+u_yy=f(x,y)");
    save("h04mesh050.m", tStart, &start, U[4], "u_xx+u_yy=f(x,y)");
    save("h05mesh050.m", tStart, &start, U[5], "u_xx+u_yy=f(x,y)");
    save("h10mesh050.m", tStart, &start, U[10], "u_xx+u_yy=f(x,y)");
    save("h15mesh050.m", tStart, &start, U[15], "u_xx+u_yy=f(x,y)");
    save("h20mesh050.m", tStart, &start, U[20], "u_xx+u_yy=f(x,y)");
    save("hFinalmesh050.m", tStart, &start, U[N_t-1], "u_xx+u_yy=f(x,y)");*/
	save(file, tStart, &start, U[N_t-1], "u_xx+u_yy=f(x,y)"); /*Save file in .m file to plot graph*/
	/*--------------------------*/
	return(0);
}
