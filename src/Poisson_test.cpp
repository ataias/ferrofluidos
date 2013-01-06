/**
 * @file Poisson_test.cpp
 * @author Ataias Pereira Reis
 *  Created on: Jan 1, 2013
 */

#include "../include/Poisson.hpp"
#include "../include/MatrixFilesIO.hpp"

int main(){
	using namespace std;
	using namespace Eigen;

	cout << "\nEnter with matrix order:\n";
	int nMatrixOrder = 0;
	cin >> nMatrixOrder;
	cout << nMatrixOrder;

	cout << "\n\nEnter with name file for Boundary Conditions: \n";
	string strBoundaryConditions;
	cin >> strBoundaryConditions;
	cout << strBoundaryConditions;

	cout << "\n\nEnter with name file for Non Homogeneity: \n";
	string strNonHomogeneity;
	cin >> strNonHomogeneity;
	cout << strNonHomogeneity;

	cout << "\n\nEnter with name file for OUTPUT matrix: \n";
	string strOutput;
	cin >> strOutput;
	cout << strOutput;

	cout << "\n\nType 0 for BINARY or 1 for ASCII: \n";
	bool bBinaryOrASCII;
	cin >> bBinaryOrASCII;
	cout << bBinaryOrASCII;

	cout << "\n\nType 0 for Non-sparse and 1 for Sparse: \n";
	bool bSparseOrNot;
	cin >> bSparseOrNot;
	cout << bSparseOrNot;

	cout << "\n\nType 0 for DIRICHLET and 1 for NEUMANN: \n";
	bool bDirichletOrNeumann;
	cin >> bDirichletOrNeumann;
	cout << bDirichletOrNeumann << endl << endl;

	//Read matrix from files
	MatrixFilesIO *p_dBoundaryConditions = new MatrixFilesIO(nMatrixOrder, READ, bBinaryOrASCII, strBoundaryConditions);
	p_dBoundaryConditions->OperationIO();
	MatrixXd dBoundaryConditions = p_dBoundaryConditions->ReturnMatrix();
	delete p_dBoundaryConditions;

	MatrixFilesIO *p_dNonHomogeneity = new MatrixFilesIO(nMatrixOrder, READ, bBinaryOrASCII, strNonHomogeneity);
	p_dNonHomogeneity->OperationIO();
	MatrixXd dNonHomogeneity = p_dNonHomogeneity->ReturnMatrix();
	delete p_dNonHomogeneity;

	Poisson *equation = new Poisson(dBoundaryConditions,dNonHomogeneity,bDirichletOrNeumann,bSparseOrNot);
	equation->PoissonSolver();
	MatrixFilesIO *p_dSolution = new MatrixFilesIO(nMatrixOrder, WRITE, bBinaryOrASCII, strOutput);
	p_dSolution->OperationIO(equation->ReturnMatrix());
	delete p_dSolution;
	delete equation;

	return(0);
}
