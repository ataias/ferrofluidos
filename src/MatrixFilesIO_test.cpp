/**
 * @file MatrixFilesIO_test.cpp
 * @author Ataias Pereira Reis
 * Created on: Jan 1, 2013
 * This file is an unit test for reading and writing matrices
 * This was well achieve for ASCII and binary files
 *
 * The class is not finished as another function of it
 * will be job to create and save a .m file.
 *

 */

#include<MatrixFilesIO.hpp>

int main(){

	std::string strFileName1 = "teste1.txt";
	std::string strFileName2 = "teste2";
	int nMatrixOrder = 6;
	MatrixFilesIO *cTextFile;
	/*Testes com arquivos ASCII*/

	cTextFile = new MatrixFilesIO(nMatrixOrder, WRITE, ASCII, strFileName1);
	cTextFile->OperationIO(Eigen::MatrixXd::Random(nMatrixOrder,nMatrixOrder));
	cTextFile->OperationIO(Eigen::MatrixXd::Constant(nMatrixOrder,nMatrixOrder,5.5));
	delete cTextFile;

	cTextFile = new MatrixFilesIO(nMatrixOrder, READ, ASCII, strFileName1);
	cTextFile->OperationIO();
	cTextFile->ShowMatrix();
	cTextFile->OperationIO();
	cTextFile->ShowMatrix();
	delete cTextFile;

	/*Testes com arquivos Binários*/

	cTextFile = new MatrixFilesIO(nMatrixOrder, WRITE, BINARY, strFileName2);
	Eigen::MatrixXd matriz = Eigen::MatrixXd::Random(nMatrixOrder,nMatrixOrder);
	matriz = matriz*2 + Eigen::MatrixXd::Constant(nMatrixOrder,nMatrixOrder,0.1);

	cTextFile->OperationIO(matriz);
	cTextFile->OperationIO(Eigen::MatrixXd::Constant(nMatrixOrder,nMatrixOrder,4.5));
	cTextFile->OperationIO(Eigen::MatrixXd::Constant(nMatrixOrder,nMatrixOrder,3.5));
	delete cTextFile;

	cTextFile = new MatrixFilesIO(nMatrixOrder, READ, BINARY, strFileName2);
	cTextFile->OperationIO();
	cTextFile->ShowMatrix();
	cTextFile->OperationIO();
	cTextFile->ShowMatrix();
	cTextFile->OperationIO();
	cTextFile->ShowMatrix();
	delete cTextFile;
	return(SUCCESS);
}
