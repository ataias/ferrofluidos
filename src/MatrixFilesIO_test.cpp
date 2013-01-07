/**
 * @file MatrixFilesIO_test.cpp
 * @author Ataias Pereira Reis
 * Created on: Jan 1, 2013
 * Este arquivo é um teste unitário para leitura e escrita em matrizes.
 * Ele testa a capacidade de leitura e escrita em arquivos texto, codificados em
 * ASCII e também em arquivos binários contendo números em double.
 * A ideia de arquivos binários é interessante pois salvando os números em
 * binário pode economizar espaço se fosse salvar tudo em texto.
 *
 * A classe responsável pelo IO deste projeto ainda não está pronta.
 * Falta criar uma função para criar arquivos .m. Mas ainda está sendo
 * pensando se isto vai em outra classe ou nesta mesmo.
 * Poderia ser criada uma classe responsável pela criação de arquivos que vão
 * criar gráficos bi e tridimensionais.
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
