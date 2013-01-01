/*
 * MatrixFilesIO.cpp
 *
 *  Created on: Jan 1, 2013
 *      Author: ataias
 */

#include<stdheader.hpp>
#include "../include/MatrixFilesIO.hpp"

MatrixFilesIO::MatrixFilesIO(int nMatrixOrder, bool bReadOrWrite, bool bBinaryOrASCII, std::string strFileName) {
	m_nMatrixOrder = nMatrixOrder;
	m_bReadOrWrite = bReadOrWrite;
	m_bBinaryOrASCII = bBinaryOrASCII;
	m_dMatrix = Eigen::MatrixXd::Zero(nMatrixOrder,nMatrixOrder);
	m_pdConvertedMatrix = NULL;
	m_nNumberOfMatrices = 0;

	switch(bReadOrWrite){
	case(READ): 	INITIALIZE_FILE_READ_BINARY_OR_ASCII 	break;
	case(WRITE):	INITIALIZE_FILE_WRITE_BINARY_OR_ASCII 	break;
	}
} /*MatrixFilesIO::MatrixFilesIO*/

MatrixFilesIO::~MatrixFilesIO() {
	CHECK_AND_CLOSE_OPENED_FILES
	CHECK_AND_CLEAN_POINTER_DOUBLE
}

void MatrixFilesIO::ReadBinary() {
	ErrorMessageIfWrongParameters(READ,BINARY);
	READ_DOUBLE_AND_CONVERT_TO_EIGEN_MATRIX
	CHECK_AND_CLEAN_POINTER_DOUBLE
	std::cout << "Successfuly read from binary file.\n";
}
void MatrixFilesIO::WriteBinary(Eigen::MatrixXd dMatrixToBeWritten){
	ErrorMessageIfWrongParameters(WRITE,BINARY);
	ConvertEigenToMatrix(dMatrixToBeWritten);
	m_WriteToFile.write((char *)m_pdConvertedMatrix, sizeof(double)*m_nMatrixOrder*m_nMatrixOrder);
	std::cout << "Successfuly written in binary file.\n";
}
void MatrixFilesIO::ReadASCII(){
	ErrorMessageIfWrongParameters(READ,ASCII);
	for (int i = 0; i < m_nMatrixOrder; i++) {
		for (int j = 0; j < m_nMatrixOrder; j++) {
			m_ReadFromFile >> m_dMatrix(i,j);
	    }
	  }
	std::cout << "Successfuly read from file in ASCII.\n";
}
void MatrixFilesIO::WriteASCII(Eigen::MatrixXd dMatrixToBeWritten){
	ErrorMessageIfWrongParameters(WRITE,ASCII);
	m_WriteToFile << dMatrixToBeWritten << std::endl;
	if(m_nMatrixOrder!=dMatrixToBeWritten.cols()){
		WARNING_WRITE_ASCII
	}
	std::cout << "Successfuly written in file in ASCII.\n";
}

void MatrixFilesIO::ShowMatrix(){
	std::cout << m_dMatrix << std::endl;
}

void MatrixFilesIO::OperationIO(){
		CHECK_READ_SYSTEM_BINARY_OR_ASCII
		UPDATE_NUMBER_OF_MATRICES_READ_OR_WRITTEN
}

void MatrixFilesIO::OperationIO(Eigen::MatrixXd dMatrixToBeWritten){
		CHECK_WRITE_SYSTEM_BINARY_OR_ASCII
		UPDATE_NUMBER_OF_MATRICES_READ_OR_WRITTEN
}

void MatrixFilesIO::ConvertEigenToMatrix(Eigen::MatrixXd dMatrixToBeConverted){
		CHECK_AND_CLEAN_POINTER_DOUBLE
		CONVERT_MATRIX_FROM_EIGEN_TO_DOUBLE_POINTER
}
int MatrixFilesIO::ErrorMessageIfWrongParameters(bool bReadOrWrite, bool bBinaryOrASCII){
	if(WRONG_PARAMETERS){
		ERROR_MESSAGE_WRONG_PARAMETERS
		return(EXIT_FAILURE);
	}
	return(SUCCESS);
}

void CreateAndSaveMFile(){
	//TODO Precisa implementar uma função para criar e salvar arquivos .m
	// Já existe uma no arquivo ownmath.hpp.
	// Tem te preparar para fazer vídeo se necessário
	// Criar mais de uma opção... ver quais as melhores pra colocar aqui
}
