/**
 * @file MatrixFilesIO.cpp
 * @author Ataias Pereira Reis
 *  Created on: Jan 1, 2013
 */

#include "../include/stdheader.hpp"
#include "../include/MatrixFilesIO.hpp"

MatrixFilesIO::MatrixFilesIO(int nMatrixOrder, bool bReadOrWrite, bool bBinaryOrASCII, std::string strFileName) {
	/** Este é o construtor da classe de IO para matrizes, MatrixFilesIO. O construtor deve ser utilizado sempre.
	 * Uma vez criado um objeto, ele só poderá executar o método específico com o qual foi construído.
	 * Vamos falar dos parâmetros então.
	 * @param nMatrixOrder é um inteiro que pede a ordem da matriz que será lida ou salva.
	 * @param bReadOrWrite pode tomar os valores READ ou WRITE, que são 0 e 1 respectivamente.
	 * isso limite então se o objeto servirá para ler dados ou salvá-los.
	 * @param bBinaryOrASCII pode tomar os valores BINARY ou ASCII que indica se
	 * será salvo ou lido arquivos em binário double ou em caracteres ASCII.
	 * @param strFileName é o nome do arquivo que será criado ou lido dele.
	 * */
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
	/** Esta função mostra a matriz m_dMatrix que é uma matriz lida ou a última matriz salva em arquivo.
	 * */
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
	/** Este método tem um propósito bem específico. Converter uma matriz do formato Eigen para uma matriz em C++ de doubles.
	 * Com isto é possível salvar sem problemas no arquivo binário mesmo lidando com o formato MatrixXd.
	 * Formatos estáticos não precisam de conversão, mas este formato dinâmico apresentou muitos problemas. Daí a necessidade
	 * desta função.
	 * */
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

Eigen::MatrixXd MatrixFilesIO::ReturnMatrix(){
	return(m_dMatrix);
}
