/*
 * MatrixFilesIO.hpp
 *
 *  Created on: Jan 1, 2013
 *      Author: ataias
 */

#ifndef MATRIXFILESIO_HPP_
#define MATRIXFILESIO_HPP_

#include<stdheader.hpp>

#define INITIALIZE_FILE_READ_BINARY_OR_ASCII 		\
		switch(bBinaryOrASCII){\
		case(BINARY):\
			m_ReadFromFile.open(strFileName.c_str(), std::ios::binary);\
			break;\
		case(ASCII):\
			m_ReadFromFile.open(strFileName.c_str());\
			break;\
		}

#define INITIALIZE_FILE_WRITE_BINARY_OR_ASCII		\
		switch(bBinaryOrASCII){\
		case(BINARY):\
			m_WriteToFile.open(strFileName.c_str(), std::ios::binary);\
			break;\
		case(ASCII):\
			m_WriteToFile.open(strFileName.c_str());\
			break;\
		}

#define ERROR_MESSAGE_WRONG_PARAMETERS \
		std::cerr << "Wrong parameters given, method can only execute ";\
		if(bReadOrWrite) std::cerr << "WRITE and ";\
		else std::cerr << "READ and ";\
		\
		if(bBinaryOrASCII) std::cerr << "ASCII operations.\n";\
		else std::cerr << "Binary operations.\n";

#define WARNING_WRITE_ASCII \
		std::cerr << "Warning, given matrix order is different than real.\n";\
		std::cerr << "This can cause problems to reading, but none for writing in an ASCII file.\n";

#define WRONG_PARAMETERS \
	(m_bReadOrWrite!=bReadOrWrite) || (m_bBinaryOrASCII!=bBinaryOrASCII)

#define CHECK_AND_CLEAN_POINTER_DOUBLE \
		if(m_pdConvertedMatrix!=NULL) delete m_pdConvertedMatrix;\
			m_pdConvertedMatrix = NULL;

#define CONVERT_MATRIX_FROM_EIGEN_TO_DOUBLE_POINTER \
		m_pdConvertedMatrix = new double[m_nMatrixOrder*m_nMatrixOrder];\
		for(int i=0; i<m_nMatrixOrder; i++){\
			for(int j=0; j<m_nMatrixOrder; j++){\
				m_pdConvertedMatrix[i*m_nMatrixOrder+j]=dMatrixToBeConverted(i,j);\
			}\
		}

#define CHECK_AND_CLOSE_OPENED_FILES \
	if(m_ReadFromFile.is_open()) m_ReadFromFile.close();\
	else if(m_WriteToFile.is_open()) m_WriteToFile.close();

#define READ_DOUBLE_AND_CONVERT_TO_EIGEN_MATRIX \
	m_pdConvertedMatrix = new double[m_nMatrixOrder*m_nMatrixOrder];\
	m_ReadFromFile.read((char*)m_pdConvertedMatrix, sizeof(double)*m_nMatrixOrder*m_nMatrixOrder);\
	m_dMatrix = Eigen::Map<Eigen::MatrixXd>(m_pdConvertedMatrix,m_nMatrixOrder,m_nMatrixOrder);

#define CHECK_WRITE_SYSTEM_BINARY_OR_ASCII \
		if(!m_bBinaryOrASCII) WriteBinary(dMatrixToBeWritten);\
		else WriteASCII(dMatrixToBeWritten);

#define CHECK_READ_SYSTEM_BINARY_OR_ASCII \
		if(!m_bBinaryOrASCII) ReadBinary();\
		else ReadASCII();

#define UPDATE_NUMBER_OF_MATRICES_READ_OR_WRITTEN \
			m_nNumberOfMatrices++;

class MatrixFilesIO {
private:
	int m_nMatrixOrder;
	std::ofstream m_WriteToFile;
	std::ifstream m_ReadFromFile;
	bool m_bReadOrWrite;
	bool m_bBinaryOrASCII;
	Eigen::MatrixXd m_dMatrix;
	double *m_pdConvertedMatrix;
	int m_nNumberOfMatrices;

	int ErrorMessageIfWrongParameters(bool bReadOrWrite, bool bBinaryOrASCII);
	void ReadBinary();
	void WriteBinary(Eigen::MatrixXd dMatrixToBeWritten);
	void ReadASCII();
	void WriteASCII(Eigen::MatrixXd dMatrixToBeWritten);
	void ConvertEigenToMatrix(Eigen::MatrixXd dMatrixToBeConverted);

public:
	MatrixFilesIO(int nMatrixOrder, bool bReadOrWrite, bool bBinaryOrASCII , std::string strFileName);
	virtual ~MatrixFilesIO();
	void OperationIO();
	void OperationIO(Eigen::MatrixXd dMatrixToBeWritten);
	void CreateAndSaveMFile();
	void ShowMatrix();
	Eigen::MatrixXd ReturnMatrix();
};

#endif /* MATRIXFILESIO_HPP_ */
