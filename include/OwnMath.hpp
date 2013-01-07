/**
 * @file OwnMath.hpp
 * @author  Ataias Pereira Reis <ataiasreis@gmail.com>
 * @version 1.0.1
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 * A classe Time representa um momento de tempo.
 */

#include<stdheader.hpp>



dMatrix dReadConditions(char strFile[])
{
	/**
	 * The method is made to read a matrix of doubles from an ASCII file
	 * It requires the number of columns and lines
	 * It will check if the read matrix is of right size
	 *
	 * @param  strFile[] is the name of text file
	 * @return dMatrixRead is returned, with the contents of the file
	 */
	ifstream mFile;
	mFile.open(strFile);

	/*Loop to discover number of rows and columns while reading file*/
	int nRowNumber = 0;
	string strLineAux;
	while(getline(mFile, strLineAux)){
		nRowNumber++;
		}
	mFile.close();

	mFile.open(strFile);
	int nTotalElements(0);
	while(!mFile.eof()){
		double dAux;
		mFile >> dAux;
		nTotalElements++;
	}
	int nColumnNumber = nTotalElements/nRowNumber;
	mFile.close();
	/*Verify if the file has the right size to be worked*/
	dMatrix dMatrixRead;
	if((nRowNumber == MATRIX_ORDER) && (nColumnNumber == MATRIX_ORDER)){
	mFile.open(strFile);
	for (int i = 0; i < nRowNumber; i++) {
		for (int j = 0; j < nColumnNumber; j++) {
			mFile >> dMatrixRead(i,j);
	    }
	  }
	mFile.close();
	} else {
		 cout << "Error with size of matrix." << endl;
		 exit (EXIT_FAILURE);
	}
	return(dMatrixRead);
}

dMatrix PoissonDirichlet(dMatrix dNonHomogeneity,
					  dMatrix dBoundaryConditions
					   )
{
	/**
	 * Method to solve the poisson equation by a sparse method
	 *
	 * @param dNonHomogeneity is the non homogeneity in the right-hand side of poisson equation
	 * @param dBoundaryConditions contains the boundary and initial conditions
	 * @return dPoissonSolution is returned
	 */
	std::vector<T> tripletList; 	//triplet, needed to fill sparse matrix
	tripletList.reserve(PENTADIAGONAL*nMatrixSystemOrder);
	double dDeltaX = 1.0/(MATRIX_ORDER-1); 			/*Mesh has domain 0 < x < 1 and 0 < y < 1*/
	dBiVector dAuxiliaryIndexes = dBiVector::Zero();
	dVector b;

	/*Starting to build A*/
	long unsigned i,j;
	for(i = 0; i<nMatrixSystemOrder; i++)
		tripletList.push_back(T(i,i,-4));
	for(i = 0; i<= nMatrixSystemOrder-2; i++)
		if( (i+1) % (MATRIX_ORDER-2) != 0) tripletList.push_back(T(i,i+1,1));
	for(i = 0; i<= nMatrixSystemOrder-2; i++)
	    if( (i+1) % (MATRIX_ORDER-2) != 0) tripletList.push_back(T(i+1,i,1));
	for(i = 0; i < (nMatrixSystemOrder-(MATRIX_ORDER-2)); i++)
		tripletList.push_back(T(i,i+MATRIX_ORDER-2,1));
	for(i = 0; i<(nMatrixSystemOrder-(MATRIX_ORDER-2)); i++)
		tripletList.push_back(T(i+MATRIX_ORDER-2,i,1));
	SpMat A(nMatrixSystemOrder,nMatrixSystemOrder);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	/*Ending A*/

	for(i = 1; i<= (MATRIX_ORDER-2); i++)
	    for(j=1; j<=(MATRIX_ORDER-2); j++){

	    	/*ind has the indexes of U in a matrix with two columns*/
	    	dAuxiliaryIndexes(j+(i-1)*(MATRIX_ORDER-2)-1,0)= j+1;
	    	dAuxiliaryIndexes(j+(i-1)*(MATRIX_ORDER-2)-1,1)= i+1;
	    }

	for(i = 0; i<nMatrixSystemOrder; i++){
	    double dAuxToBVector = 0;
	    if(dAuxiliaryIndexes(i,0)-1 == 1)
	    	dAuxToBVector=dAuxToBVector+dBoundaryConditions(dAuxiliaryIndexes(i,1)-1,0);
	    if(dAuxiliaryIndexes(i,0)+1 == MATRIX_ORDER)
	    	dAuxToBVector=dAuxToBVector+dBoundaryConditions(dAuxiliaryIndexes(i,1)-1,MATRIX_ORDER-1);
	    if(dAuxiliaryIndexes(i,1)-1 == 1)
	    	dAuxToBVector=dAuxToBVector+dBoundaryConditions(MATRIX_ORDER-1,dAuxiliaryIndexes(i,0)-1);
	    if(dAuxiliaryIndexes(i,1)+1 == MATRIX_ORDER)
	    	dAuxToBVector=dAuxToBVector+dBoundaryConditions(0,dAuxiliaryIndexes(i,0)-1);

	    b(i) = dDeltaX*dDeltaX*dNonHomogeneity(dAuxiliaryIndexes(i,0),dAuxiliaryIndexes(i,1))-dAuxToBVector;
	}

	Eigen::SimplicialCholesky<SpMat> chol(A);
	dVector u = chol.solve(b);

	dMatrix dPoissonSolution;
	for(i=2; i<=(MATRIX_ORDER-1); i++)
	    for(j=2; j<=(MATRIX_ORDER-1); j++)
	        dPoissonSolution(MATRIX_ORDER+1-i-1,j-1) = u((i-2)*(MATRIX_ORDER-2)+j-1-1);
	return(dPoissonSolution);
}

void fileName(char* file)
{
	/**
	 * Method for creating a file name according to matrix size
	 *
	 * Matrix size is not a parameter of function, it is set as a global variable
	 * @param strFile is a char pointer where the string will be saved
	 * @return No return value.
	 */
    file[0] = 'm';
    file[1] = 'e';
    file[2] = 's';
    file[3] = 'h';
    int m = MATRIX_ORDER/100;
    file[4] = '0' + m;
    m = MATRIX_ORDER - 100*m;
    file[5] = '0' + m/10;
    m = m%10;
    file[6] = '0' + m;
    file[7] = '.';
    file[8] = 'm';
    file[9] = '\0';
}
void save(char file[],
	      clock_t tStart,
	      time_t *inicio,
	      MatrixXd U,
	      string equation
	      )
{
	/**
	 * Method for saving a matrix in a file .m
	 * After save, it is possible to open easily with octave
	 * and print the graph for file
	 *
	 * @param  file[] name of file where matrix will be save
	 * @param  tStart has the value of cpu time in the start, in this method
	 * the "tFinal" is calculated and subtracted from tStart
	 * @param  *inicio is a pointer that indicates the time of start, and date
	 * @param U is the matrix that will be saved
	 * @param equation is a string that receives the equation being solved,
	 * just to show in the .m file. Ex: Laplace Equation - u_xx + u_yy = f(x,y)
	 * @return No return value.
	 */
  ofstream m_file;
  m_file.open(file);
  m_file.unsetf(ios::floatfield);
  m_file.precision(15);
  m_file << "%_ Start time: " << ctime(inicio) << endl;
  m_file << "%_ Solving equation " << equation << endl;
  m_file << "%_ Square Mesh" << endl;
  m_file << "%_ Domain in x and y: 0<x<1 e 0<y<1 " << endl;
  m_file << "%_ MATRIX_ORDER = " << MATRIX_ORDER << ". Spent time: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s." << endl;
  m_file << "MATRIX_ORDER = " << MATRIX_ORDER << ";" << endl;
  m_file << endl;
  m_file << "U = [" << endl;
  m_file << U << endl;
  m_file << "];" << endl;
  m_file << "X = 0:1/(MATRIX_ORDER-1):1;" << endl;
  m_file << "Y = 0:1/(MATRIX_ORDER-1):1;" << endl;
  m_file << "figure;" << endl;
  m_file << "mesh(X,Y,U);" << endl;
  m_file << "hold on" << endl;
  m_file << "title('Solution of L(U) = f(x,y) with ";
  m_file << MATRIX_ORDER << "x" << MATRIX_ORDER << " points')" << endl;
  m_file << "xlabel('x')" << endl;
  m_file << "ylabel('y')" << endl;
  m_file << "zlabel('u(x,y)');" << endl;
  m_file << "mkdir png;" << endl;
  m_file << "hold off;" << endl;
  int a = MATRIX_ORDER/100;
  int b = MATRIX_ORDER - 100*a;
  int c = b%10;
  m_file << "print -dpng 'png/mesh" << a << b/10 << c << ".png';" << endl;
  time_t fim = time(NULL);
  m_file << "%_ End of operation: " << ctime(&fim) << endl;
  m_file.close();
  cout << "Mesh "<< MATRIX_ORDER << "x" << MATRIX_ORDER <<". Succeed to write in file. Operation time: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s." << endl;
}


void saveTmEqn(char file[],
			   clock_t tStart,
			   time_t *inicio,
			   string equation,
			   MatrixXd *U,
			   int N
			   )
{
	/**
	 * Method for saving various matrices in a file .m
	 * After save, it is possible to open easily with octave
	 * and make a video of the matrices
	 *
	 * Not implemented yet
	 *
	 * @param  file[] name of file where matrix will be save
	 * @param  tStart has the value of cpu time in the start, in this method
	 * the "tFinal" is calculated and subtracted from tStart
	 * @param  *inicio is a pointer that indicates the time of start, and date
	 * @param equation is a string that receives the equation being solved,
	 * just to show in the .m file. Ex: Laplace Equation - u_xx + u_yy = f(x,y)
	 * @param *U is a pointer to the 3d matrix U, that will be saved
	 * @param N indicates how many two-dimensional matrices there are in matrix U
	 * @return No return value.
	 */
	/*
  ofstream m_file;
  m_file.open(file);
  m_file.unsetf(ios::floatfield);
  m_file.precision(15);
  m_file << "%_ Start time: " << ctime(inicio) << endl;
  m_file << "%_ Solving equation u_xx + u_yy = u_t" << endl;
  m_file << "%_ Square Mesh" << endl;
  m_file << "%_ Domain in x and y: 0<x<1 e 0<y<1 " << endl;
  m_file << "%_ N = " << N << ". Spent time: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s." << endl;
  m_file << "N = " << N << ";" << endl;
  m_file << endl;

  int frames = time / dt;
  for(int i=1; i<frames; i++){
	  m_file << "U(:,:," << i << ") = [" << endl;
	  m_file << U.block<> << endl;
  }
  m_file << "];" << endl;
  m_file << endl;

  m_file << "X = 0:(1/N):(1-(1/N)/8);" << endl;
  m_file << "Y = 0:(1/N):(1-(1/N)/8);" << endl;
  m_file << "figure;" << endl;
  m_file << "mesh(X,Y,U);" << endl;
  m_file << "hold on" << endl;
  m_file << "title('Solution of L(U) = f(x,y) with ";
  m_file << N << "x" << N << " points')" << endl;
  m_file << "xlabel('x')" << endl;
  m_file << "ylabel('y')" << endl;
  m_file << "zlabel('u(x,y)');" << endl;
  m_file << "mkdir png;" << endl;
  int a = N/100;
  int b = N - 100*a;
  int c = b%10;
  m_file << "print -dpng 'png/mesh" << a << b/10 << c << ".png';" << endl;
  time_t fim = time(NULL);
  m_file << "%_ End of operation: " << ctime(&fim) << endl;
  m_file.close();
  cout << "Succeed to write in file. Operation time: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s." << endl;
  */
}
