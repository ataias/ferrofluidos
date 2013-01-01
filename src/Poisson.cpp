/*
 * Poisson.cpp
 *
 *  Created on: Jan 1, 2013
 *      Author: ataias
 */

#include "../include/Poisson.hpp"

Poisson::Poisson(Eigen::MatrixXd dBoundaryConditions,
				 Eigen::MatrixXd dNonHomogeneity
				 ) {
	/**
		 * Constructor to the Poisson Class
		 * If the matrix of NonHomogeneity or Boundary Conditions is not given
		 * They will be assumed as zero
		 *
		 * @param dNonHomogeneity is the non homogeneity in the right-hand side of poisson equation
		 * @param dBoundaryConditions contains the boundary and initial conditions
		 * @return dPoissonSolution is returned
		 */

		bool bCheckSquareBoundaryConditions = dBoundaryConditions.rows() == dBoundaryConditions.cols();

		bool bCheckSquareNonHomogeneity = dNonHomogeneity.rows() == dNonHomogeneity.cols();

		bool bCheckSameOrderRow = dBoundaryConditions.rows()== dNonHomogeneity.rows();

		bool bCheckSameOrderColumn = dBoundaryConditions.cols() == dNonHomogeneity.cols();

		bool bCompatibleMatrices = bCheckSquareBoundaryConditions &&
				                   bCheckSquareNonHomogeneity &&
				                   bCheckSameOrderRow &&
				                   bCheckSameOrderColumn;
		if(!bCompatibleMatrices)
		{
			exit(EXIT_FAILURE);
		}

		/*It could have been used rows() or columns() of any of the matrix of parameters*/
		m_nMatrixOrder = dNonHomogeneity.rows();

		m_dNonHomogeneity = dNonHomogeneity;
		m_dBoundaryConditions = dBoundaryConditions;

}

Poisson::~Poisson() {
	// TODO Auto-generated destructor stub
}

Eigen::MatrixXd Poisson::PoissonDirichlet()
{
	std::vector<T> tripletList; 	//triplet, needed to fill sparse matrix
	int nMatrixSystemOrder = (m_nMatrixOrder-2)*(m_nMatrixOrder-2);
	tripletList.reserve(PENTADIAGONAL*nMatrixSystemOrder);
	double dDeltaX = 1.0/(m_nMatrixOrder-1); 			/*Mesh has domain 0 < x < 1 and 0 < y < 1*/
	Eigen::MatrixXd dAuxiliaryIndexes = Eigen::MatrixXd::Zero(nMatrixSystemOrder,2);;
	Eigen::VectorXd b(nMatrixSystemOrder);

	/*Starting to build A*/
	long unsigned i,j;
	for(i = 0; i<nMatrixSystemOrder; i++)
		tripletList.push_back(T(i,i,-4));
	for(i = 0; i<= nMatrixSystemOrder-2; i++)
		if( (i+1) % (m_nMatrixOrder-2) != 0) tripletList.push_back(T(i,i+1,1));
	for(i = 0; i<= nMatrixSystemOrder-2; i++)
	    if( (i+1) % (m_nMatrixOrder-2) != 0) tripletList.push_back(T(i+1,i,1));
	for(i = 0; i < (nMatrixSystemOrder-(m_nMatrixOrder-2)); i++)
		tripletList.push_back(T(i,i+m_nMatrixOrder-2,1));
	for(i = 0; i<(nMatrixSystemOrder-(m_nMatrixOrder-2)); i++)
		tripletList.push_back(T(i+m_nMatrixOrder-2,i,1));

	SpMat A(nMatrixSystemOrder,nMatrixSystemOrder);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	/*Ending A*/

	for(i = 1; i<= (m_nMatrixOrder-2); i++)
	    for(j=1; j<=(m_nMatrixOrder-2); j++){

	    	/*ind has the indexes of U in a matrix with two columns*/
	    	dAuxiliaryIndexes(j+(i-1)*(m_nMatrixOrder-2)-1,0)= j+1;
	    	dAuxiliaryIndexes(j+(i-1)*(m_nMatrixOrder-2)-1,1)= i+1;
	    }

	for(i = 0; i<nMatrixSystemOrder; i++){
	    double dAuxToBVector = 0;
	    if(dAuxiliaryIndexes(i,0)-1 == 1)
	    	dAuxToBVector=dAuxToBVector+m_dBoundaryConditions(dAuxiliaryIndexes(i,1)-1,0);
	    if(dAuxiliaryIndexes(i,0)+1 == m_nMatrixOrder)
	    	dAuxToBVector=dAuxToBVector+m_dBoundaryConditions(dAuxiliaryIndexes(i,1)-1,m_nMatrixOrder-1);
	    if(dAuxiliaryIndexes(i,1)-1 == 1)
	    	dAuxToBVector=dAuxToBVector+m_dBoundaryConditions(m_nMatrixOrder-1,dAuxiliaryIndexes(i,0)-1);
	    if(dAuxiliaryIndexes(i,1)+1 == m_nMatrixOrder)
	    	dAuxToBVector=dAuxToBVector+m_dBoundaryConditions(0,dAuxiliaryIndexes(i,0)-1);

	    b(i) = dDeltaX*dDeltaX*m_dNonHomogeneity(dAuxiliaryIndexes(i,0),dAuxiliaryIndexes(i,1))-dAuxToBVector;
	}

	Eigen::SimplicialCholesky<SpMat> chol(A);
	Eigen::VectorXd u(nMatrixSystemOrder);
	u = chol.solve(b);
	Eigen::MatrixXd dPoissonSolution(m_nMatrixOrder,m_nMatrixOrder);
	for(i=2; i<=(m_nMatrixOrder-1); i++)
	    for(j=2; j<=(m_nMatrixOrder-1); j++)
	        dPoissonSolution(m_nMatrixOrder+1-i-1,j-1) = u((i-2)*(m_nMatrixOrder-2)+j-1-1);
	return(dPoissonSolution);
}

Eigen::MatrixXd Poisson::PoissonNeumann()
{
	std::vector<T> tripletList; 	//triplet, needed to fill sparse matrix
	int nMatrixSystemOrder = m_nMatrixOrder*m_nMatrixOrder;
	tripletList.reserve(PENTADIAGONAL*nMatrixSystemOrder);
	double dDeltaX = 1.0/(m_nMatrixOrder-1); 			/*Mesh has domain 0 < x < 1 and 0 < y < 1*/
	Eigen::MatrixXd dAuxiliaryIndexes = Eigen::MatrixXd::Zero(nMatrixSystemOrder,2);;
	Eigen::VectorXd b(nMatrixSystemOrder);

	/*Starting to build A*/
	long unsigned i,j;
	for(i = 0; i<nMatrixSystemOrder; i++)
		tripletList.push_back(T(i,i,-4));
	for(i = 0; i<= nMatrixSystemOrder-2; i++)
		if( (i+1) % (m_nMatrixOrder-2) != 0) tripletList.push_back(T(i,i+1,1));
	for(i = 0; i<= nMatrixSystemOrder-2; i++)
	    if( (i+1) % (m_nMatrixOrder-2) != 0) tripletList.push_back(T(i+1,i,1));
	for(i = 0; i < (nMatrixSystemOrder-(m_nMatrixOrder-2)); i++)
		tripletList.push_back(T(i,i+m_nMatrixOrder-2,1));
	for(i = 0; i<(nMatrixSystemOrder-(m_nMatrixOrder-2)); i++)
		tripletList.push_back(T(i+m_nMatrixOrder-2,i,1));

	SpMat A(nMatrixSystemOrder,nMatrixSystemOrder);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	/*Ending A*/

	for(i = 1; i<= (m_nMatrixOrder-2); i++)
	    for(j=1; j<=(m_nMatrixOrder-2); j++){
	    	/*ind has the indexes of U in a matrix with two columns*/
	    	dAuxiliaryIndexes(j+(i-1)*(m_nMatrixOrder-2)-1,0)= j+1;
	    	dAuxiliaryIndexes(j+(i-1)*(m_nMatrixOrder-2)-1,1)= i+1;
	    }

	for(i = 0; i<nMatrixSystemOrder; i++){
	    double dAuxToBVector = 0;
	    if(dAuxiliaryIndexes(i,0)-1 == 1)
	    	dAuxToBVector=dAuxToBVector+m_dBoundaryConditions(dAuxiliaryIndexes(i,1)-1,0);
	    if(dAuxiliaryIndexes(i,0)+1 == m_nMatrixOrder)
	    	dAuxToBVector=dAuxToBVector+m_dBoundaryConditions(dAuxiliaryIndexes(i,1)-1,m_nMatrixOrder-1);
	    if(dAuxiliaryIndexes(i,1)-1 == 1)
	    	dAuxToBVector=dAuxToBVector+m_dBoundaryConditions(m_nMatrixOrder-1,dAuxiliaryIndexes(i,0)-1);
	    if(dAuxiliaryIndexes(i,1)+1 == m_nMatrixOrder)
	    	dAuxToBVector=dAuxToBVector+m_dBoundaryConditions(0,dAuxiliaryIndexes(i,0)-1);

	    b(i) = dDeltaX*dDeltaX*m_dNonHomogeneity(dAuxiliaryIndexes(i,0),dAuxiliaryIndexes(i,1))-dAuxToBVector;
	}

	Eigen::SimplicialCholesky<SpMat> chol(A);
	Eigen::VectorXd u(nMatrixSystemOrder);
	u = chol.solve(b);
	Eigen::MatrixXd dPoissonSolution(m_nMatrixOrder,m_nMatrixOrder);
	for(i=2; i<=(m_nMatrixOrder-1); i++)
	    for(j=2; j<=(m_nMatrixOrder-1); j++)
	        dPoissonSolution(m_nMatrixOrder+1-i-1,j-1) = u((i-2)*(m_nMatrixOrder-2)+j-1-1);
	return(dPoissonSolution);
}
