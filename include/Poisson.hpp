/*
 * Poisson.hpp
 *
 *  Created on: Jan 1, 2013
 *      Author: ataias
 */

#include<stdheader.hpp>

#ifndef POISSON_HPP_
#define POISSON_HPP_

#define DIRICHLET 0
#define NEUMANN 1

#define CORRECT_CORNERS_DIRICHLET \
	m_dDirichletSolution(0,0) = 0.5*(m_dDirichletSolution(0,1)+m_dDirichletSolution(1,0));\
	m_dDirichletSolution(m_nMatrixOrder-1,0) = 0.5*(m_dDirichletSolution(m_nMatrixOrder-2,0)+\
		m_dDirichletSolution(m_nMatrixOrder-1,1));\
	m_dDirichletSolution(m_nMatrixOrder-1,m_nMatrixOrder-1) = 0.5*(m_dDirichletSolution(m_nMatrixOrder-2,m_nMatrixOrder-1)+\
			m_dDirichletSolution(m_nMatrixOrder-1,m_nMatrixOrder-2));\
	m_dDirichletSolution(0,m_nMatrixOrder-1) = 0.5*(m_dDirichletSolution(0,m_nMatrixOrder-2)+\
			m_dDirichletSolution(m_nMatrixOrder-1,1));

#define POISSON_NOSPARSE_INTERNAL_POINTS \
	dPoissonNoSparse(i,j) = -0.25*dDeltaX*dDeltaX*m_dNonHomogeneity(i,j) +\
	0.25*(dPoissonNoSparse(i-1,j)+dPoissonNoSparse(i+1,j)+\
	dPoissonNoSparse(i,j-1)+dPoissonNoSparse(i,j+1));

#define POISSON_NOSPARSE_WEST_POINTS \
	dPoissonNoSparse(i,0) = 0.25*(2*dPoissonNoSparse(i,1)+\
	2*dDeltaX*m_dBoundaryConditions(i,0)+dPoissonNoSparse(i+1,0)\
	+dPoissonNoSparse(i-1,0))\
	-0.25*dDeltaX*dDeltaX*m_dNonHomogeneity(i,0);

#define POISSON_NOSPARSE_EAST_POINTS \
	dPoissonNoSparse(i,m_nMatrixOrder-1) = 0.25*(2*dPoissonNoSparse(i,m_nMatrixOrder-2)+\
	2*dDeltaX*m_dBoundaryConditions(i,m_nMatrixOrder-1)+\
	dPoissonNoSparse(i+1,m_nMatrixOrder-1)+\
	dPoissonNoSparse(i-1,m_nMatrixOrder-1))-\
	0.25*dDeltaX*dDeltaX*m_dNonHomogeneity(i,m_nMatrixOrder-1);

#define POISSON_NOSPARSE_NORTH_POINTS \
	dPoissonNoSparse(0,j) = 0.25*(2*dPoissonNoSparse(1,j)+\
	2*dDeltaX*m_dBoundaryConditions(0,j)+\
	dPoissonNoSparse(0,j-1)+\
	dPoissonNoSparse(0,j+1))-\
	0.25*dDeltaX*dDeltaX*m_dNonHomogeneity(0,j);

#define POISSON_NOSPARSE_SOUTH_POINTS \
	dPoissonNoSparse(m_nMatrixOrder-1,j) = 0.25*(2*dPoissonNoSparse(m_nMatrixOrder-2,j)-\
	2*dDeltaX*m_dBoundaryConditions(m_nMatrixOrder-1,j)+\
	dPoissonNoSparse(m_nMatrixOrder-1,j-1)+\
	dPoissonNoSparse(m_nMatrixOrder-1,j+1))-\
	0.25*dDeltaX*dDeltaX*m_dNonHomogeneity(m_nMatrixOrder-1,j);

#define INTERNAL_POINT \
	(i>0) && (i<(m_nMatrixOrder-1)) && \
	(j>0) && (j<(m_nMatrixOrder-1))

#define WEST_POINT \
		(j==0) && (i>0) &&\
		(i<(m_nMatrixOrder-1))

#define EAST_POINT \
	(j==(m_nMatrixOrder-1)) && \
	(i>0) && (i<(m_nMatrixOrder-1))

#define NORTH_POINT \
	(i==0) && (j>0) &&\
	(j<(m_nMatrixOrder-1))

#define SOUTH_POINT \
		(i==(m_nMatrixOrder-1)) && (j>0) \
		&& (j<(m_nMatrixOrder-1))

class Poisson {
private:
	int m_nMatrixOrder;
	Eigen::MatrixXd m_dNonHomogeneity;
	Eigen::MatrixXd m_dBoundaryConditions;
	Eigen::MatrixXd m_dDirichletSolution;
	Eigen::MatrixXd m_dNeumannSolution;
	bool m_bDirichletOrNeumann;
	bool m_bCheckIfSolved;
	bool m_bSparseOrNot;
	void PoissonDirichletSparseSolver();
	void PoissonNeumannSparseSolver();
	void PoissonDirichletNoSparseSolver();
	void PoissonNeumannNoSparseSolver();
public:
	Poisson(Eigen::MatrixXd dBoundaryConditions,
			Eigen::MatrixXd dNonHomogeneity,
			bool bDirichletOrNeumann,
			bool bSparseOrNot
			);
	virtual ~Poisson();
	void PoissonSolver();
	Eigen::MatrixXd ReturnMatrix();
};

#endif /* POISSON_HPP_ */
