/*
 * Poisson.hpp
 *
 *  Created on: Jan 1, 2013
 *      Author: ataias
 */


#ifndef POISSON_HPP_
#define POISSON_HPP_

#include <stdheader.hpp>
#include <boost/python/detail/wrap_python.hpp>
#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
using namespace boost::python;

#ifndef DIRICHLET
#define DIRICHLET 0
#endif /*DIRICHLET*/

#ifndef NEUMANN
#define NEUMANN 1
#endif /*NEUMANN*/

#define CORRECT_CORNERS_DIRICHLET \
	m_dDirichletSolution(0,0) = 0.5*(m_dDirichletSolution(0,1)+m_dDirichletSolution(1,0));\
	m_dDirichletSolution(m_nMatrixOrder-1,0) = 0.5*(m_dDirichletSolution(m_nMatrixOrder-2,0)+\
		m_dDirichletSolution(m_nMatrixOrder-1,1));\
	m_dDirichletSolution(m_nMatrixOrder-1,m_nMatrixOrder-1) = 0.5*(m_dDirichletSolution(m_nMatrixOrder-2,m_nMatrixOrder-1)+\
			m_dDirichletSolution(m_nMatrixOrder-1,m_nMatrixOrder-2));\
	m_dDirichletSolution(0,m_nMatrixOrder-1) = 0.5*(m_dDirichletSolution(0,m_nMatrixOrder-2)+\
			m_dDirichletSolution(m_nMatrixOrder-1,1));

#define CORRECT_CORNERS_NEUMANN \
	m_dNeumannSolution(0,0) = 0.5*(m_dNeumannSolution(0,1)+m_dNeumannSolution(1,0));\
	m_dNeumannSolution(m_nMatrixOrder-1,0) = 0.5*(m_dNeumannSolution(m_nMatrixOrder-2,0)+\
		m_dNeumannSolution(m_nMatrixOrder-1,1));\
	m_dNeumannSolution(m_nMatrixOrder-1,m_nMatrixOrder-1) = 0.5*(m_dNeumannSolution(m_nMatrixOrder-2,m_nMatrixOrder-1)+\
			m_dNeumannSolution(m_nMatrixOrder-1,m_nMatrixOrder-2));\
	m_dNeumannSolution(0,m_nMatrixOrder-1) = 0.5*(m_dNeumannSolution(0,m_nMatrixOrder-2)+\
			m_dNeumannSolution(m_nMatrixOrder-1,1));

#define POISSON_NOSPARSE_INTERNAL_POINTS \
	dPoissonNoSparse(i,j) = -0.25*dDeltaX2*m_dNonHomogeneity(i,j) +\
	0.25*(dPoissonNoSparse(i-1,j)+dPoissonNoSparse(i+1,j)+\
	dPoissonNoSparse(i,j-1)+dPoissonNoSparse(i,j+1));

#define POISSON_NOSPARSE_WEST_POINTS \
		dPoissonNoSparse(i,0) =   (2.*dDeltaX*m_dBoundaryConditions(i,0)\
								+4.*dPoissonNoSparse(i,1)\
								-dPoissonNoSparse(i,2))/3.;

#define POISSON_NOSPARSE_EAST_POINTS \
		dPoissonNoSparse(i,m_nMatrixOrder-1) =   (-2.*dDeltaX*m_dBoundaryConditions(i,m_nMatrixOrder-1)\
												+4.*dPoissonNoSparse(i,m_nMatrixOrder-2)\
												-dPoissonNoSparse(i,m_nMatrixOrder-3))/3.;

#define POISSON_NOSPARSE_NORTH_POINTS \
		dPoissonNoSparse(0,j) =   (2.*dDeltaX*m_dBoundaryConditions(0,j)\
								+4.*dPoissonNoSparse(1,j)\
								-dPoissonNoSparse(2,j))/3.;

#define POISSON_NOSPARSE_SOUTH_POINTS \
		dPoissonNoSparse(m_nMatrixOrder-1,j) = (-2.*dDeltaX*m_dBoundaryConditions(m_nMatrixOrder-1,j)\
											  +4.*dPoissonNoSparse(m_nMatrixOrder-2,j)\
											  -dPoissonNoSparse(m_nMatrixOrder-3,j))/3.;

#define NORTH_DERIVATIVE_POISSON \
		dError(i,j) =  \
				              (-1.5*dPoissonNoSparse(i,j)+\
				              2*dPoissonNoSparse(i+1,j)-\
				              0.5*dPoissonNoSparse(i+2,j))/dDeltaX\
				              +m_dBoundaryConditions(i,j);

#define SOUTH_DERIVATIVE_POISSON \
		dError(i,j) =  \
				              ( 1.5*dPoissonNoSparse(i,j)\
				              -2*dPoissonNoSparse(i-1,j)\
				              +0.5*dPoissonNoSparse(i-2,j))/dDeltaX\
				              +m_dBoundaryConditions(i,j);

#define WEST_DERIVATIVE_POISSON \
		dError(i,j) =  \
		              (-1.5*dPoissonNoSparse(i,j)\
		              +2*dPoissonNoSparse(i,j+1)\
		              -0.5*dPoissonNoSparse(i,j+2))/dDeltaX\
		              +m_dBoundaryConditions(i,j);

#define EAST_DERIVATIVE_POISSON\
		dError(i,j) =  \
				              ( 1.5*dPoissonNoSparse(i,j)\
				              -2*dPoissonNoSparse(i,j-1)\
				              +0.5*dPoissonNoSparse(i,j-2))/dDeltaX\
				              +m_dBoundaryConditions(i,j);

#define POISSON_EQUATION_INTERNAL_POINT \
		dError(i,j) = \
				              (dPoissonNoSparse(i-1,j)\
				              +dPoissonNoSparse(i+1,j)\
				              +dPoissonNoSparse(i,j-1)\
				              +dPoissonNoSparse(i,j+1)\
				              -4*dPoissonNoSparse(i,j))/dDeltaX2\
				              - m_dNonHomogeneity(i,j);

#define INTERNAL_POINT \
	(i>0) && (i<(m_nMatrixOrder-1)) && \
	(j>0) && (j<(m_nMatrixOrder-1))

#ifndef WEST_POINT
#define WEST_POINT \
		(j==0) && (i>0) &&\
		(i<(m_nMatrixOrder-1))
#endif /*WEST_POINT*/

#ifndef LEFT_POINT
#define LEFT_POINT WEST_POINT
#endif /*LEFT POINT*/

#ifndef EAST_POINT
#define EAST_POINT \
	(j==(m_nMatrixOrder-1)) && \
	(i>0) && (i<(m_nMatrixOrder-1))
#endif /*EAST_POINT*/

#ifndef RIGHT_POINT
#define RIGHT_POINT EAST_POINT
#endif /*RIGHT_POINT*/

#ifndef NORTH_POINT
#define NORTH_POINT \
	(i==0) && (j>0) &&\
	(j<(m_nMatrixOrder-1))
#endif /*NORTH_POINT*/

#ifndef TOP_POINT
#define TOP_POINT NORTH_POINT
#endif /*TOP_POINT*/

#ifndef SOUTH_POINT
#define SOUTH_POINT \
		(i==(m_nMatrixOrder-1)) && (j>0) \
		&& (j<(m_nMatrixOrder-1))
#endif /*SOUTH_POINT*/

#ifndef BOTTOM_POINT
#define BOTTOM_POINT SOUTH_POINT
#endif /*BOTTOM_POINT*/

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
	Poisson(Eigen::MatrixXd dBoundaryConditions_,
				 	Eigen::MatrixXd dNonHomogeneity_,
				 	bool bDirichletOrNeumann,
				 	bool bSparseOrNot
				 	);

	Poisson();

	int PoissonPython(PyObject* dBoundaryConditions,
					  PyObject* dNonHomogeneity,
					  bool bDirichletOrNeumann,
					  bool bSparseOrNot,
					  const int nMatrixOrder
					  );
	template<typename Derived>
	int PoissonInit(const Eigen::MatrixBase<Derived>& dBoundaryConditions_,
			 	 	const Eigen::MatrixBase<Derived>& dNonHomogeneity_,
			 	 	bool bDirichletOrNeumann,
			 	 	bool bSparseOrNot
			 	 	);

	virtual ~Poisson();
	void PoissonSolver();
	void saveSolution(PyObject* pyArraySolution);

};

#endif /* POISSON_HPP_ */
