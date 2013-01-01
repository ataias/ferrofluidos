/*
 * Poisson.hpp
 *
 *  Created on: Jan 1, 2013
 *      Author: ataias
 */

#include<stdheader.hpp>

#ifndef POISSON_HPP_
#define POISSON_HPP_

class Poisson {
private:
	int m_nMatrixOrder;
	Eigen::MatrixXd m_dNonHomogeneity;
	Eigen::MatrixXd m_dBoundaryConditions;

public:
	Poisson();
	virtual ~Poisson();
	Poisson(Eigen::MatrixXd dBoundaryConditions,
			Eigen::MatrixXd dNonHomogeneity
			);
	Eigen::MatrixXd PoissonDirichlet();
	Eigen::MatrixXd PoissonNeumann();
};

#endif /* POISSON_HPP_ */
