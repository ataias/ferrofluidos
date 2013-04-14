/**
 * @file main.cpp
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
 * This file is the main program that makes calls of all other functions
 * in a project at University of Brasilia dealing with ferrofluids.
 * It makes use of Eigen library and implements some algorithms to solve
 * partial differential equations, including Navier-Stokes equation.
 * \f$\rho\left( \frac{\partial \textbf{v}}{\partial t}+
 \textbf{v}\cdot\nabla\textbf{v}\right)=
 -\nabla p+\mu\nabla^2\textbf{v}+\textbf{f}\f$
 */

#include <../include/NavierStokes.hpp>

#if WRAP_PYTHON_NS
template<typename Derived>
int NavierStokes::NavierStokesInit(
		const Eigen::MatrixBase<Derived>& dVelocityXBoundaryCondition_,
 	 	const Eigen::MatrixBase<Derived>& dVelocityYBoundaryCondition_,
 	 	const Eigen::MatrixBase<Derived>& dExternalForceX_,
 	 	const Eigen::MatrixBase<Derived>& dExternalForceY_,
 	 	double dMi, double dRho, double dDeltaT
 	 	)
{
#else
	NavierStokes::NavierStokes(
			Eigen::MatrixXd dVelocityXBoundaryCondition,
			Eigen::MatrixXd dVelocityYBoundaryCondition,
			Eigen::MatrixXd dExternalForceX,
			Eigen::MatrixXd dExternalForceY,
			double dMi, double dRho
			);
#endif

#if WRAP_PYTHON_NS
		Eigen::MatrixXd dVelocityXBoundaryCondition = dVelocityXBoundaryCondition_;
		Eigen::MatrixXd dVelocityYBoundaryCondition = dVelocityYBoundaryCondition_;
		Eigen::MatrixXd dExternalForceX = dExternalForceX_;
		Eigen::MatrixXd dExternalForceY = dExternalForceY_;
#endif
		bool bCheckSquareVelocityX = dVelocityXBoundaryCondition.rows() == dVelocityXBoundaryCondition.cols();
		bool bCheckSquareVelocityY = dVelocityYBoundaryCondition.rows() == dVelocityYBoundaryCondition.cols();
		bool bCheckSquareForceX    = dExternalForceX.rows() == dExternalForceX.cols();
		bool bCheckSquareForceY    = dExternalForceY.rows() == dExternalForceY.cols();

		bool bCheckSameOrderRow    = (dVelocityXBoundaryCondition.rows() == dExternalForceX.rows()) &&
									 (dVelocityXBoundaryCondition.rows() == dExternalForceY.rows()) &&
									 (dVelocityYBoundaryCondition.rows() == dExternalForceX.rows()) &&
									 (dVelocityYBoundaryCondition.rows() == dExternalForceY.rows());

		bool bCheckSameOrderColumn = (dVelocityXBoundaryCondition.cols() == dExternalForceX.cols()) &&
				 	 	 	 	 	 (dVelocityXBoundaryCondition.cols() == dExternalForceY.cols()) &&
				 	 	 	 	 	 (dVelocityYBoundaryCondition.cols() == dExternalForceX.cols()) &&
				 	 	 	 	 	 (dVelocityYBoundaryCondition.cols() == dExternalForceY.cols());

		bool bCompatibleMatrices   = bCheckSquareVelocityX &&
								     bCheckSquareVelocityY &&
								     bCheckSquareForceX &&
								     bCheckSquareForceY;
		if(!bCompatibleMatrices)
		{
			exit(EXIT_FAILURE);
		}

		/*It could have been used rows() or columns() of any of the matrix of parameters*/
		m_nMatrixOrder = dVelocityXBoundaryCondition.rows();
		m_dVelocityXBoundaryCondition = dVelocityXBoundaryCondition;
		m_dVelocityYBoundaryCondition = dVelocityYBoundaryCondition;
		m_dExternalForceX = dExternalForceX;
		m_dExternalForceY = dExternalForceY;
		m_dMi = dMi;
		m_dRho = dRho;
		m_dDeltaT = dDeltaT;
#if WRAP_PYTHON_NS
		return(0);
#endif
}

void NavierStokes::VelocityNoPressure()
	//Computes velocity ignoring pressure
	//only for time t+DeltaT
	{
/*	double dVelocityYAverage(0); 		//Auxiliary variable in the calculation of u*
	double dVelocityYSum(0); 		//Auxiliary variable in the calculation of v*
	double dVelocityXAverage(0);
	double dVelocityXSum(0);

	for(int i=1; i<m_nMatrixOrder-1; i++)
		{
			for(int j=1; j<m_nMatrixOrder-1; j++)
			{
			//These are additional variables to help the calculation of velocities
			dVelocityXSum = VELOCITY_X_SUM;
			dVelocityYSum = VELOCITY_Y_SUM;
			dVelocityXAverage = VELOCITY_X_AVERAGE;
			dVelocityYAverage = VELOCITY_Y_AVERAGE;
			//m_dNu = m_dMi/m_dRho = Kinematic Viscosity
			double m_dNu = m_dMi/m_dRho;

			//Calculation of Velocity in X and Y without considering the term of pressure
			dVelocityXNoPressure = VELOCITY_X_NO_PRESSURE;
			dVelocityYNoPressure = VELOCITY_Y_NO_PRESSURE;
			}
		}*/
	} /*VelocityNoPressure()*/

#if WRAP_PYTHON_NS
int NavierStokes::NavierStokesPython(
		PyObject* dVelocityXBoundaryCondition,
		PyObject* dVelocityYBoundaryCondition,
		PyObject* dExternalForceX,
		PyObject* dExternalForceY,
		double dMi, double dRho,
		const int nMatrixOrder,
		const double dDeltaT
		)
{
	/**
		 * NavierStokes::NavierStokesPython
		 * Esta função tem só como função mapear matrizes numéricas do numpy em MatrixXd da Eigen.
		 * Ela já está totalmente funcional.
		 */
	Eigen::Map<Eigen::MatrixXd> _dVelocityXBoundaryCondition((double *) PyArray_DATA(dVelocityXBoundaryCondition),nMatrixOrder,nMatrixOrder);
	Eigen::Map<Eigen::MatrixXd> _dVelocityYBoundaryCondition((double *) PyArray_DATA(dVelocityYBoundaryCondition),nMatrixOrder,nMatrixOrder);
	Eigen::Map<Eigen::MatrixXd> _dExternalForceX((double *) PyArray_DATA(dExternalForceX),nMatrixOrder,nMatrixOrder);
	Eigen::Map<Eigen::MatrixXd> _dExternalForceY((double *) PyArray_DATA(dExternalForceY),nMatrixOrder,nMatrixOrder);

    return(
    		NavierStokesInit(_dVelocityXBoundaryCondition,
    						_dVelocityYBoundaryCondition,
    						_dExternalForceX,
    						_dExternalForceY,
    						dMi, dRho, dDeltaT
    	                    )
    	);
}
#endif

void NavierStokes::NavierStokesSolver(){
	std::cout << "Em construção.\n";
}

#if WRAP_PYTHON_NS
void NavierStokes::move(PyObject* pyArraySolutionX, PyObject* pyArraySolutionY){
	Eigen::Map<Eigen::MatrixXd> _pyArraySolutionX((double *) PyArray_DATA(pyArraySolutionX),m_nMatrixOrder,m_nMatrixOrder);
	Eigen::Map<Eigen::MatrixXd> _pyArraySolutionY((double *) PyArray_DATA(pyArraySolutionY),m_nMatrixOrder,m_nMatrixOrder);
	_pyArraySolutionX = m_dNavierStokesSolutionX;
	_pyArraySolutionY = m_dNavierStokesSolutionY;
}
#endif

#if WRAP_PYTHON_NS
BOOST_PYTHON_MODULE(libnavierstokes)
{
	class_<NavierStokes>("NavierStokes")
		.def("set", &NavierStokes::NavierStokesPython,
				(arg("dVelocityXBoundaryCondition"),
				 arg("dVelocityYBoundaryCondition"),
				 arg("dExternalForceX"),
				 arg("dExternalForceY"),
				 arg("dMi"),
				 arg("dRho"),
				 arg("nMatrixOrder"),
				 arg("dDeltaT")),
				 "This function takes important values required to solve the Navier Stokes problem on a square.")
		.def("solve", &NavierStokes::NavierStokesSolver, "Solves the Poisson problem on a square, with previous given conditions.")
		.def("move", &NavierStokes::move, (arg("pyArraySolutionX"),arg("pyArraySolutionY")), "Saves the solution of NS problem in pyArraySolutionX and Y objects.")
	;
}
#endif
