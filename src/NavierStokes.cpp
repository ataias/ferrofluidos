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

#include "../include/NavierStokes.hpp"

template<typename Derived>
int NavierStokes::NavierStokesInit(
		const Eigen::MatrixBase<Derived>& dVelocityXBoundaryCondition_,
 	 	const Eigen::MatrixBase<Derived>& dVelocityYBoundaryCondition_,
 	 	const Eigen::MatrixBase<Derived>& dExternalForceX_,
 	 	const Eigen::MatrixBase<Derived>& dExternalForceY_,
 	 	double dMi, double dRho, double dDeltaT
 	 	)
{
		Eigen::MatrixXd dVelocityXBoundaryCondition = dVelocityXBoundaryCondition_;
		Eigen::MatrixXd dVelocityYBoundaryCondition = dVelocityYBoundaryCondition_;
		Eigen::MatrixXd dExternalForceX = dExternalForceX_;
		Eigen::MatrixXd dExternalForceY = dExternalForceY_;

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
								     bCheckSquareForceY &&
                                                                     bCheckSameOrderRow &&
                                                                     bCheckSameOrderColumn;
		if(!bCompatibleMatrices)
		{
			exit(EXIT_FAILURE);
		}

		/*It could have been used rows() or columns() of any of the matrix of parameters*/
		m_nMatrixOrder = dVelocityXBoundaryCondition.rows();
		m_dVelocityXBoundaryCondition = dVelocityXBoundaryCondition;
		m_dVelocityYBoundaryCondition = dVelocityYBoundaryCondition;
		m_dVelocityX = dVelocityXBoundaryCondition;
		m_dVelocityY = dVelocityYBoundaryCondition;
		m_dExternalForceX = dExternalForceX;
		m_dExternalForceY = dExternalForceY;
		m_dMi = dMi;
		m_dRho = dRho;
		m_dNu = dMi/dRho;
		m_dDeltaT = dDeltaT;
		m_nTime = 0;
		m_dDeltaX = 1.0/(m_nMatrixOrder-1);
		m_dDeltaX2 = m_dDeltaX*m_dDeltaX;
		//---------------------Ordem da matriz de pressão é uma a mais que de Navier stokes -----------------------
		const int n_ExtraOrder = 2;
		m_dPressureNonHomogeneity = Eigen::MatrixXd::Zero(m_nMatrixOrder+n_ExtraOrder,m_nMatrixOrder+n_ExtraOrder);
		m_dPressureBoundaryConditions = Eigen::MatrixXd::Zero(m_nMatrixOrder+n_ExtraOrder,m_nMatrixOrder+n_ExtraOrder);
		m_dPressure = Eigen::MatrixXd::Zero(m_nMatrixOrder+n_ExtraOrder,m_nMatrixOrder+n_ExtraOrder);
		m_dError = Eigen::MatrixXd::Zero(m_nMatrixOrder+n_ExtraOrder,m_nMatrixOrder+n_ExtraOrder);
		//---------------------------------------------end pressure-------------------------------------------------
	    m_dVelocityXNextStep = Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder);
	    m_dVelocityYNextStep = Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder);

		return(0);

} /*NavierStokesInit()*/

void NavierStokes::VelocityNoPressure()
	//Computes velocity ignoring pressure
	//only for time t+DeltaT
	{
	//Auxiliary variable in the calculation of u*
	double dVelocityYAverage(0);
	//Auxiliary variable in the calculation of v*
	double dVelocityYSum(0);
	double dVelocityXAverage(0);
	double dVelocityXSum(0);

	Eigen::MatrixXd dVelocityXNoPressure = m_dVelocityX;
	Eigen::MatrixXd dVelocityYNoPressure = m_dVelocityY;

	for(int i=1; i<m_nMatrixOrder-1; i++)
		{
			for(int j=1; j<m_nMatrixOrder-1; j++)
			{
			//These are additional variables to help the calculation of velocities
			dVelocityXSum = VELOCITY_X_SUM; /*!< VELOCITY_X_SUM
			 * This part stands for the sum of velocities in the X axis.
			 * @f[\textbf{v}^s=(u^s, v^s) @f]
			 * The 's' stands for 'sum'
			 * @f[ u_{ij}^s=u_{i+1j}+u_{i-1j}+u_{ij+1}+u_{ij-1}@f]*/

			dVelocityYSum = VELOCITY_Y_SUM; /*!<This part stands for the sum of velocities in the Y axis.
			 * @f[\textbf{v}^s=(u^s, v^s) @f]
			 * The 's' stands for 'sum'
			 * @f[ v_{ij}^s=v_{i+1j}+v_{i-1j}+v_{ij+1}+v_{ij-1}@f]*/
			dVelocityXAverage = VELOCITY_X_AVERAGE;
			dVelocityYAverage = VELOCITY_Y_AVERAGE;
			//m_dNu = m_dMi/m_dRho = Kinematic Viscosity
			double m_dNu = m_dMi/m_dRho;

			//Calculation of Velocity in X and Y without considering the term of pressure
			dVelocityXNoPressure(i,j) = VELOCITY_X_NO_PRESSURE;
			dVelocityYNoPressure(i,j) = VELOCITY_Y_NO_PRESSURE;
			}
		}
	m_dVelocityXNoPressure = dVelocityXNoPressure;
	m_dVelocityYNoPressure = dVelocityYNoPressure;
	std::cout << "VelocityNoPressure\n";
	} /*VelocityNoPressure()*/

void NavierStokes::MakePressureConditions(){
	const int n_ExtraOrder = 2;
    Eigen::MatrixXd dPressureNonHomogeneity = Eigen::MatrixXd::Zero(m_nMatrixOrder+n_ExtraOrder,m_nMatrixOrder+n_ExtraOrder);
    Eigen::MatrixXd dPressureBoundaryConditions = Eigen::MatrixXd::Zero(m_nMatrixOrder+n_ExtraOrder,m_nMatrixOrder+n_ExtraOrder);
        
//  INTERNAL POINTS - NON HOMOGENEITY
	for(int i=1; i<=m_nMatrixOrder; i++)
		for(int j=1; j<=m_nMatrixOrder; j++){
            NON_HOMOGENEITY_NAVIER_PRESSURE;
		}
        
//  LEFT BOUNDARY CONDITIONS
    for(int j=1; j<=m_nMatrixOrder; j++){
        POISSON_NAVIER_LEFT_BOUNDARY_CONDITION;
    }
//  RIGHT BOUNDARY CONDITIONS
    for(int i=1; i<m_nMatrixOrder-1; i++){
        POISSON_NAVIER_RIGHT_BOUNDARY_CONDITION;
    }
//  TOP BOUNDARY CONDITIONS
    for(int j=1; j<m_nMatrixOrder-1; j++){
        POISSON_NAVIER_TOP_BOUNDARY_CONDITION;
    }
//  BOTTOM BOUNDARY CONDITIONS
    for(int j=1; j<m_nMatrixOrder-1; j++){
        POISSON_NAVIER_BOTTOM_BOUNDARY_CONDITION;
    }

    m_dPressureNonHomogeneity = dPressureNonHomogeneity;
    m_dPressureBoundaryConditions = dPressureBoundaryConditions;
	std::cout << "MakePressureConditions\n";
} /*MakePressureConditions()*/

void NavierStokes::compute_dError(){
	Eigen::MatrixXd dError = Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder);
	//Derivative tests
	for(int i = 0; i<m_nMatrixOrder; i++)
		for(int j = 0; j<m_nMatrixOrder; j++){
			//QUINAS
			if(i==0 && j==0)
				continue;
			if(i==0 && j==m_nMatrixOrder-1)
				continue;
			if(i==m_nMatrixOrder-1 && j==0)
				continue;
			if(i==m_nMatrixOrder-1 && j==m_nMatrixOrder-1)
				continue;

			if(i==0){        //NORTH POINT
				PRESSURE_NORTH_DERIVATIVE
			} else if(i==m_nMatrixOrder-1){ //SOUTH POINT
				PRESSURE_SOUTH_DERIVATIVE
			} else if (j==0){      //WEST POINT
				PRESSURE_WEST_DERIVATIVE
			} else if(j==m_nMatrixOrder-1){ //EAST POINT
				PRESSURE_EAST_DERIVATIVE
			} else { //INTERNAL POINT
				PRESSURE_EQUATION_INTERNAL_POINT
			}
		}
	m_dError = dError;
	//End derivative tests
}
void NavierStokes::PressureSolver(){
	bool bStopCondition=false;
	Eigen::MatrixXd dPressure = Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder);
    Eigen::MatrixXd dPressureOld = Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder);
	Eigen::MatrixXd dError = Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder);
	int k=0;
	int nMatrixIterations=0;
	do{
      if (nMatrixIterations == ITERATION_LIMIT) bStopCondition=true;
	  nMatrixIterations++;

    for(int i = -1; i<=m_nMatrixOrder; i++)
        for(int j = -1; j<=m_nMatrixOrder; j++){
        	//Agora vai de -1 até N, por causa da malha escalonada!
            if(LEFT_POINT){
                POISSON_NAVIER_LEFT;
            } else if(RIGHT_POINT) {
                POISSON_NAVIER_RIGHT;
            } else if(TOP_POINT) {
                POISSON_NAVIER_TOP;			
            } else if(BOTTOM_POINT){
                POISSON_NAVIER_BOTTOM;				
            }
		}

    for(int i = 0; i<m_nMatrixOrder; i++)
        for(int j = 0; j<m_nMatrixOrder; j++){
        	if(INTERNAL_POINT){
        		PRESSURE_INTERNAL_POINTS;
        	 }
		}

    for(int i = -1; i<=m_nMatrixOrder; i++)
        for(int j = -1; j<=m_nMatrixOrder; j++){
            if(LEFT_POINT){
                POISSON_NAVIER_LEFT;
            } else if(RIGHT_POINT) {
                POISSON_NAVIER_RIGHT;
            } else if(TOP_POINT) {
                POISSON_NAVIER_TOP;
            } else if(BOTTOM_POINT){
                POISSON_NAVIER_BOTTOM;
            }
		}

		  double dMinimumValue = dPressure.minCoeff();
		  dPressure = dPressure - Eigen::MatrixXd::Constant(m_nMatrixOrder, m_nMatrixOrder, dMinimumValue);
          NavierStokes::compute_dError();
		  double d_error = m_dError.norm();
		  if(d_error < 1e-13) bStopCondition=true; //Condição de parada

		  dPressureOld = dPressure;
		  std::cout << "N = " << nMatrixIterations << " " << d_error << "\n";
	}while(!bStopCondition);

	//Corners
	m_dPressure = dPressure;
	CORRECT_CORNERS_PRESSURE
}

void NavierStokes::VelocityNextStep(){
	Eigen::MatrixXd dPressureDiffX = Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder);
	Eigen::MatrixXd dPressureDiffY = Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder);
	//Computing gradient of pressure
	//in X
	for(int i=1; i<m_nMatrixOrder-1; i++){
		for(int j=1; j<m_nMatrixOrder-1; j++){
			dPressureDiffX(i,j)=PRESSURE_DIFF_X;
			dPressureDiffY(i,j)=PRESSURE_DIFF_Y;
		}
	}
	m_dVelocityXNextStep = m_dVelocityXNoPressure - (m_dDeltaT / m_dRho)*dPressureDiffX;
	m_dVelocityYNextStep = m_dVelocityYNoPressure - (m_dDeltaT / m_dRho)*dPressureDiffY;
	//Compute next step

	std::cout << "VelocityNextStep\n";
}

void NavierStokes::NextStep(){
	m_dVelocityX = m_dVelocityXNextStep;
	m_dVelocityY = m_dVelocityYNextStep;
	std::cout << "NextStep\n";
}

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
    		NavierStokesInit((Eigen::MatrixXd)_dVelocityXBoundaryCondition,
    						 (Eigen::MatrixXd)_dVelocityYBoundaryCondition,
    						 (Eigen::MatrixXd)_dExternalForceX,
    						 (Eigen::MatrixXd)_dExternalForceY,
    						dMi, dRho, dDeltaT
    	                    )
    	);
}


void NavierStokes::NavierStokesSolver(){
	NavierStokes::VelocityNoPressure();
	NavierStokes::MakePressureConditions();
	NavierStokes::PressureSolver();
	NavierStokes::VelocityNextStep();
	std::cout << "Velocity in X given" << std::endl;
	std::cout << m_dVelocityX << "\n\n";
	std::cout << "Velocity in Y given" << std::endl;
	std::cout << m_dVelocityY << std::endl;
	std::cout << "Pressure Boundary Conditions" << "\n";
	std::cout << m_dPressureBoundaryConditions << "\n\n";
	std::cout << "Pressure Non-Homogeneity" << "\n";
    std::cout << m_dPressureNonHomogeneity << "\n\n";
    std::cout << "VelocityXNoPressure" << "\n";
    std::cout << m_dVelocityXNoPressure << "\n\n";
    std::cout << "VelocityYNoPressure" << "\n";
    std::cout << m_dVelocityYNoPressure << "\n\n";
	std::cout << "Em construção.\n";
}


void NavierStokes::move(PyObject* pyArraySolutionX, PyObject* pyArraySolutionY){
	Eigen::Map<Eigen::MatrixXd> _pyArraySolutionX((double *) PyArray_DATA(pyArraySolutionX),m_nMatrixOrder,m_nMatrixOrder);
	Eigen::Map<Eigen::MatrixXd> _pyArraySolutionY((double *) PyArray_DATA(pyArraySolutionY),m_nMatrixOrder,m_nMatrixOrder);
	_pyArraySolutionX = m_dVelocityXNextStep;
	_pyArraySolutionY = m_dVelocityYNextStep;
}

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
