/*
 * OwnMath.cpp
 *
 *  Created on: Dec 29, 2012
 *      Author: ataias
 */

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

#include<stdheader.hpp>
#include<NavierStokes.hpp>

class NavierStokes
{

private:
		dMatrix *dVelocityX;
		dMatrix *dVelocityY;

		dMatrix *dVelocityXNoPressure; 	/*term of velocity* in x axis*/
		dMatrix *dVelocityYNoPressure; 	/*term of velocity* in y axis*/
		/*f = (fx, fy)*/
		dMatrix *dExternalForceX;
		dMatrix *dExternalForceY;

		int nTime;
		const double dDeltaX = 1.0/(MATRIX_ORDER-1);
		double dDeltaT;

		/*nu = mi/rho = kinematic viscosity*/
		double m_dMi; /*mi-> dynamic viscosity coefficient*/
		double m_dRho; /*rho-> fluid density*/

		void VelocityNoPressure()
		{
			double dVelocityYAverage(0); 		/*Auxiliary variable in the calculation of u* */
			double dVelocityYSum(0); 		/*Auxiliary variable in the calculation of v* */
			double dVelocityXAverage(0);
			double dVelocityXSum(0);

			for(int i=1; i<MATRIX_ORDER-1; i++)
			{
				for(int j=1; j<MATRIX_ORDER-1; j++)
				{

				/*These are additional variables to help the calculation of velocities*/
				dVelocityXSum = VELOCITY_X_SUM;
				dVelocityYSum = VELOCITY_Y_SUM;
				dVelocityXAverage = VELOCITY_X_AVERAGE;
				dVelocityYAverage = VELOCITY_Y_AVERAGE;

				/*m_dNu = m_dMi/m_dRho = Kinematic Viscosity*/
				double m_dNu = m_dMi/m_dRho;

				/*Calculation of Velocity in X and Y without considering the term of pressure*/
				dVelocityXNoPressure = VELOCITY_X_NO_PRESSURE;
				dVelocityYNoPressure = VELOCITY_Y_NO_PRESSURE;
				}
			}
		} /*VelocityNoPressure()*/

public:
		NavierStokes(double dMi=1.0, double dRho=1.0){
			m_dMi = dMi;
			m_dRho = dRho;
		}
};

