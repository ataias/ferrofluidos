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

#ifndef IOSTREAM_H
#define IOSTREAM_H
#include<iostream>
#endif

#ifndef STRING_H
#define STRING_H
#include<string>
#endif

#ifndef CSTDLIB_H
#define CSTDLIB_H
#include<cstdlib>
#endif

#ifndef CTIME_H
#define CTIME_H
#include<ctime>
#endif

#ifndef OWNMATH_HPP
#define OWNMATH_HPP
#include<OwnMath.hpp>
#endif

using namespace std;
using namespace Eigen;

int main(){
	/*Starting to implement Navier Stokes over here*/
	/*Variables*/
	/*V=(u,v)*/ /*dVelocityXdirection*/
	MatrixATd *dVelocityX; 	/*term of velocity in x axis*/
	MatrixATd *dVelocityY; 	/*term of velocity in y axis*/
	/*V*=(u*,v*)*/
	MatrixATd *dVelocityXNoPressure; 	/*term of velocity* in x axis*/
	MatrixATd *dVelocityYNoPressure; 	/*term of velocity* in y axis*/
	/*f = (fx, fy)*/
	MatrixATd *dExternalForceX;
	MatrixATd *dExternalForceY;
	double dVelocityYAverage; 		/*Auxiliary variable in the calculation of u* */
	double dVelocityYSum; 		/*Auxiliary variable in the calculation of v* */
	double dVelocityXAverage;
	double dVelocityXSum;
	int nTime;

	/*Here, needs to implement malloc of the pointers above*/
	//------------------------
	//Implement here
	//-----------------------
	/*First part of Navier Stokes algorithm*/
	// Considering fixed time for now
	int i,j; /*i,j-> varíaveis índice*/
	double dx, dt; /*dDeltaX, dDeltaT-> space and time steps*/

	for(i=1; i<N-1; i++){
		for(j=1; j<N-1; j++){
	dVelocityXSum = dVelocityX(i+1,j)[nTime]+dVelocityX(i-1,j)[nTime]+dVelocityX(i,j+1)[nTime]+dVelocityX(i,j-1)[nTime];
	dVelocityYSum = dVelocityY(i+1,j)[nTime]+dVelocityY(i-1,j)[nTime]+dVelocityY(i,j+1)[nTime]+dVelocityY(i,j-1)[nTime];

	dVelocityXAverage = 0.25*(dVelocityX(i,j)[nTime]+dVelocityX(i-1,j)[nTime]+dVelocityX(i,j+1)[nTime]+dVelocityX(i-1,j-1)[nTime]);
	dVelocityYAverage = 0.25*(dVelocityY(i,j)[nTime]+dVelocityY(i-1,j)[nTime]+dVelocityY(i,j+1)[nTime]+dVelocityY(i-1,j-1)[nTime]);

	double mi; /*mi-> dynamic viscosity coefficient*/
	double rho; /*rho-> fluid density*/
	/*nu = mi/rho = kinematic viscosity*/
	dVelocityXNoPressure = ((mi/rho)*(dVelocityXSum-4*dVelocityX(i,j)[nTime])/(dx*dx)+dExternalForceX(i,j)/rho
		 -dVelocityX(i,j)*(dVelocityX(i+1,j)[nTime]-dVelocityX(i-1,j)[nTime])/(2*dx)
		 -dVelocityYAverage*(dVelocityX(i,j+1)[nTime]-dVelocityX(i,j-1)[nTime])/(2*dx))*dt
	     +dVelocityX(i,j)[nTime];

	dVelocityYNoPressure = ((mi/rho)*(dVelocityYSum-4*dVelocityY(i,j)[nTime])/(dx*dx)+dExternalForceY(i,j)/rho
		 -dVelocityXAverage*(dVelocityY(i+1,j)[nTime]-dVelocityY(i-1,j)[nTime])/(2*dx)
		 -dVelocityY(i,j)[nTime]*(dVelocityY(i,j+1)[nTime]-dVelocityY(i,j-1)[nTime])/(2*dx))*dt
	     +dVelocityY(i,j)[nTime];
		}
	}
	/*Second Part of Calculations*/
	return(0);
}
