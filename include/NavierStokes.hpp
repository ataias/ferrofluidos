/**
 * @file NavierStokes.hpp
 * @author Ataias Pereira Reis
 *
 *  Created on: Dec 29, 2012
 *      Author: ataias
 */

#ifndef NAVIERSTOKES_HPP_
#define NAVIERSTOKES_HPP_

#include "stdheader.hpp"
#include <boost/python/detail/wrap_python.hpp>
#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#include <boost/python/module.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/def.hpp>
using namespace boost::python;

/** @defgroup group1 The First Group
 *  This is the first group
 *  @{
 */

#ifndef VELOCITY_X_SUM
#define VELOCITY_X_SUM /*!< VELOCITY_X_SUM
 * This part stands for the sum of velocities in the X axis.
 * @f[\textbf{v}^s=(u^s, v^s) @f]
 * The 's' stands for 'sum'
 * @f[ u_{ij}^s=u_{i+1j}+u_{i-1j}+u_{ij+1}+u_{ij-1}@f]*/ \
m_dVelocityX(i+1,j)+\
m_dVelocityX(i-1,j)\
+m_dVelocityX(i,j+1)+\
m_dVelocityX(i,j-1)
#endif /* VELOCITY_X_SUM */

#ifndef VELOCITY_Y_SUM
#define VELOCITY_Y_SUM /*!<This part stands for the sum of velocities in the Y axis.
 * @f[\textbf{v}^s=(u^s, v^s) @f]
 * The 's' stands for 'sum'
 * @f[ v_{ij}^s=v_{i+1j}+v_{i-1j}+v_{ij+1}+v_{ij-1}@f]*/ \
m_dVelocityY(i+1,j)+\
m_dVelocityY(i-1,j)\
+m_dVelocityY(i,j+1)+\
m_dVelocityY(i,j-1)
#endif /* VELOCITY_Y_SUM */

/**This part stands for the average of velocities in the X axis.
 * @f[\textbf{v}^t=(u^t, v^t) @f]
 * The 't' stands for nothing special
 * @f[ u_{ij}^t=\frac{1}{4}(u_{ij}+u_{i-1j}+u_{i-1j-1}+u_{ij+1})@f]
 * */
#ifndef VELOCITY_X_AVERAGE
#define VELOCITY_X_AVERAGE \
0.25*(m_dVelocityX(i,j)+\
m_dVelocityX(i+1,j)\
+m_dVelocityX(i+1,j-1)+\
m_dVelocityX(i,j-1))
#endif /* VELOCITY_X_AVERAGE */

/**This part stands for the average of velocities in the X axis.
 * @f[\textbf{v}^t=(u^t, v^t) @f]
 * The 't' stands for nothing special
 * @f[ v_{ij}^t=\frac{1}{4}(v_{ij}+v_{i-1j}+v_{i-1j-1}+v_{ij+1})@f]
 * */
#ifndef VELOCITY_Y_AVERAGE
#define VELOCITY_Y_AVERAGE \
0.25*(m_dVelocityY(i,j)+\
m_dVelocityY(i-1,j)\
+m_dVelocityY(i-1,j+1)+\
m_dVelocityY(i,j+1))
#endif /* VELOCITY_Y_AVERAGE */

/** This part stand for the velocity obtained by Navier Stokes equation without considering the pressure
 * @f[u_{ij}^{*}=\left[\frac{\mu}{\rho}\left(\frac{u_{ij}^s-4u_{ij}}{\Delta x^2}\right)+\frac{f_{x,ij}}{\rho}-
 u_{ij}\frac{u_{i+1j}-u_{i-1j}}{2\Delta x}-v_{ij}^t\frac{u_{ij+1}-u_{ij-1}}{2\Delta x}\right]\Delta t + u_{ij}@f]
 * */
#ifndef VELOCITY_X_NO_PRESSURE
#define VELOCITY_X_NO_PRESSURE \
(m_dNu*(dVelocityXSum-4*m_dVelocityX(i,j))/m_dDeltaX2+((m_dExternalForceX(i,j)+m_dExternalForceX(i-1,j))/(2.*m_dRho))\
-m_dVelocityX(i,j)*(m_dVelocityX(i+1,j)-m_dVelocityX(i-1,j))/(2*m_dDeltaX)\
-dVelocityYAverage*(m_dVelocityX(i,j+1)-m_dVelocityX(i,j-1))/(2*m_dDeltaX))*m_dDeltaT\
+m_dVelocityX(i,j)
#endif /* VELOCITY_X_NO_PRESSURE */

/** This part stand for the velocity obtained by Navier Stokes equation without considering the pressure
 * @f[v_{ij}^{*}=\left[\frac{\mu}{\rho}\left(\frac{v_{ij}^s-4v_{ij}}{\Delta x^2}\right)+
\frac{f_{y,ij}}{\rho}-u_{ij}^t\frac{v_{i+1j}-v_{i-1j}}{2\Delta x}-v_{ij}\frac{v_{ij+1}-v_{ij-1}}{2\Delta x}\right]\Delta t + v_{ij}@f]
 * */
#ifndef VELOCITY_Y_NO_PRESSURE
#define VELOCITY_Y_NO_PRESSURE \
(m_dNu*(dVelocityYSum-4*m_dVelocityY(i,j))/m_dDeltaX2+((m_dExternalForceY(i,j)+m_dExternalForceY(i,j-1))/(2.*m_dRho))\
-dVelocityXAverage*(m_dVelocityY(i+1,j)-m_dVelocityY(i-1,j))/(2*m_dDeltaX)\
-m_dVelocityY(i,j)*(m_dVelocityY(i,j+1)-m_dVelocityY(i,j-1))/(2*m_dDeltaX))*m_dDeltaT\
+m_dVelocityY(i,j)
#endif /* VELOCITY_Y_NO_PRESSURE */

/** @}*/

#ifndef NON_HOMOGENEITY_NAVIER_PRESSURE
#define NON_HOMOGENEITY_NAVIER_PRESSURE \
dPressureNonHomogeneity(i,j) =\
m_dRho*(\
(m_dVelocityXNoPressure(i+1,j)-m_dVelocityXNoPressure(i,j))/m_dDeltaX\
+(m_dVelocityYNoPressure(i,j+1)-m_dVelocityYNoPressure(i,j))/m_dDeltaX\
)/m_dDeltaT
#endif /*NON_HOMOGENEITY_NAVIER*/

#ifndef POISSON_NAVIER_LEFT_BOUNDARY_CONDITION
#define POISSON_NAVIER_LEFT_BOUNDARY_CONDITION \
int j = 0; \
dPressureBoundaryConditions(i,0) =	\
-m_dRho*m_dNu*(-5*m_dVelocityXNoPressure(i,j+1)\
+4*m_dVelocityXNoPressure(i,j+2)\
-m_dVelocityXNoPressure(i,j+3))/m_dDeltaX2 \
-m_dRho*(m_dExternalForceX(i,j)+m_dExternalForceX(i,j+1))/2.0
#endif /*NON_HOMOGENEITY_NAVIER_LEFT*/

#ifndef POISSON_NAVIER_RIGHT_BOUNDARY_CONDITION
#define POISSON_NAVIER_RIGHT_BOUNDARY_CONDITION \
int j = m_nMatrixOrder - 1; \
dPressureBoundaryConditions(i,m_nMatrixOrder-1)= \
-m_dRho*m_dNu*(-5*m_dVelocityXNoPressure(i,j-1)\
+4*m_dVelocityXNoPressure(i,j-2)\
-m_dVelocityXNoPressure(i,j-3))/m_dDeltaX2 \
-m_dRho*(m_dExternalForceX(i,j)+m_dExternalForceX(i,j-1))/2.0
#endif /*NON_HOMOGENEITY_NAVIER_RIGHT*/

#ifndef POISSON_NAVIER_TOP_BOUNDARY_CONDITION
#define POISSON_NAVIER_TOP_BOUNDARY_CONDITION \
int i = 0;\
dPressureBoundaryConditions(0,j)= \
-m_dRho*m_dNu*(-5*m_dVelocityXNoPressure(i+1,j)\
+4*m_dVelocityXNoPressure(i+2,j)\
-m_dVelocityXNoPressure(i+3,j))/m_dDeltaX2 \
-m_dRho*(m_dExternalForceX(i,j)+m_dExternalForceX(i+1,j))/2.0
#endif /*NON_HOMOGENEITY_NAVIER_TOP*/

#ifndef POISSON_NAVIER_BOTTOM_BOUNDARY_CONDITION
#define POISSON_NAVIER_BOTTOM_BOUNDARY_CONDITION \
int i = m_nMatrixOrder - 1; \
dPressureBoundaryConditions(i,j)= \
-m_dRho*m_dNu*(-5*m_dVelocityXNoPressure(i-1,j)\
+4*m_dVelocityXNoPressure(i-2,j)\
-m_dVelocityXNoPressure(i-3,j))/m_dDeltaX2 \
-m_dRho*(m_dExternalForceX(i,j)+m_dExternalForceX(i-1,j))/2.0
#endif /*NON_HOMOGENEITY_NAVIER_RIGHT*/

#ifndef PRESSURE_INTERNAL_POINTS
#define PRESSURE_INTERNAL_POINTS \
dPressure(i,j) = -0.25*m_dDeltaX2*m_dPressureNonHomogeneity(i,j) +\
0.25*(dPressureOld(i-1,j)+dPressureOld(i+1,j)+\
dPressureOld(i,j-1)+dPressureOld(i,j+1))
#endif

/*nessa parte, eu estou na dúvida se tem sinais errados*/
#ifndef POISSON_NAVIER_LEFT
#define POISSON_NAVIER_LEFT \
dPressure(i,j) = dPressureOld(i,j+1) - m_dDeltaX*m_dPressureBoundaryConditions(i,j)
#endif

#ifndef POISSON_NAVIER_RIGHT
#define POISSON_NAVIER_RIGHT \
dPressure(i,j) = dPressureOld(i,j-1) - m_dDeltaX*m_dPressureBoundaryConditions(i,j)
#endif

#ifndef POISSON_NAVIER_TOP
#define POISSON_NAVIER_TOP \
dPressure(i,j) = dPressureOld(i+1,j) - m_dDeltaX*m_dPressureBoundaryConditions(i,j)
#endif

#ifndef POISSON_NAVIER_BOTTOM
#define POISSON_NAVIER_BOTTOM \
dPressure(i,j) = dPressureOld(i-1,j) - m_dDeltaX*m_dPressureBoundaryConditions(i,j)
#endif

#ifndef PRESSURE_NORTH_DERIVATIVE //assert i == 0
#define PRESSURE_NORTH_DERIVATIVE \
dError(i,j) = (\
m_dPressure(i+1,j)\
-m_dPressure(i,j))\
/m_dDeltaX\
+m_dPressureBoundaryConditions(i,j);
#endif

#ifndef PRESSURE_SOUTH_DERIVATIVE //assert i == m_nMatrixOrder-1
#define PRESSURE_SOUTH_DERIVATIVE \
dError(i,j) =  (\
m_dPressure(i,j)\
-m_dPressure(i-1,j))\
/m_dDeltaX\
-m_dPressureBoundaryConditions(i,j);
#endif

#ifndef PRESSURE_WEST_DERIVATIVE // assert j == 0
#define PRESSURE_WEST_DERIVATIVE \
dError(i,j) =  (\
m_dPressure(i,j)\
-m_dPressure(i,j+1))\
/m_dDeltaX\
+m_dPressureBoundaryConditions(i,j);
#endif

#ifndef PRESSURE_EAST_DERIVATIVE // assert j == m_nMatrixOrder - 1
#define PRESSURE_EAST_DERIVATIVE \
dError(i,j) = (\
m_dPressure(i,j)\
-m_dPressure(i,j-1))\
/m_dDeltaX\
-m_dPressureBoundaryConditions(i,j);
#endif

#ifndef PRESSURE_EQUATION_INTERNAL_POINT
#define PRESSURE_EQUATION_INTERNAL_POINT \
dError(i,j) = \
(m_dPressure(i-1,j)\
+m_dPressure(i+1,j)\
+m_dPressure(i,j-1)\
+m_dPressure(i,j+1)\
-4*m_dPressure(i,j))/m_dDeltaX2\
-m_dPressureNonHomogeneity(i,j);
#endif

#ifndef CORRECT_CORNERS_PRESSURE
#define CORRECT_CORNERS_PRESSURE \
m_dPressure(0,0) = 0.5*(m_dPressure(0,1)+m_dPressure(1,0));\
m_dPressure(m_nMatrixOrder-1,0) = 0.5*(m_dPressure(m_nMatrixOrder-2,0)+\
m_dPressure(m_nMatrixOrder-1,1));\
\
m_dPressure(m_nMatrixOrder-1,m_nMatrixOrder-1) = 0.5*(m_dPressure(m_nMatrixOrder-2,m_nMatrixOrder-1)+\
m_dPressure(m_nMatrixOrder-1,m_nMatrixOrder-2));\
\
m_dPressure(0,m_nMatrixOrder-1) = 0.5*(m_dPressure(0,m_nMatrixOrder-2)+\
m_dPressure(m_nMatrixOrder-1,1));
#endif

#ifndef PRESSURE_DIFF_X
#define PRESSURE_DIFF_X (m_dPressure(i+1,j)-m_dPressure(i-1,j))/(2*m_dDeltaX)
#endif

#ifndef PRESSURE_DIFF_Y
#define PRESSURE_DIFF_Y (m_dPressure(i,j+1)-m_dPressure(i,j-1))/(2*m_dDeltaX)
#endif

#ifndef ITERATION_LIMIT
#define ITERATION_LIMIT 30000
#endif

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


class NavierStokes {
private:

    Eigen::MatrixXd m_dVelocityXNoPressure; 	/*term of velocity* in x axis*/
    Eigen::MatrixXd m_dVelocityYNoPressure; 	/*term of velocity* in y axis*/
    /*f = (fx, fy)*/
    Eigen::MatrixXd m_dExternalForceX;
    Eigen::MatrixXd m_dExternalForceY;

    Eigen::MatrixXd m_dPressureNonHomogeneity;
    Eigen::MatrixXd m_dPressureBoundaryConditions;

    Eigen::MatrixXd m_dVelocityX;
    Eigen::MatrixXd m_dVelocityY;

    Eigen::MatrixXd m_dVelocityXNextStep;
    Eigen::MatrixXd m_dVelocityYNextStep;

    Eigen::MatrixXd m_dVelocityXBoundaryCondition;
    Eigen::MatrixXd m_dVelocityYBoundaryCondition;

    Eigen::MatrixXd m_dPressure;
    Eigen::MatrixXd m_dError;

    int m_nMatrixOrder;

    int m_nTime;
    double m_dDeltaX;
    double m_dDeltaX2;
    double m_dDeltaT;

    /*nu = mi/rho = kinematic viscosity*/
    double m_dMi; /*mi-> dynamic viscosity coefficient*/
    double m_dRho; /*rho-> fluid density*/
    double m_dNu;

    void VelocityNoPressure();
    void MakePressureConditions();
    void PressureSolver();
    void compute_dError();
    void VelocityNextStep();

public:
    void NavierStokesSolver();
    int NavierStokesPython(
        PyObject* dVelocityXBoundaryCondition,
        PyObject* dVelocityYBoundaryCondition,
        PyObject* dExternalForceX,
        PyObject* dExternalForceY,
        double dMi, double dRho,
        const int nMatrixOrder,
        const double dDeltaT
    );

    template<typename Derived>
    int NavierStokesInit(
        const Eigen::MatrixBase<Derived>& dVelocityXBoundaryCondition_,
        const Eigen::MatrixBase<Derived>& dVelocityYBoundaryCondition_,
        const Eigen::MatrixBase<Derived>& dExternalForceX_,
        const Eigen::MatrixBase<Derived>& dExternalForceY_,
        double dMi, double dRho, double dDeltaT
    );
//    virtual ~NavierStokes();

    void move(PyObject* pyArraySolutionX, PyObject* pyArraySolutionY); //function to move to python variable
    void NextStep();
};

#endif /* NAVIERSTOKES_HPP_ */
