/**
 * @file Poisson.cpp
 * @author Ataias Pereira Reis
 *  Created on: Jan 1, 2013
 */

#include "../include/Poisson.hpp"

#if WRAP_PYTHON
template<typename Derived>
int Poisson::PoissonInit(const Eigen::MatrixBase<Derived>& dBoundaryConditions_,
						 const Eigen::MatrixBase<Derived>& dNonHomogeneity_,
						 bool bDirichletOrNeumann,
						 bool bSparseOrNot)
{
#else
Poisson::Poisson(Eigen::MatrixXd dBoundaryConditions,
				 Eigen::MatrixXd dNonHomogeneity,
				 bool bDirichletOrNeumann,
				 bool bSparseOrNot
				 ) {
#endif
	/**
		 * Construtor da classe Poisson
		 * É obrigatória a entrada da matriz de condições de contorno,
		 * pode ter uma condição inicial na matriz, mas ela será desconsiderada.
		 * A matriz de não-homogeneidade é a matriz que segue a equação do lado direito.
		 * A equação em questão é a seguinte:
		 * \f$ \nabla^2 u = g \f$
		 * @param dNonHomogeneity é matriz de não homogeneidade, a função g
		 * @param dBoundaryConditions contém a condição inicial e de contorno
		 *
		 * É importante que as matrizes sejam do mesmo tamanho. Caso contrário o programa
		 * será finalizado, com erro.
		 * Outro detalhe, as matrizes devem ser quadradas!
		 */

#if WRAP_PYTHON
		Eigen::MatrixXd dBoundaryConditions = dBoundaryConditions_;
		Eigen::MatrixXd dNonHomogeneity = dNonHomogeneity_;
#endif
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
		m_bDirichletOrNeumann = bDirichletOrNeumann;
		m_dDirichletSolution = Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder);
		m_dNeumannSolution = Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder);
		m_bCheckIfSolved = false;
		m_bSparseOrNot = bSparseOrNot;
#if WRAP_PYTHON
		return(0);
#endif
}

Poisson::~Poisson() {
	// TODO Auto-generated destructor stub
}

void Poisson::PoissonDirichletSparseSolver()
{
	std::vector<T> tripletList; 	//triplet, needed to fill sparse matrix
	int nMatrixSystemOrder = (m_nMatrixOrder-2)*(m_nMatrixOrder-2);
	tripletList.reserve(PENTADIAGONAL*nMatrixSystemOrder);
	double dDeltaX = 1.0/(m_nMatrixOrder-1); 			/*Mesh has domain 0 < x < 1 and 0 < y < 1*/
	Eigen::MatrixXd dAuxiliaryIndexes = Eigen::MatrixXd::Zero(nMatrixSystemOrder,2);;
	Eigen::VectorXd b(nMatrixSystemOrder);

	/*Starting to build A*/
	long unsigned i,j;
	//Construir diagonal principal
	for(i = 0; i<nMatrixSystemOrder; i++)
		tripletList.push_back(T(i,i,-4));

	// Construir diagonal superior 1
	for(i = 0; i<= nMatrixSystemOrder-2; i++)
		if( (i+1) % (m_nMatrixOrder-2) != 0) tripletList.push_back(T(i,i+1,1));

	// Construir diagonal inferior 1
	for(i = 0; i<= nMatrixSystemOrder-2; i++)
	    if( (i+1) % (m_nMatrixOrder-2) != 0) tripletList.push_back(T(i+1,i,1));

	//Construir diagonal superior 2
	for(i = 0; i < (nMatrixSystemOrder-(m_nMatrixOrder-2)); i++)
		tripletList.push_back(T(i,i+m_nMatrixOrder-2,1));

	//Construir diagonal inferior 2
	for(i = 0; i<(nMatrixSystemOrder-(m_nMatrixOrder-2)); i++)
		tripletList.push_back(T(i+m_nMatrixOrder-2,i,1));

	SpMat A(nMatrixSystemOrder,nMatrixSystemOrder);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	/*Ending A*/

	//---------------------
	for(i = 1; i<= (m_nMatrixOrder-2); i++)
	    for(j=1; j<=(m_nMatrixOrder-2); j++){
	    	/*ind has the indexes of U in a matrix with two columns*/
	    	dAuxiliaryIndexes(j+(i-1)*(m_nMatrixOrder-2)-1,0)= j;
	    	dAuxiliaryIndexes(j+(i-1)*(m_nMatrixOrder-2)-1,1)= i;
	    }

	for(i = 0; i<nMatrixSystemOrder; i++){
	    double dAuxToBVector = 0;
	    if(dAuxiliaryIndexes(i,0) == 1)
	    	dAuxToBVector=dAuxToBVector+m_dBoundaryConditions(0,dAuxiliaryIndexes(i,1));
	    if(dAuxiliaryIndexes(i,0)+2 == m_nMatrixOrder)
	    	dAuxToBVector=dAuxToBVector+m_dBoundaryConditions(m_nMatrixOrder-1,dAuxiliaryIndexes(i,1));
	    if(dAuxiliaryIndexes(i,1) == 1)
	    	dAuxToBVector=dAuxToBVector+m_dBoundaryConditions(dAuxiliaryIndexes(i,0),0);
	    if(dAuxiliaryIndexes(i,1)+2 == m_nMatrixOrder)
	    	dAuxToBVector=dAuxToBVector+m_dBoundaryConditions(dAuxiliaryIndexes(i,0),m_nMatrixOrder-1);

	    b(i) = dDeltaX*dDeltaX*m_dNonHomogeneity(dAuxiliaryIndexes(i,0),dAuxiliaryIndexes(i,1))-dAuxToBVector;
	}
	//--------------------
	Eigen::SimplicialLDLT<SpMat> chol(A);
	Eigen::VectorXd u(nMatrixSystemOrder);
	u = chol.solve(b);

	Eigen::MatrixXd dPoissonSolution = m_dBoundaryConditions;
	for(i = 0; i<nMatrixSystemOrder; i++){
		dPoissonSolution(dAuxiliaryIndexes(i,0),dAuxiliaryIndexes(i,1)) = u(i);
	}
	m_dDirichletSolution = dPoissonSolution;
	CORRECT_CORNERS_DIRICHLET
}

void Poisson::PoissonDirichletNoSparseSolver(){
	bool STOP=false;
	Eigen::MatrixXd dPoissonNoSparse = m_dBoundaryConditions;
	Eigen::MatrixXd dPoissonNoSparseOld = m_dBoundaryConditions;
	Eigen::MatrixXd ERROR = Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder);
	double dDeltaX = 1.0/(m_nMatrixOrder-1);
	do{
		for(int i = 1; i<m_nMatrixOrder-1; i++)
		  for(int j = 1; j<m_nMatrixOrder-1; j++)
		  dPoissonNoSparse(i,j) = -0.25*dDeltaX*dDeltaX*m_dNonHomogeneity(i,j) + 0.25*(dPoissonNoSparse(i-1,j)+dPoissonNoSparse(i+1,j)+dPoissonNoSparse(i,j-1)+dPoissonNoSparse(i,j+1));
		  ERROR = dPoissonNoSparse-dPoissonNoSparseOld;
		  if(ERROR.norm() < 0.0000001) STOP = true;
		  dPoissonNoSparseOld = dPoissonNoSparse;
	}while(!STOP);

	//Quinas
	m_dDirichletSolution = dPoissonNoSparse;
	CORRECT_CORNERS_DIRICHLET
}

void Poisson::PoissonNeumannSparseSolver()
{
	std::vector<T> tripletList; 	//triplet, needed to fill sparse matrix
	int nMatrixSystemOrder = m_nMatrixOrder*m_nMatrixOrder;
	tripletList.reserve(PENTADIAGONAL*nMatrixSystemOrder);
	double dDeltaX = 1.0/(m_nMatrixOrder-1); 			/*Mesh has domain 0 < x < 1 and 0 < y < 1*/
	Eigen::MatrixXd dAuxiliaryIndexes = Eigen::MatrixXd::Zero(nMatrixSystemOrder,2);;
	Eigen::VectorXd b(nMatrixSystemOrder);

	/*Starting to build A*/
	long unsigned i,j;

	//Construir diagonal principal
	for(i = 0; i<nMatrixSystemOrder; i++)
		tripletList.push_back(T(i,i,-4));

	// Construir diagonal superior 1
	for(i = 0; i<= nMatrixSystemOrder-2; i++)
		if( (i+1) % (m_nMatrixOrder-2) != 0) tripletList.push_back(T(i,i+1,1));

	// Construir diagonal inferior 1
	for(i = 0; i<= nMatrixSystemOrder-2; i++)
	    if( (i+1) % (m_nMatrixOrder-2) != 0) tripletList.push_back(T(i+1,i,1));

	//Construir diagonal superior 2
	for(i = 0; i < (nMatrixSystemOrder-(m_nMatrixOrder-2)); i++)
		tripletList.push_back(T(i,i+m_nMatrixOrder-2,1));

	//Construir diagonal inferior 2
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
	m_dNeumannSolution = dPoissonSolution;
}

void Poisson::PoissonNeumannNoSparseSolver(){
	bool STOP=false;
	Eigen::MatrixXd dPoissonNoSparse = m_dBoundaryConditions;
	Eigen::MatrixXd dPoissonNoSparseOld = m_dBoundaryConditions;
	Eigen::MatrixXd ERROR = Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder);
	int k=0;
	double dDeltaX = 1.0/(m_nMatrixOrder-1);
	do{
		for(int i = 0; i<m_nMatrixOrder; i++)
		  for(int j = 0; j<m_nMatrixOrder; j++){
		   if(INTERNAL_POINT){
			   POISSON_NOSPARSE_INTERNAL_POINTS
		   }
		   else if(WEST_POINT) {
			   POISSON_NOSPARSE_WEST_POINTS
		   }
		   else if(EAST_POINT) {
			   POISSON_NOSPARSE_EAST_POINTS
		   }
		   else if(NORTH_POINT) {
			   POISSON_NOSPARSE_NORTH_POINTS
		   }
		   else if(SOUTH_POINT) {
			   POISSON_NOSPARSE_SOUTH_POINTS
		   }
		  }
		  ERROR = dPoissonNoSparse-dPoissonNoSparseOld;
		  double dERROR = ERROR.norm();
		  std::cout << "Erro: " << dERROR << " " << k <<std::endl;
		  k = k + 1;
		  std::cout << std::endl << dPoissonNoSparse << std::endl;
		  if(dERROR < 1e-15) STOP = true;
		  dPoissonNoSparseOld = dPoissonNoSparse;
	}while(!STOP);

	//Quinas
	m_dNeumannSolution = dPoissonNoSparse;
	CORRECT_CORNERS_NEUMANN
}
void Poisson::PoissonSolver(){
	if(m_bSparseOrNot){
		if(m_bDirichletOrNeumann==DIRICHLET) PoissonDirichletSparseSolver();
		if(m_bDirichletOrNeumann == NEUMANN) PoissonNeumannSparseSolver();
	} else {
		if(m_bDirichletOrNeumann==DIRICHLET) PoissonDirichletNoSparseSolver();
		if(m_bDirichletOrNeumann == NEUMANN) PoissonNeumannNoSparseSolver();
	}
	m_bCheckIfSolved = true;
}
#if WRAP_PYTHON
void Poisson::saveSolution(PyObject* pyArraySolution){
	Eigen::Map<Eigen::MatrixXd> _pyArraySolution((double *) PyArray_DATA(pyArraySolution),m_nMatrixOrder,m_nMatrixOrder);
	if(!m_bDirichletOrNeumann) _pyArraySolution = m_dDirichletSolution;
	else if(m_bDirichletOrNeumann) _pyArraySolution = m_dNeumannSolution;
}
#else
Eigen::MatrixXd Poisson::ReturnMatrix(){
	if(!m_bDirichletOrNeumann) return(m_dDirichletSolution);
	else if(m_bDirichletOrNeumann) return(m_dNeumannSolution);

	return(Eigen::MatrixXd::Zero(m_nMatrixOrder,m_nMatrixOrder));
}
#endif
#if WRAP_PYTHON
int Poisson::PoissonPython(PyObject* dBoundaryConditions,
		  	  	  	  	   PyObject* dNonHomogeneity,
		  	  	  	  	   bool bDirichletOrNeumann,
		  	  	  	  	   bool bSparseOrNot,
		  	  	  	  	   const int nMatrixOrder
		  	  	  	  	  )
{
		Eigen::Map<Eigen::MatrixXd> _dBoundaryConditions((double *) PyArray_DATA(dBoundaryConditions),nMatrixOrder,nMatrixOrder);
		Eigen::Map<Eigen::MatrixXd> _dNonHomogeneity((double *) PyArray_DATA(dNonHomogeneity),nMatrixOrder,nMatrixOrder);
        return PoissonInit(_dBoundaryConditions, _dNonHomogeneity, bDirichletOrNeumann, bSparseOrNot);
}
#endif

#if WRAP_PYTHON
BOOST_PYTHON_MODULE(libpoisson)
{
	class_<Poisson>("Poisson")
		.def("conditions", &Poisson::PoissonPython,
				(arg("dBoundaryConditions"),
				 arg("dNonHomogeneity"),
				 arg("bDirichletOrNeumann"),
				 arg("bSparseOrNot"),
				 arg("nMatrixOrder")),
				 "This function takes important values required to solve the Poisson problem on a square.")
		.def("poissonSolver", &Poisson::PoissonSolver, "Solves the Poisson problem on a square, with previous given conditions.")
		.def("saveSolution", &Poisson::saveSolution, (arg("pyArraySolution")), "Saves the solution of Poisson problem in pyArraySolution object.")
	;
}
#endif
