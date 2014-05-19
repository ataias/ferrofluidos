#include <FFH.hh>
#include <vector>
/* Routines to solve Poisson problem in square (1,1) with n x n points with Dirichlet boundary conditions
if you don't want square (1,1), you need to change dx, instead of 1/(n-1) you multiply dx*y then you solve in (y,y)
TODO: if performance with those codes is critical, make them void instead of returning! Then add one extra argument that will be changed
*/
/*-----------------------------------------------------------------------------------*/

SpMat genSMPD(int n){
	//Generate Sparse Matrix Poisson Dirichlet
	int ns = (n-2)*(n-2); int i;
	vector<T> tripletList;
	tripletList.reserve(5*ns);

	for(i = 0; i<ns; i++) //Diagonal principal
		tripletList.push_back(T(i,i,-4));
	
	for(i = 0; i< ns-1; i++){
		if( (i+1) % (n-2) != 0) tripletList.push_back(T(i,i+1,1)); // Diagonal superior 1
		if( (i+1) % (n-2) != 0) tripletList.push_back(T(i+1,i,1)); // Diagonal inferior 1
	}

	for(i = 0; i < (ns-(n-2)); i++){
		tripletList.push_back(T(i,i+n-2,1)); //Diagonal superior 2
		tripletList.push_back(T(i+n-2,i,1)); //Diagonal inferior 2
	} 

	SpMat A(ns,ns);
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	return A;
}
/*-----------------------------------------------------------------------------------*/

MatrixXd gen_index(int n){
	int ns = (n-2)*(n-2); int i,j;
	MatrixXd I = MatrixXd::Zero(ns,2);

	for(i = 1; i< (n-1); i++)
	    for(j=1; j< (n-1); j++){
	    	/*I has the indexes of the internal points of matrix of order n x n*/
	    	I(j+(i-1)*(n-2)-1,0)= j;
	    	I(j+(i-1)*(n-2)-1,1)= i;
	    }
	return I;
}
/*-----------------------------------------------------------------------------------*/

VectorXd gen_b(MatrixXd B, MatrixXd F, int n){
	int ns = (n-2)*(n-2); VectorXd b(ns); int i;
	double dx = 1.0/(n-1); double dx2 = dx*dx; double aux = 0;
	MatrixXd I = gen_index(n);

	for(i = 0; i<ns; i++){
		aux  = (I(i,0) == 1)   ? B(0,I(i,1))   : 0.0;
		aux += (I(i,0)+2 == n) ? B(n-1,I(i,1)) : 0.0;
		aux += (I(i,1) == 1)   ? B(I(i,0),0)   : 0.0;
		aux += (I(i,1)+2 == n) ? B(I(i,0),n-1) : 0.0;
    b(i) = dx2*F(I(i,0),I(i,1))-aux;
	}
	return b;
}
/*-----------------------------------------------------------------------------------*/
MatrixXd vectorToMatrix(MatrixXd B, VectorXd u, int n){
	MatrixXd S; S = B; int i;
	int ns = (n-2)*(n-2);
	MatrixXd I = gen_index(n);
	for(i = 0; i<ns; i++){
		S(I(i,0),I(i,1)) = u(i);
	}
	return S;
}