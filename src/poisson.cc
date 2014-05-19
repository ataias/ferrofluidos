#include <FFH.hh>
#include <sparsePD.hh>

int main(){
	int n = 5; int ns = (n-2)*(n-2);
	SpMat A; A = genSMPD(n);

	MatrixXd B(5,5); 
	B << 1, 1, 1, 1, 1,
	     0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0;

	MatrixXd F(5,5); 
	F << 0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0;

	VectorXd b; b = gen_b(B,F,n);


	SimplicialLDLT<SpMat> chol(A);
	VectorXd u(ns);
	u = chol.solve(b);

	MatrixXd S = vectorToMatrix(B, u, n);
	cout << S << endl;
	return 0;
}