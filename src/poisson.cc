#include <FFH.hh>
#include <sparsePD.hh>

int main(){
	MatrixXd S;
	S = exPoissonLDLT();
	cout << S << endl << endl;

	MatrixXd F(5,5);
	F << 0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0;

	cout << F << endl;
	return 0;
}