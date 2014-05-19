#include <FFH.hh>
#include <sparsePD.hh>

int main(){
	MatrixXd S;
	S = exPoissonLDLT();
	cout << S << endl;
	return 0;
}