#include <FFH.hh>
#include <sparsePD.hh>

void hey(MatrixXd *A){
	(*A)(0,0)=4;
}

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

	hey(&F);
	cout << F << endl;
	return 0;
}