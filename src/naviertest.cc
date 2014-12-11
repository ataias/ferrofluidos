#include <FFH.hh>
#include <navierstokes.hh>
#include <chrono>
#include <limits>

using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::chrono::steady_clock;

double f(double i, double n);
void initialize(int n, MatrixXd *u, MatrixXd *v, 
					   MatrixXd *u_old, MatrixXd *v_old, 
					   MatrixXd *fx, MatrixXd *fy);
bool isdXok(double Re, double n);
double getDt(double n, double divFactor, double Re);

int main(){
	const int n = 10; //tamanho da malha escalonada, mas dx é obtido pelo tamanho da malha normal (n-1)
	double mu = 1;
	const double rho = 1; //deixe constante
	double Re = 1/mu;
	double divFactor = 5; //must be greater than 1
	MatrixXd u, v, u_old, v_old, fx, fy;
	double dt = getDt(n, divFactor, Re);
	initialize(n,&u,&v, &u_old, &v_old, &fx, &fy);

	if(!isdXok(1/mu, n)) {
		cout << "There is a problem with your dx. Increase n.\n";
		return 1;
	}


	VectorXd uB(n);
/*
	//Let's measure the time!
	steady_clock::time_point start = steady_clock::now();
	navier_stokes_step1(n, dt, mu, rho, &u, &v, &u_old, &v_old, &fx, &fy, uB);
	steady_clock::time_point end = steady_clock::now();

	std::cout << "Computing took "
          // duration_cast is required to avoid accidentally losing precision.
          << duration_cast<microseconds>(end - start).count()
          << "us.\n";
          
	cout << "\n u_old" << endl << u_old << endl;
	cout << "\n v_old" << endl << v_old << endl;

	cout << "\n u_new" << endl << u << endl;
	cout << "\n v_new" << endl << v << endl;
*/
	//Deallocating memory
	u.resize(0,0); v.resize(0,0);
	u_old.resize(0,0); v_old.resize(0,0);
	fx.resize(0,0); fy.resize(0,0);
	return(0);
}

double f(double i, double n){
	return (i-1)/(n-1);
}

void initialize(int n, MatrixXd *u, MatrixXd *v, 
					   MatrixXd *u_old, MatrixXd *v_old, 
					   MatrixXd *fx, MatrixXd *fy){
	//Inicializar matriz u --------------------------------------------------
	*u = MatrixXd::Zero(n,n);
	for(int i=1; i<n-1; i++) (*u)(i,n-1) = 2*f(i,n)-(*u)(i,n-2); // em cima
	for(int i=1; i<n-1; i++) (*u)(i,0) = -(*u)(i,1); //embaixo
	for(int j=0; j<n; j++) (*u)(0,j) = std::numeric_limits<double>::quiet_NaN();
	(*u)(n-1,n-1) = std::numeric_limits<double>::quiet_NaN();
	(*u)(n-1,0) = std::numeric_limits<double>::quiet_NaN();
	//------------------------------------------------------------------------

	//Inicializar matriz v ---------------------------------------------------
	*v = MatrixXd::Zero(n,n);
	for(int j=1; j<n-1; j++) (*v)(0,j) = -(*v)(0,j); //esquerda
	for(int j=1; j<n-1; j++) (*v)(n-1,j) = -(*v)(n-2,j); //direita
	for(int i=0; i<n; i++) (*v)(i,0) = std::numeric_limits<double>::quiet_NaN();
	(*v)(n-1,n-1) = std::numeric_limits<double>::quiet_NaN();
	(*v)(0,n-1) = std::numeric_limits<double>::quiet_NaN();
	//------------------------------------------------------------------------

	(*u_old) = (*u);	
	(*v_old) = (*v);
	*fx = MatrixXd::Zero(n,n);
	*fy = MatrixXd::Zero(n,n);


	//Imprimir matrizes
	cout << "\n u" << endl << *u << endl;
	cout << "\n v" << endl << *v << endl;
}

bool isdXok(double Re, double n){
	double dx = 1/(n-2); //because n = n' + 1
	if(dx < (1/Re)) return true;
	else return false;

}

double getDt(double n, double divFactor, double Re){
	double dx = 1/(n-2); //n é o tamanho da malha escalonada, por isso está n-2 ao invés de n-1
	double dt1 = 0.25*dx*dx/Re;
	double dt2 = dx;
	double dt;
	if(dt1 < dt2) dt = dt1/divFactor;
	else dt = dt2/divFactor;
	cout << "dt = " << dt << std::endl;
	return dt;
}