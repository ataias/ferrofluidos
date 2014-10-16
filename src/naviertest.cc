#include <FFH.hh>
#include <navierstokes.hh>
#include <chrono>

using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::chrono::steady_clock;

double f(int i, int n);

int main(){
	const int n = 6;
	double dt = 0.002;
	double mu = 1;
	double rho = 1;
	MatrixXd u(n,n); 
	u << 999, 999, 999, 999, 999, 999,
	     0, 0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0, 0;
//	cout << "u(0,0)=" << u(0,0) << endl;

	MatrixXd v(n,n);
	v << 999, 0, 0, 0, 0, 0,
	    999, 0, 0, 0, 0, 0,
	    999, 0, 0, 0, 0, 0,
	    999, 0, 0, 0, 0, 0,
	    999, 0, 0, 0, 0, 0,
	    999, 0, 0, 0, 0, 0;
//	cout << "v(2,0)=" << v(2,0) << endl;
//	cout << "v(5,1)=" << v(5,1) << endl;

	MatrixXd u_old;
	u_old = u;

	MatrixXd v_old;	
	v_old = v;

	MatrixXd fx(n,n);
	fx << 	 0, 0, 0, 0, 0, 0,
		     0, 0, 0, 0, 0, 0,
		     0, 0, 0, 0, 0, 0,
		     0, 0, 0, 0, 0, 0,
		     0, 0, 0, 0, 0, 0,
		     0, 0, 0, 0, 0, 0;

	MatrixXd fy(n,n);
	fy = fx;

	VectorXd uB(n);
	uB << 10, 10, 10, 10, 10, 10;

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

	//Deallocating memory
	u.resize(0,0); v.resize(0,0);
	u_old.resize(0,0); v_old.resize(0,0);
	fx.resize(0,0); fy.resize(0,0);
	return(0);
}

double f(int i, int n){
	return (i-1)/n;
}