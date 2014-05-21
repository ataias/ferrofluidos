#include <FFH.hh>
#include <navierstokes.hh>
#include <chrono>

using std::chrono::duration_cast;
using std::chrono::microseconds;
using std::chrono::steady_clock;

int main(){
	const int n = 5;
	double dt = 0.2;
	double mu = 1;
	double rho = 1;
	MatrixXd u(n,n); 
	u << 1, 1, 1, 1, 1,
	     0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0,
	     0, 0, 0, 0, 0;

	MatrixXd v;
	v = u;

	MatrixXd u_old;
	u_old = u;

	MatrixXd v_old;
	v_old = v;

	MatrixXd fx(n,n);
	fx << 	0, 0, 0, 0, 0,
	     	0, 0, 0, 0, 0,
	     	0, 0, 0, 0, 0,
	     	0, 0, 0, 0, 0,
	     	0, 0, 0, 0, 0;

	MatrixXd fy(n,n);
	fy = fx;

	//Let's measure the time!
	steady_clock::time_point start = steady_clock::now();
	navier_stokes_step1(n, dt, mu, rho, &u, &v, &u_old, &v_old, &fx, &fy);
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