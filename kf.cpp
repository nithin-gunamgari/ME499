#include <math.h>
#include <fstream>
#include <assert.h>
#include <armadillo>  // using armadillo for matrix operations

using namespace std;
using namespace arma; 

int main(){
    mat A(1,1), B(1,1), H(1,1), Q(1,1), R(1,1);
    colvec u(1),x0(1);
    colvec v(1),w(1); // noise 
    colvec x(1),z(1); 
    colvec x_m(1),x_p(1),z_m(1);
    mat P_p(1,1),P_m(1,1);
    int n = A.n_cols;
    int m = H.n_rows;
    A << 0;
    B << 1;
    Q << 2*2;
    H << 1;
    R << 2*2;
    u << 12.0;
    x =  x.zeros();
    P_m = P_m.eye();
    x_m = x_m.zeros();
    x = x0;
    x_m = x0;

    ofstream myfile;
    myfile.open ("kalman.csv");

    for (int k = 0; k < 100; k++) {
	v.randn(1);
	w.randn(1);
	v = sqrt(Q) * v;
	w = sqrt(R) * w;
	x = A * x + B * u + v;
	z = H * x + w;

	x_p = A * x_m + B * u;
	P_p_ = A * P_m * trans(A) + Q;

	mat K = P_p * trans(H) * inv(H * P_p * trans(H) + R);
	x_m = x_p + K * (z - H * x_p);
	P_m = P_p - K * H * P_p;

	z_m = H * x_m;
	cout << x_m;

	myfile << x_m<<x<<z_m;
    }

    myfile.close();
    return 0;
}

