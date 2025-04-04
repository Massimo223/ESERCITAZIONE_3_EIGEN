#include <iostream>
#include <iomanip>
#include <vector>
#include "Eigen/Eigen"

using namespace Eigen;
using namespace std;

Vector2d PALU(Matrix2d A, Vector2d b) {
    
	PartialPivLU<Matrix2d> lu(A);
	Vector2d x = A.partialPivLu().solve(b);          // Funzione che, dati A e b, restituisce la soluzione del sistema Ax = b
	return x;                                        // applicando la fattorizzazione PA = LU ad A.
}

Vector2d QR(Matrix2d A, Vector2d b) {
	
	Vector2d x = A.householderQr().solve(b);         // Funzione che, dati A e b, restituisce la soluzione del sistema Ax = b
	return x;                                        // applicando la fattorizzazione A = QR ad A.
}


int main() {
    
	Matrix2d A_1;
	A_1 << 5.547001962252291e-01, -3.770900990025203e-02, 8.320502943378437e-01, -9.992887623566787e-01;
	Vector2d b_1;
    b_1 << -5.169911863249772e-01, 1.672384680188350e-01;
	Matrix2d A_2; 
	A_2 << 5.547001962252291e-01, -5.540607316466765e-01, 8.320502943378437e-01, -8.324762492991313e-01;
	Vector2d b_2;
    b_2	<< -6.394645785530173e-04, 4.259549612877223e-04;
	Matrix2d A_3;
	A_3 << 5.547001962252291e-01, -5.547001955851905e-01, 8.320502943378437e-01, -8.320502947645361e-01;
	Vector2d b_3;
    b_3 << -6.400391328043042e-10, 4.266924591433963e-10;
	Vector2d x_sol;
	x_sol << -1.0e+00, -1.0e+00; 
	vector<Vector2d> x_palu;
	vector<Vector2d> x_qr;
	x_palu.push_back(PALU(A_1,b_1));
	x_palu.push_back(PALU(A_2,b_2));
	x_palu.push_back(PALU(A_3,b_3));
	x_qr.push_back(QR(A_1,b_1));
	x_qr.push_back(QR(A_2,b_2));
	x_qr.push_back(QR(A_3,b_3));
	int n_di_sistemi = x_palu.size();
	
	for (int j = 0; j < n_di_sistemi; j++) {
		double errore_relativo = ((x_sol - x_palu[j]).norm()) / (x_sol.norm());
		cout << "L'errore relativo del sistema " << j + 1 << " risolto con la fattorizzazione PALU è: " << scientific << setprecision(16) << errore_relativo << "." << endl;
	}
	
	for (int j = 0; j < n_di_sistemi; j++) {
		double errore_relativo = ((x_sol - x_qr[j]).norm()) / (x_sol.norm());
		cout << "L'errore relativo del sistema " << j + 1 << " risolto con la fattorizzazione QR è: " << scientific << setprecision(16) << errore_relativo << "." << endl;
	}
		
}
